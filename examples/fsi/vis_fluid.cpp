// ============================================================================
// vis_fluid.cpp
//
// This routine provides visualization of the fluid sub-domain in FSI problems.
//
// Date: Jan 17 2022
// ============================================================================
#include "AGlobal_Mesh_Info.hpp"
#include "ANL_Tools.hpp"
#include "QuadPtsFactory.hpp"
#include "FEAElementFactory.hpp"
#include "VisDataPrep_ALE_NS.hpp"
#include "VTK_Writer_FSI.hpp"

int main( int argc, char * argv[] )
{
  const std::string element_part_file = "epart.h5";
  const std::string an_v_mapping_file = "node_mapping_v.h5";
  const std::string an_p_mapping_file = "node_mapping_p.h5";
  const std::string pn_v_mapping_file = "post_node_mapping_v.h5";
  const std::string pn_p_mapping_file = "post_node_mapping_p.h5";
  const std::string part_v_file="./ppart/postpart_v";
  const std::string part_p_file="./ppart/postpart_p";

  std::string disp_sol_bname("SOL_disp_");
  std::string velo_sol_bname("SOL_velo_");
  std::string pres_sol_bname("SOL_pres_");
  std::string out_bname("VIS_F_");

  int time_start = 0;
  int time_step = 1;
  int time_end = 1;
  bool isXML = true;
  bool isClean = true;

  // Load analysis code parameter from solver_cmd.h5 file
  hid_t prepcmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

  double dt = cmd_h5r -> read_doubleScalar("/","init_step");

  const int sol_rec_freq = cmd_h5r -> read_intScalar("/", "sol_record_freq");

  delete cmd_h5r; H5Fclose(prepcmd_file);

  // ===== PETSc Initialization =====
#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionReal("-dt", dt);
  SYS_T::GetOptionString("-disp_sol_bname", disp_sol_bname);
  SYS_T::GetOptionString("-velo_sol_bname", velo_sol_bname);
  SYS_T::GetOptionString("-pres_sol_bname", pres_sol_bname);
  SYS_T::GetOptionString("-out_bname", out_bname);
  SYS_T::GetOptionBool("-xml", isXML);
  SYS_T::GetOptionBool("-clean", isClean);

  // Correct time_step if it does not match with sol_rec_freq
  if( time_step % sol_rec_freq != 0 ) time_step = sol_rec_freq;

  SYS_T::cmdPrint("-disp_sol_bname:", disp_sol_bname);
  SYS_T::cmdPrint("-velo_sol_bname:", velo_sol_bname);
  SYS_T::cmdPrint("-pres_sol_bname:", pres_sol_bname);
  SYS_T::cmdPrint("-out_bname:", out_bname);
  SYS_T::cmdPrint("-time_start:", time_start);
  SYS_T::cmdPrint("-time_step:", time_step);
  SYS_T::cmdPrint("-time_end:", time_end);
  SYS_T::cmdPrint("-dt:",dt);
  if(isXML) SYS_T::commPrint("-xml: true \n");
  else SYS_T::commPrint("-xml: false \n");

  if(isClean) SYS_T::commPrint("-clean: true \n");
  else SYS_T::commPrint("-clean: false \n");

  // If demand cleaning, remove all previous visualization files
  if( isClean )
  {
    SYS_T::execute("rm -rf VIS_F_*_p*.vtu");
    SYS_T::execute("rm -rf VIS_F_*.pvtu");
    SYS_T::execute("rm -rf VIS_F_.pvd");
  }
  
  SYS_T::print_fatal_if(size != ANL_T::get_cpu_size(part_v_file, 0), "Error: number of processors does not match with prepost! \n");

  SYS_T::commPrint("===> %d processor(s) are assigned.", size);

  auto fNode = SYS_T::make_unique<FEANode>(part_v_file, rank);

  auto locIEN_v = SYS_T::make_unique<ALocal_IEN>(part_v_file, rank);
  auto locIEN_p = SYS_T::make_unique<ALocal_IEN>(part_p_file, rank);

  auto GMIptr_v = SYS_T::make_unique<AGlobal_Mesh_Info>(part_v_file, rank);
  auto GMIptr_p = SYS_T::make_unique<AGlobal_Mesh_Info>(part_p_file, rank);

  auto locElem = SYS_T::make_unique<ALocal_Elem>(part_v_file, rank);

  auto pNode_v = SYS_T::make_unique<APart_Node>(part_v_file, rank);
  auto pNode_p = SYS_T::make_unique<APart_Node>(part_p_file, rank);

  // Allocate the quadrature rule and element container
  // We assume that the same element type is used for pressure and velocity
  const auto elemType = GMIptr_v->get_elemType();
  std::unique_ptr<IQuadPts> quad = QuadPtsFactory::createVisQuadrature(elemType);
  std::unique_ptr<FEAElement> element = ElementFactory::createVolElement(elemType,
    quad->get_num_quadPts());

  quad -> print_info();

  // For the fluid subdomain, we need to prepare a mapping from the
  // FSI nodal index to the fluid subdomain nodal index
  std::vector<int> subdomain_nodes; subdomain_nodes.clear();
  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    if( locElem -> get_elem_tag(ee) == 0 )
      VEC_T::insert_end( subdomain_nodes, locIEN_v -> get_LIEN(ee) );
  }
  VEC_T::sort_unique_resize( subdomain_nodes );
  const int num_subdomain_nodes = static_cast<int>( subdomain_nodes.size() );

  const int nelem = locElem -> get_nlocalele();
  const int nlocbas = locIEN_v -> get_stride();
  std::vector<int> fIEN( nelem * nlocbas );
  for( int ee=0; ee<nelem; ++ee )
    for(int ii=0; ii<nlocbas; ++ii)
      fIEN[ee*nlocbas+ii] = VEC_T::get_pos( subdomain_nodes, locIEN_v->get_LIEN(ee, ii) );

  // Visualization prepration
  std::unique_ptr<IVisDataPrep> visprep = SYS_T::make_unique<VisDataPrep_ALE_NS>();  
  visprep->print_info();

  double ** pointArrays = new double * [3];
  pointArrays[0] = new double [pNode_v->get_nlocghonode() * 3];
  pointArrays[1] = new double [pNode_p->get_nlocghonode() * 1];
  pointArrays[2] = new double [pNode_v->get_nlocghonode() * 3];

  auto vtk_w = SYS_T::make_unique<VTK_Writer_FSI>( GMIptr_v->get_nElem(), 
      element->get_nLocBas(), element_part_file );

  std::ostringstream time_index;

  for(int time = time_start; time<=time_end; time += time_step)
  {
    std::string disp_name_to_read(disp_sol_bname);
    std::string velo_name_to_read(velo_sol_bname);
    std::string pres_name_to_read(pres_sol_bname);
    std::string name_to_write(out_bname);
    time_index.str("");
    time_index<< 900000000 + time;
    disp_name_to_read.append(time_index.str());
    velo_name_to_read.append(time_index.str());
    pres_name_to_read.append(time_index.str());
    name_to_write.append(time_index.str());

    SYS_T::commPrint("Time %d: Read %s %s %s and Write %s \n",
        time, disp_name_to_read.c_str(), pres_name_to_read.c_str(),
        velo_name_to_read.c_str(), name_to_write.c_str() );

    visprep->get_pointArray(disp_name_to_read, pres_name_to_read,
        velo_name_to_read, an_v_mapping_file, an_p_mapping_file,
        pn_v_mapping_file, pn_p_mapping_file,
        pNode_v.get(), pNode_p.get(), GMIptr_v->get_nFunc(), GMIptr_p->get_nFunc(),
        pointArrays);

    vtk_w->writeOutput_fluid( fNode.get(), locIEN_v.get(), locIEN_p.get(), fIEN, locElem.get(),
        visprep.get(), element.get(), quad.get(), pointArrays, rank, size,
        num_subdomain_nodes,
        time * dt, out_bname, name_to_write, isXML );
  }


  MPI_Barrier(PETSC_COMM_WORLD);

  // Clean up memory
  delete [] pointArrays[0]; delete [] pointArrays[1]; delete [] pointArrays[2];
  delete [] pointArrays;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF 
