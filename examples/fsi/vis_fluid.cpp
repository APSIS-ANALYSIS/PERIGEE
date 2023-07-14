// ============================================================================
// vis_fluid.cpp
//
// This routine provides visualization of the fluid sub-domain in FSI problems.
//
// Date: Jan 17 2022
// ============================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "APart_Node_FSI.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "FEAElement_Tet4.hpp"
#include "VisDataPrep_ALE_NS.hpp"
#include "VTK_Writer_FSI_Tet4.hpp"

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
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

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
  if(isXML) PetscPrintf(PETSC_COMM_WORLD, "-xml: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-xml: false \n");

  if(isClean) PetscPrintf(PETSC_COMM_WORLD, "-clean: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-clean: false \n");

  // If demand cleaning, remove all previous visualization files
  if( isClean )
  {
    SYS_T::execute("rm -rf VIS_F_*_p*.vtu");
    SYS_T::execute("rm -rf VIS_F_*.pvtu");
    SYS_T::execute("rm -rf VIS_F_.pvd");
  }
  
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_v_file, 0);

  SYS_T::print_fatal_if(size != PartBasic->get_cpu_size(), "Error: number of processors does not match with prepost! \n");

  SYS_T::commPrint("===> %d processor(s) are assigned.", size);

  FEANode * fNode = new FEANode(part_v_file, rank);

  ALocal_IEN * locIEN_v = new ALocal_IEN(part_v_file, rank);
  ALocal_IEN * locIEN_p = new ALocal_IEN(part_p_file, rank);

  IAGlobal_Mesh_Info * GMIptr_v = new AGlobal_Mesh_Info_FEM_3D(part_v_file, rank);
  IAGlobal_Mesh_Info * GMIptr_p = new AGlobal_Mesh_Info_FEM_3D(part_p_file, rank);

  ALocal_Elem * locElem = new ALocal_Elem(part_v_file, rank);

  APart_Node * pNode_v = new APart_Node(part_v_file, rank);
  APart_Node * pNode_p = new APart_Node(part_p_file, rank);

  IQuadPts * quad = new QuadPts_vis_tet4();

  quad -> print_info();

  FEAElement * element = new FEAElement_Tet4( quad-> get_num_quadPts() );

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
  IVisDataPrep * visprep = new VisDataPrep_ALE_NS();
  visprep->print_info();

  double ** pointArrays = new double * [3];
  pointArrays[0] = new double [pNode_v->get_nlocghonode() * 3];
  pointArrays[1] = new double [pNode_p->get_nlocghonode() * 1];
  pointArrays[2] = new double [pNode_v->get_nlocghonode() * 3];

  VTK_Writer_FSI_Tet4 * vtk_w = new VTK_Writer_FSI_Tet4(
      GMIptr_v->get_nElem(), element_part_file );

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
        pNode_v, pNode_p, GMIptr_v->get_nFunc(), GMIptr_p->get_nFunc(),
        pointArrays);

    vtk_w->writeOutput_fluid( fNode, locIEN_v, locIEN_p, fIEN, locElem,
        visprep, element, quad, pointArrays, rank, size,
        num_subdomain_nodes,
        time * dt, out_bname, name_to_write, isXML );
  }


  MPI_Barrier(PETSC_COMM_WORLD);

  // Clean up memory
  delete quad; delete element; delete visprep;
  delete fNode; delete locIEN_v; delete locIEN_p; delete GMIptr_v; delete GMIptr_p;
  delete PartBasic; delete locElem; delete pNode_v; delete pNode_p;
  delete [] pointArrays[0]; delete [] pointArrays[1]; delete [] pointArrays[2];
  delete [] pointArrays; delete vtk_w;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF 
