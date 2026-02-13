// ============================================================================
// vis_elastodynamics.cpp
//
// Visualization driver for the elastodynamics solver.
//
// Date Created: Nov. 4 2023
// ============================================================================
#include "AGlobal_Mesh_Info.hpp"
#include "ANL_Tools.hpp"
#include "QuadPtsFactory.hpp"
#include "FEAElementFactory.hpp"
#include "VisDataPrep_Elastodynamics.hpp"
#include "VTK_Writer_Elastodynamics.hpp"

int main( int argc, char * argv[] ) 
{ 
  const std::string element_part_file = "epart.h5";
  const std::string anode_mapping_file = "node_mapping.h5";
  const std::string pnode_mapping_file = "post_node_mapping.h5";
  const std::string part_file="./ppart/part";
  constexpr int dof = 3;

  std::string sol_bname("SOL_disp_");
  std::string out_bname = sol_bname;
  int time_start = 0, time_step = 1, time_end = 1;
  bool isXML = true, isRestart = false;

  // Read analysis code parameter if the solver_cmd.h5 exists
  hid_t prepcmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

  double dt       = cmd_h5r -> read_doubleScalar("/","init_step");
  double module_E = cmd_h5r -> read_doubleScalar("/","youngs_module");
  double nu       = cmd_h5r -> read_doubleScalar("/","poissons_ratio");

  delete cmd_h5r; H5Fclose(prepcmd_file);

  // ===== Initialize the MPI run =====
#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  PetscMPIInt rank = SYS_T::get_MPI_rank();
  PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionReal("-dt", dt);
  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionString("-out_bname", out_bname);
  SYS_T::GetOptionBool("-xml", isXML);
  SYS_T::GetOptionBool("-restart", isRestart);

  SYS_T::commPrint("=== Command line arguments ===\n");
  SYS_T::cmdPrint("-sol_bname:", sol_bname);
  SYS_T::cmdPrint("-out_bname:", out_bname);
  SYS_T::cmdPrint("-time_start:", time_start);
  SYS_T::cmdPrint("-time_step:", time_step);
  SYS_T::cmdPrint("-time_end:", time_end);
  SYS_T::cmdPrint("-dt:",dt);
  if(isXML) SYS_T::commPrint("-xml: true \n");
  else SYS_T::commPrint("-xml: false \n");

  if(isRestart) SYS_T::commPrint("-restart: true \n");
  else SYS_T::commPrint("-restart: false \n");
  SYS_T::commPrint("==============================\n");
  
  // Clean the visualization files if not restart
  if( !isRestart )
  {
    SYS_T::execute("rm -rf *_p*.vtu");
    SYS_T::execute("rm -rf *.pvtu");
    SYS_T::execute("rm -rf *_.pvd");
  }

  auto fNode = SYS_T::make_unique<FEANode>(part_file, rank);

  auto locIEN = SYS_T::make_unique<ALocal_IEN>(part_file, rank);

  auto GMIptr = SYS_T::make_unique<AGlobal_Mesh_Info>(part_file,rank);

  auto locElem = SYS_T::make_unique<ALocal_Elem>(part_file, rank);

  auto pNode = SYS_T::make_unique<APart_Node>(part_file, rank);

  SYS_T::print_fatal_if(size != ANL_T::get_cpu_size(part_file, rank), "Error: number of processors does not match with prepost! \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for:", size);

  // Allocate the quadrature rule and element container
  const auto elemType = GMIptr->get_elemType();
  std::unique_ptr<IQuadPts> quad = QuadPtsFactory::createVisQuadrature(elemType);
  std::unique_ptr<FEAElement> element = ElementFactory::createVolElement(elemType,
    quad->get_num_quadPts());

  // Print the sampling points on screen
  quad -> print_info();

  // Create the visualization data object
  std::unique_ptr<IVisDataPrep> visprep = SYS_T::make_unique<VisDataPrep_Elastodynamics>();

  visprep->print_info();

  // Allocate the container to store the solution values into physical fields.
  double ** solArrays = new double * [visprep->get_ptarray_size()];
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    solArrays[ii] = new double [pNode->get_nlocghonode() * visprep->get_ptarray_comp_length(ii)];

  // VTK writer 
  auto vtk_w = SYS_T::make_unique<VTK_Writer_Elastodynamics>( GMIptr->get_nElem(), 
      GMIptr->get_nLocBas(), element_part_file, module_E, nu );
  
  std::ostringstream time_index;

  const auto anode_mapping = VIS_T::readNodeMapping(anode_mapping_file, "old_2_new");
  const auto pnode_mapping = VIS_T::readNodeMapping(pnode_mapping_file, "new_2_old");

  for(int time = time_start; time<=time_end; time+= time_step)
  {
    std::string name_to_read(sol_bname);
    std::string name_to_write(out_bname);
    time_index.str("");
    time_index<< 900000000 + time;
    name_to_read.append(time_index.str());
    name_to_write.append(time_index.str());

    SYS_T::commPrint("Time %d: Read %s and Write %s \n",
        time, name_to_read.c_str(), name_to_write.c_str() );

    visprep->get_pointArray(name_to_read, anode_mapping, pnode_mapping,
        pNode.get(), dof, solArrays);

    vtk_w->writeOutput( fNode.get(), locIEN.get(), locElem.get(),
        visprep.get(), element.get(), quad.get(), solArrays,
        rank, size, time * dt, sol_bname, out_bname, name_to_write, isXML );
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== Clean the memory =====
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    delete [] solArrays[ii];
  delete [] solArrays;

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
