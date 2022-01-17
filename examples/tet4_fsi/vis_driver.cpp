// ============================================================================
// vis_driver.cpp
// This is the visualization driver for the whole FSI system.
//
// Author: Ju Liu
// Date: Jan 17 2021
// ============================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "FEAElement_Tet4.hpp"
#include "VisDataPrep_FSI.hpp"

#include "APart_Node.hpp"

int main( int argc, char * argv[] )
{
  const std::string element_part_file = "epart.h5";
  const std::string anode_mapping_file = "node_mapping.h5";
  const std::string pnode_mapping_file = "post_node_mapping.h5";
  const std::string part_v_file="./ppart/postpart_v";
  const std::string part_p_file="./ppart/postpart_p";

  std::string sol_bname("SOL_");
  std::string out_bname("VIS_FSI_");
  
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
  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionString("-out_bname", out_bname);
  SYS_T::GetOptionBool("-xml", isXML);
  SYS_T::GetOptionBool("-clean", isClean);

  // Correct time_step if it does not match with sol_rec_freq
  if( time_step % sol_rec_freq != 0 ) time_step = sol_rec_freq;

  SYS_T::cmdPrint("-sol_bname:", sol_bname);
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
    SYS_T::execute("rm -rf VIS_FSI_*_p*.vtu");
    SYS_T::execute("rm -rf VIS_FSI_*.pvtu");
    SYS_T::execute("rm -rf VIS_FSI_.pvd");
  }

  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_v_file, 0);
  
  SYS_T::print_fatal_if(size != PartBasic->get_cpu_size(), "Error: number of processors does not match with prepost! \n");
  
  SYS_T::commPrint("===> %d processor(s) are assigned.", size);

  FEANode * fNode = new FEANode(part_v_file, rank);

  ALocal_IEN * locIEN_v = new ALocal_IEN(part_v_file, rank);
  ALocal_IEN * locIEN_p = new ALocal_IEN(part_p_file, rank);
  
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_v_file, rank);
  
  ALocal_Elem * locElem = new ALocal_Elem(part_v_file, rank);
  
  APart_Node * pNode_v = new APart_Node(part_v_file, rank);
  APart_Node * pNode_p = new APart_Node(part_p_file, rank);

  IQuadPts * quad = new QuadPts_vis_tet4();

  quad -> print_info();

  FEAElement * element = new FEAElement_Tet4( quad-> get_num_quadPts() );

  IVisDataPrep * visprep = new VisDataPrep_FSI();

  visprep->print_info();











  delete quad; delete element; delete visprep;
  delete fNode; delete locIEN_v; delete locIEN_p; delete GMIptr; delete PartBasic;
  delete locElem; delete pNode_v; delete pNode_p;  
  PetscFinalize();
  return EXIT_SUCCESS;
}

// 
