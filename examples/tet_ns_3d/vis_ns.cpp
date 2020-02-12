// ==================================================================
// vis_ns.cpp
//
// Date Created: Feb. 12 2020
// ==================================================================
#include "Sys_Tools.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "APart_Node.hpp"
#include "FEANode.hpp"
#include "ALocal_IEN.hpp"

int main( int argc, char * argv[] )
{
  const std::string element_part_file = "epart.h5";
  const std::string anode_mapping_file = "node_mapping.h5";
  const std::string pnode_mapping_file = "post_node_mapping.h5";
  
  std::string part_file="postpart";
  
  std::string sol_bname("SOL_");

  std::string out_bname = sol_bname;

  const int dof = 4;

  int time_start = 0;
  int time_step = 1;
  int time_end = 1;
  double dt = 0.1;
  bool isXML = true;
  bool isRestart = false;

  PetscMPIInt rank, size;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionReal("-dt", dt);
  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionString("-out_bname", out_bname);
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionBool("-xml", isXML);
  SYS_T::GetOptionBool("-restart", isRestart);

  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-sol_bname:", sol_bname);
  SYS_T::cmdPrint("-out_bname:", out_bname);
  SYS_T::cmdPrint("-time_start:", time_start);
  SYS_T::cmdPrint("-time_step:", time_step);
  SYS_T::cmdPrint("-time_end:", time_end);
  SYS_T::cmdPrint("-dt:",dt);
  if(isXML) PetscPrintf(PETSC_COMM_WORLD, "-xml: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-xml: false \n");

  if(isRestart) PetscPrintf(PETSC_COMM_WORLD, "-restart: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-restart: false \n");
  
  // Clean the visualization files if not restart
  if( !isRestart )
  {
    int sysret = system("rm -rf *_p*.vtu");
    SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
    sysret = system("rm -rf *.pvtu");
    SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
    sysret = system("rm -rf *_.pvd");
    SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  }
  
  SYS_T::commPrint("===> Reading mesh files ... ");
  
  FEANode * fNode = new FEANode(part_file, rank);
  
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);
  
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);
  
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
 
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);
  
  APart_Node * pNode = new APart_Node(part_file, rank);
  
  SYS_T::commPrint("Done! \n");
  
  if(size != PartBasic->get_cpu_size()) SYS_T::print_fatal(
      "Error: number of processors does not match with prepost! \n");

  PetscPrintf(PETSC_COMM_WORLD,
      "\n===> %d processor(s) are assigned for:", size);
  PetscPrintf(PETSC_COMM_WORLD, "Postprocessing - visualization.\n");








  // ===== Clean the memory =====
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic; delete locElem;
  delete pNode;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
