// ============================================================================
// wall_solver.cpp
//
// Triangle element based finite element code for wall mechanics.
//
// ============================================================================
#include "HDF5_Reader.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_EBC_wall.hpp"
#include "FEAElement_Triangle3_membrane.hpp"
#include "FEAElement_Triangle6_membrane.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "PLocAssem_Tet_CMM_GenAlpha.hpp"
#include "PGAssem_Tet_CMM_GenAlpha.hpp"

int main( int argc, char *argv[] )
{
  // Prestress solver parameters
  const bool   prestress_flag = false;

  // Generalized-alpha rho_inf
  double genA_rho_inf = 0.5;
  bool is_backward_Euler = false;

  // We assume that a 3D solver has been called (to generate the wall traction)
  // and a suite of command line arguments has been saved to disk
  hid_t solver_cmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  
  HDF5_Reader * cmd_h5r = new HDF5_Reader( solver_cmd_file );

  const double wall_density = cmd_h5r -> read_doubleScalar("/", "wall_density");
  const double wall_poisson = cmd_h5r -> read_doubleScalar("/", "wall_poisson");
  const double wall_kappa   = cmd_h5r -> read_doubleScalar("/", "wall_kappa");
  const int nqp_tri = cmd_h5r -> read_intScalar("/", "nqp_tri");
  const int nqp_tet = cmd_h5r -> read_intScalar("/", "nqp_tet");
  const double fl_density = cmd_h5r -> read_doubleScalar("/", "fl_density");
  const double fl_mu = cmd_h5r -> read_doubleScalar("/", "fl_mu");

  delete cmd_h5r; H5Fclose(solver_cmd_file);

  // Partition filename prefix
  std::string part_file("part");

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  
  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();
  
  // ===== Read Command Line Arguments =====
  SYS_T::GetOptionString("-part_file",       part_file);
  SYS_T::GetOptionReal(  "-rho_inf",         genA_rho_inf);
  SYS_T::GetOptionBool(  "-is_backward_Euler", is_backward_Euler);

  // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint(      "-part_file:",       part_file);
  if( is_backward_Euler )
    SYS_T::commPrint(     "-is_backward_Euler: true \n");
  else
    SYS_T::cmdPrint(      "-rho_inf:",         genA_rho_inf);

  // ===== Load Analysis Data Structure =====
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file);

  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(), "Error: Assigned CPU number does not match the partition. \n");

  delete PartBasic;

   // Global mesh info
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file, rank);

  // Local sub-domain's element indices
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  // Local sub-domain's nodal indices
  APart_Node * pNode = new APart_Node(part_file, rank);

  // ===== Quadrature rules =====
  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );

  // ===== Finite element containers =====
  SYS_T::commPrint("===> Set up volumetric and surface element containers. \n");
  FEAElement * elementw = nullptr;

  if( GMIptr->get_elemType() == 501 )          // linear tet
    elementw = new FEAElement_Triangle3_membrane( nqp_tri );
  else if( GMIptr->get_elemType() == 502 )     // quadratic tet
    elementw = new FEAElement_Triangle6_membrane( nqp_tri );
  else SYS_T::print_fatal("Error: Element type not supported.\n");

  // Local sub-domain's wall elemental (Neumann) BC for CMM
  ALocal_EBC * locebc_wall = new ALocal_EBC_wall(part_file, rank, quads, "ebc_wall", prestress_flag);

  // ===== Generalized-alpha =====
  SYS_T::commPrint("===> Set up the generalized-alpha time integration scheme.\n");
  TimeMethod_GenAlpha * tm_galpha_ptr = nullptr;

  if( is_backward_Euler )
    tm_galpha_ptr = new TimeMethod_GenAlpha( 1.0, 1.0, 1.0 );
  else
    tm_galpha_ptr = new TimeMethod_GenAlpha( genA_rho_inf, false );

  tm_galpha_ptr->print_info();

  // ===== Local Assembly Routine =====
  const double bs_beta = 0.0;
  const double c_tauc = 0.0;
  IPLocAssem * locAssem_ptr = new PLocAssem_Tet_CMM_GenAlpha(
      tm_galpha_ptr, nqp_tet, fl_density, fl_mu, bs_beta,
      wall_density, wall_poisson, wall_kappa,
      c_tauc, GMIptr->get_elemType() );

  // ===== Solution vector =====
  PDNSolution * base = new PDNSolution_NS( pNode, 0 );

  PDNSolution * sol = new PDNSolution_NS( pNode, 0 );

  PDNSolution * dot_sol = new PDNSolution_NS( pNode, 0 );

  PDNSolution * sol_wall_disp = new PDNSolution_Wall_Disp( pNode, 0 );

  PDNSolution * dot_sol_wall_disp = new PDNSolution_Wall_Disp( pNode, 0 );

  int    initial_index = 0;          // restart solution time index
  double initial_time = 0.0;         // restart time
  double initial_step = 1.0e-3;      // restart simulation time step size
  std::string restart_name = "SOL_re"; // restart solution base name

  // Read in [pres velo]
  SYS_T::file_check(restart_name.c_str());
  sol->ReadBinary(restart_name.c_str());

  // Read in [dot_pres dot_velo]
  std::string restart_dot_name = "dot_" + restart_name;
  SYS_T::file_check(restart_dot_name.c_str());
  dot_sol->ReadBinary(restart_dot_name.c_str());

  // ===== Time step info =====
  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);







  delete locElem;
  delete GMIptr; delete quads; delete elementw; delete locebc_wall; delete tm_galpha_ptr;
  delete pNode; delete locAssem_ptr; delete base; delete sol; delete dot_sol;
  delete sol_wall_disp; delete dot_sol_wall_disp; delete timeinfo;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
