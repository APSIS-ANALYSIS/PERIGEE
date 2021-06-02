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
#include "PGAssem_Tet_Wall.hpp"
#include "PTime_CMM_Solver.hpp"

int main( int argc, char *argv[] )
{
  // Prestress solver parameters
  // Generalized-alpha rho_inf
  double genA_rho_inf = 0.5;
  bool is_backward_Euler = false;
  
  // Estimate of num nonzeros per row for the sparse tangent matrix
  int nz_estimate = 300;

  // Prestress tolerance
  double prestress_disp_tol = 1.0e-6;

  // Nonlinear solver parameters
  double nl_rtol = 1.0e-3;           // convergence criterion relative tolerance
  double nl_atol = 1.0e-6;           // convergence criterion absolute tolerance
  double nl_dtol = 1.0e3;             // divergence criterion
  int    nl_maxits = 20;             // maximum number if nonlinear iterations
  int    nl_refreq = 4;              // frequency of tangent matrix renewal
  int    nl_threshold = 4;           // threshold of tangent matrix renewal

  // Time stepping parameters
  double initial_time = 0.0;         // time of initial condition
  double initial_step = 0.1;         // time step size
  int    initial_index = 0;          // index of initial condition
  double final_time = 1.0;           // end time of simulation
  std::string sol_bName("SOL_");     // base name of the solution file
  int    ttan_renew_freq = 1;        // frequency of tangent matrix renewal
  int    sol_record_freq = 1;        // frequency for recording the solution

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
  SYS_T::GetOptionReal(  "-nl_rtol",         nl_rtol);
  SYS_T::GetOptionReal(  "-nl_atol",         nl_atol);
  SYS_T::GetOptionReal(  "-nl_dtol",         nl_dtol);
  SYS_T::GetOptionInt(   "-nl_maxits",       nl_maxits);
  SYS_T::GetOptionInt(   "-nl_refreq",       nl_refreq);
  SYS_T::GetOptionInt(   "-nl_threshold",    nl_threshold);
  SYS_T::GetOptionBool(  "-is_backward_Euler", is_backward_Euler);
  SYS_T::GetOptionReal(  "-init_time",       initial_time);
  SYS_T::GetOptionReal(  "-fina_time",       final_time);
  SYS_T::GetOptionReal(  "-init_step",       initial_step);
  SYS_T::GetOptionInt(   "-init_index",      initial_index);
  SYS_T::GetOptionInt(   "-ttan_freq",       ttan_renew_freq);
  SYS_T::GetOptionInt(   "-sol_rec_freq",    sol_record_freq);

  // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint(      "-part_file:",       part_file);
  SYS_T::cmdPrint(      "-nl_rtol:",         nl_rtol);
  SYS_T::cmdPrint(      "-nl_atol:",         nl_atol);
  SYS_T::cmdPrint(      "-nl_dtol:",         nl_dtol);
  SYS_T::cmdPrint(      "-nl_maxits:",       nl_maxits);
  SYS_T::cmdPrint(      "-nl_refreq:",       nl_refreq);
  SYS_T::cmdPrint(      "-nl_threshold:",    nl_threshold);
  
  if( is_backward_Euler )
    SYS_T::commPrint(     "-is_backward_Euler: true \n");
  else
    SYS_T::cmdPrint(      "-rho_inf:",         genA_rho_inf);
  
  SYS_T::cmdPrint(      "-init_time:",       initial_time);
  SYS_T::cmdPrint(      "-init_step:",       initial_step);
  SYS_T::cmdPrint(      "-init_index:",      initial_index);
  SYS_T::cmdPrint(      "-fina_time:",       final_time);
  SYS_T::cmdPrint(      "-ttan_freq:",       ttan_renew_freq);
  SYS_T::cmdPrint(      "-sol_rec_freq:",    sol_record_freq);

  // ===== Load Analysis Data Structure =====
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file);

  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(), "Error: Assigned CPU number does not match the partition. \n");

  delete PartBasic;

  // ===== Quadrature rules =====
  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );

   // Global mesh info
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file, rank);

  // Local sub-domain's element indices
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  // Local sub-domain's IEN array
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  // Local sub-domain's nodal indices
  APart_Node * pNode = new APart_Node(part_file, rank);

  // Local sub-domain's nodal (Dirichlet) BC
  ALocal_NodalBC * locnbc = new ALocal_NodalBC(part_file, rank);

  // Local sub-domain's wall elemental (Neumann) BC for CMM
  ALocal_EBC * locebc_wall = new ALocal_EBC_wall(part_file, rank, quads, "ebc_wall");

  // Control points' xyz coordinates
  FEANode * fNode = new FEANode(part_file, rank);

  // ===== Finite element containers =====
  SYS_T::commPrint("===> Set up volumetric and surface element containers. \n");
  FEAElement * elementv = nullptr;
  FEAElement * elements = nullptr;
  FEAElement * elementw = nullptr;

  if( GMIptr->get_elemType() == 501 )          // linear tet
    elementw = new FEAElement_Triangle3_membrane( nqp_tri );
  else if( GMIptr->get_elemType() == 502 )     // quadratic tet
    elementw = new FEAElement_Triangle6_membrane( nqp_tri );
  else SYS_T::print_fatal("Error: Element type not supported.\n");

  // ===== Generate a sparse matrix for enforcing nodal BCs ====
  Matrix_PETSc * pmat = new Matrix_PETSc(pNode, locnbc);

  pmat -> gen_perm_bc(pNode, locnbc);

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

  // ===== Global assembly =====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");

  IPGAssem * gloAssem_ptr = new PGAssem_Tet_Wall( locAssem_ptr,
      GMIptr, locElem, locIEN, pNode, locnbc, locebc_wall, nz_estimate );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locIEN, locnbc );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();

  // ===== Linear solver context =====
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();

  PC upc; lsolver->GetPC(&upc);
  const PetscInt pfield[1] = {0}, vfields[] = {1,2,3};
  PCFieldSplitSetBlockSize(upc,4);
  PCFieldSplitSetFields(upc,"u",3,vfields,vfields);
  PCFieldSplitSetFields(upc,"p",1,pfield,pfield);

  // ===== Nonlinear solver context =====
  PNonlinear_CMM_Solver * nsolver = new PNonlinear_CMM_Solver( pNode,
      nl_rtol, nl_atol, nl_dtol, nl_maxits, nl_refreq, nl_threshold );

  nsolver->print_info();

  // ===== Temporal solver context =====
  PTime_CMM_Solver * tsolver = new PTime_CMM_Solver( sol_bName,
      sol_record_freq, ttan_renew_freq, final_time );

  tsolver->print_info();

  // ===== FEM analysis =====
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");

  ICVFlowRate * inflow_rate_ptr = nullptr;
  IGenBC * gbc = nullptr;
  ALocal_Inflow_NodalBC * locinfnbc = nullptr;
  ALocal_Ring_NodalBC * locringnbc = nullptr;
  ALocal_EBC * locebc = nullptr;
  IQuadPts * quadv = nullptr;

  tsolver->TM_Prestress( prestress_disp_tol, base, dot_sol, sol, dot_sol_wall_disp, sol_wall_disp,
      tm_galpha_ptr, timeinfo, inflow_rate_ptr, locElem, locIEN, pNode, fNode,
      locnbc, locinfnbc, locringnbc, locebc, locebc_wall, gbc, pmat, elementv, elements, elementw,
      quadv, quads, locAssem_ptr, gloAssem_ptr, lsolver, nsolver );


  delete locElem;
  delete GMIptr; delete quads; delete elementw; delete locebc_wall; delete tm_galpha_ptr;
  delete pNode; delete locAssem_ptr; delete base; delete sol; delete dot_sol;
  delete sol_wall_disp; delete dot_sol_wall_disp; delete timeinfo; delete gloAssem_ptr;
  delete tsolver; delete nsolver; delete lsolver;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
