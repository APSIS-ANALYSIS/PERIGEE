// ============================================================================
// driver.cpp
// ============================================================================
#include "HDF5_Writer.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_EBC.hpp"
#include "ALocal_NBC.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss_Hex.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Hex27.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FEAElement_Quad9_3D_der0.hpp"
#include "PLocAssem_Transport_GenAlpha.hpp"
#include "PGAssem_Transport_GenAlpha.hpp"
#include "PNonlinear_LinearPDE_Solver.hpp"
#include "PTime_LinearPDE_Solver.hpp"

int main(int argc, char *argv[])
{
  // material parameters
  double rho = 1.0, cap = 1.0, kap = 1.0;

  // number of quadrature points
  int nqp_vol = 5, nqp_sur = 4;
  int nqp_vol_1D = 2, nqp_sur_1D = 2;

  // generalized-alpha rho_inf
  double genA_rho_inf = 0.5;
  bool is_backward_Euler = false;

  // Estimate of the nonzero per row for the sparse matrix
  int nz_estimate = 300;

  // part file location
  std::string part_file("./apart/part");
  
  // nonlinear solver parameters
  double nl_rtol = 1.0e-3; // convergence criterion relative tolerance
  double nl_atol = 1.0e-6; // convergence criterion absolute tolerance
  double nl_dtol = 10.0;   // divergence criterion
  int nl_maxits = 20;      // maximum number if nonlinear iterations
  int nl_refreq = 4;       // frequency of tangent matrix renewal
  int nl_threshold = 4;    // threshold of tangent matrix renewal

  // time stepping parameters
  double initial_time = 0.0; // time of the initial condition
  double initial_step = 0.1; // time step size
  int initial_index = 0;     // indiex of the initial condition
  double final_time = 1.0;   // final time
  std::string sol_bName("SOL_"); // base name of the solution file
  int ttan_renew_freq = 1;   // frequency of tangent matrix renewal
  int sol_record_freq = 1;   // frequency of recording the solution

  // Restart options
  bool is_restart = false;
  int restart_index = 0;     // restart solution time index
  double restart_time = 0.0; // restart time
  double restart_step = 1.0e-3; // restart simulation time step size
  std::string restart_name = "SOL_"; // restart solution base name

  // Yaml options
  bool is_loadYaml = true;
  std::string yaml_file("./linearPDE_input.yml");

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  
  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::print_perigee_art();

  // ===== Yaml Arguments =====
  SYS_T::GetOptionBool("-is_loadYaml", is_loadYaml);
  SYS_T::GetOptionString("-yaml_file", yaml_file);

  if (is_loadYaml) SYS_T::InsertFileYAML( yaml_file,  false );

  // ===== Read Command Line Arguments =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionReal("-rho", rho);
  SYS_T::GetOptionReal("-cap", cap);
  SYS_T::GetOptionReal("-kap", kap);
  SYS_T::GetOptionInt("-nqp_vol", nqp_vol);
  SYS_T::GetOptionInt("-nqp_sur", nqp_sur);
  SYS_T::GetOptionInt("-nqp_vol_1d", nqp_vol_1D);
  SYS_T::GetOptionInt("-nqp_sur_1d", nqp_sur_1D);
  SYS_T::GetOptionReal("-rho_inf", genA_rho_inf);
  SYS_T::GetOptionBool("-is_backward_Euler", is_backward_Euler);
  SYS_T::GetOptionInt("-nz_estimate", nz_estimate);
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionReal("-nl_rtol", nl_rtol);
  SYS_T::GetOptionReal("-nl_atol", nl_atol);
  SYS_T::GetOptionReal("-nl_dtol", nl_dtol);
  SYS_T::GetOptionInt("-nl_maxits", nl_maxits);
  SYS_T::GetOptionInt("-nl_refreq", nl_refreq);
  SYS_T::GetOptionInt("-nl_threshold", nl_threshold);
  SYS_T::GetOptionReal("-init_time", initial_time);
  SYS_T::GetOptionReal("-fina_time", final_time);
  SYS_T::GetOptionReal("-init_step", initial_step);
  SYS_T::GetOptionInt("-init_index", initial_index);
  SYS_T::GetOptionInt("-ttan_freq", ttan_renew_freq);
  SYS_T::GetOptionInt("-sol_rec_freq", sol_record_freq);
  SYS_T::GetOptionString("-sol_name", sol_bName);
  SYS_T::GetOptionBool("-is_restart", is_restart);
  SYS_T::GetOptionInt("-restart_index", restart_index);
  SYS_T::GetOptionReal("-restart_time", restart_time);
  SYS_T::GetOptionReal("-restart_step", restart_step);
  SYS_T::GetOptionString("-restart_name", restart_name);

  // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint("-nqp_vol:", nqp_vol);
  SYS_T::cmdPrint("-nqp_sur:", nqp_sur);
  SYS_T::cmdPrint("-nqp_vol_1d", nqp_vol_1D);
  SYS_T::cmdPrint("-nqp_sur_1d", nqp_sur_1D);
  if( is_backward_Euler )
    SYS_T::commPrint(   "-is_backward_Euler: true \n");
  else
    SYS_T::cmdPrint(    "-rho_inf:",         genA_rho_inf);
  
  SYS_T::cmdPrint("-nz_estimate:", nz_estimate);
  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-nl_rtol:", nl_rtol);
  SYS_T::cmdPrint("-nl_atol:", nl_atol);
  SYS_T::cmdPrint("-nl_dtol:", nl_dtol);
  SYS_T::cmdPrint("-nl_maxits:", nl_maxits);
  SYS_T::cmdPrint("-nl_refreq:", nl_refreq);
  SYS_T::cmdPrint("-nl_threshold:", nl_threshold);
  SYS_T::cmdPrint("-init_time:", initial_time);
  SYS_T::cmdPrint("-init_step:", initial_step);
  SYS_T::cmdPrint("-init_index:", initial_index);
  SYS_T::cmdPrint("-fina_time:", final_time);
  SYS_T::cmdPrint("-ttan_freq:", ttan_renew_freq);
  SYS_T::cmdPrint("-sol_rec_freq:", sol_record_freq);
  SYS_T::cmdPrint("-sol_name:", sol_bName);
  if(is_restart)
  {
    SYS_T::commPrint("-is_restart: true \n");
    SYS_T::cmdPrint("-restart_index:", restart_index);
    SYS_T::cmdPrint("-restart_time:", restart_time);
    SYS_T::cmdPrint("-restart_step:", restart_step);
    SYS_T::cmdPrint("-restart_name:", restart_name);
  }
  else SYS_T::commPrint("-is_restart: false \n");

  // ===== Record important solver options =====
  if(rank == 0)
  {
    hid_t cmd_file_id = H5Fcreate("solver_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

    cmdh5w->write_doubleScalar("init_step", initial_step);
    cmdh5w->write_intScalar("sol_record_freq", sol_record_freq);
    cmdh5w->write_intScalar("nqp_vol", nqp_vol);
    cmdh5w->write_intScalar("nqp_sur", nqp_sur);

    delete cmdh5w; H5Fclose(cmd_file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== Data from Files =====
  // Control points' xyz coordinates
  FEANode * fNode = new FEANode(part_file, rank);

  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);

  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  APart_Node * pNode = new APart_Node(part_file, rank);

  ALocal_NBC * locnbc = new ALocal_NBC(part_file, rank);

  ALocal_EBC * locebc = new ALocal_EBC(part_file, rank);

  SYS_T::commPrint("===> Data from HDF5 files are read from disk.\n");

  SYS_T::print_fatal_if( size != PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for FEM analysis.\n", size);

  // ===== Finite Element Container & Quadrature rules=====
  SYS_T::commPrint("===> Setup element container. \n");
  SYS_T::commPrint("===> Build quadrature rules. \n");
  FEAElement * elementv = nullptr;
  FEAElement * elements = nullptr;
  IQuadPts * quadv = nullptr;
  IQuadPts * quads = nullptr;

  if( GMIptr->get_elemType() == 501 )
  {
    if( nqp_vol > 5 ) SYS_T::commPrint("Warning: the tet element is linear and you are using more than 5 quadrature points.\n");
    if( nqp_sur > 4 ) SYS_T::commPrint("Warning: the tri element is linear and you are using more than 4 quadrature points.\n");

    elementv = new FEAElement_Tet4( nqp_vol ); // elem type 501
    elements = new FEAElement_Triangle3_3D_der0( nqp_sur );
    quadv = new QuadPts_Gauss_Tet( nqp_vol );
    quads = new QuadPts_Gauss_Triangle( nqp_sur );
  }
  else if( GMIptr->get_elemType() == 502 )
  {
    SYS_T::print_fatal_if( nqp_vol < 29, "Error: not enough quadrature points for tets.\n" );
    SYS_T::print_fatal_if( nqp_sur < 13, "Error: not enough quadrature points for triangles.\n" );

    elementv = new FEAElement_Tet10_v2( nqp_vol ); // elem type 502
    elements = new FEAElement_Triangle6_3D_der0( nqp_sur );
    quadv = new QuadPts_Gauss_Tet( nqp_vol );
    quads = new QuadPts_Gauss_Triangle( nqp_sur );
  }
  else if( GMIptr->get_elemType() == 601 )
  {
    SYS_T::print_fatal_if( nqp_vol_1D < 2, "Error: not enough quadrature points for hex.\n" );
    SYS_T::print_fatal_if( nqp_sur_1D < 1, "Error: not enough quadrature points for quad.\n" );

    elementv = new FEAElement_Hex8( nqp_vol_1D * nqp_vol_1D * nqp_sur_1D ); // elem type 601
    elements = new FEAElement_Quad4_3D_der0( nqp_sur_1D * nqp_sur_1D );
    quadv = new QuadPts_Gauss_Hex( nqp_vol_1D );
    quads = new QuadPts_Gauss_Quad( nqp_sur_1D );
  }
  else if( GMIptr->get_elemType() == 602 )
  {
    SYS_T::print_fatal_if( nqp_vol_1D < 4, "Error: not enough quadrature points for hex.\n" );
    SYS_T::print_fatal_if( nqp_sur_1D < 3, "Error: not enough quadrature points for quad.\n" );

    elementv = new FEAElement_Hex27( nqp_vol_1D * nqp_vol_1D * nqp_sur_1D ); // elem type 602
    elements = new FEAElement_Quad9_3D_der0( nqp_sur_1D * nqp_sur_1D );
    quadv = new QuadPts_Gauss_Hex( nqp_vol_1D );
    quads = new QuadPts_Gauss_Quad( nqp_sur_1D );
  }
  else SYS_T::print_fatal("Error: Element type not supported.\n");

  // print the information of element and quadrature rule
  elementv->print_info();
  elements->print_info();
  quadv->print_info();
  quads->print_info();

  // ===== Generate a sparse matrix for the enforcement of essential BCs
  Matrix_PETSc * pmat = new Matrix_PETSc(pNode, locnbc);

  pmat->gen_perm_bc(pNode, locnbc);

  // ===== Generalized-alpha =====
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");

  TimeMethod_GenAlpha * tm_galpha_ptr = nullptr;

  if( is_backward_Euler )
    tm_galpha_ptr = new TimeMethod_GenAlpha( 1.0, 1.0, 1.0 );
  else
    tm_galpha_ptr = new TimeMethod_GenAlpha( genA_rho_inf, false );

  tm_galpha_ptr->print_info();

  // ===== Local Assembly Routine =====
  IPLocAssem * locAssem_ptr = new PLocAssem_Transport_GenAlpha(
      rho, cap, kap, tm_galpha_ptr,
      elementv->get_nLocBas(), elements->get_nLocBas(), 
      locebc -> get_num_ebc());

  // ===== Initial condition =====
  PDNSolution * sol = new PDNSolution_LinearPDE( pNode, 0, true, "temperature solution" );
  
  PDNSolution * dot_sol = new PDNSolution_LinearPDE( pNode, 0, true, "temperature solution" );
  
  if( is_restart )
  {
    initial_index = restart_index;
    initial_time  = restart_time;
    initial_step  = restart_step;

    // Read sol file
    SYS_T::file_check(restart_name.c_str());
    sol->ReadBinary(restart_name.c_str());

    // generate the corresponding dot_sol file name
    std::string restart_dot_name = "dot_";
    restart_dot_name.append(restart_name);

    // Read dot_sol file
    SYS_T::file_check(restart_dot_name.c_str());
    dot_sol->ReadBinary(restart_dot_name.c_str());

    SYS_T::commPrint("===> Read sol from disk as a restart run... \n");
    SYS_T::commPrint("     restart_name: %s \n", restart_name.c_str());
    SYS_T::commPrint("     restart_dot_name: %s \n", restart_dot_name.c_str());
    SYS_T::commPrint("     restart_time: %e \n", restart_time);
    SYS_T::commPrint("     restart_index: %d \n", restart_index);
    SYS_T::commPrint("     restart_step: %e \n", restart_step);
  }

  // ===== Time step info =====
  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);

  // ===== Global assembly =====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  IPGAssem * gloAssem_ptr = new PGAssem_Transport_GenAlpha( locAssem_ptr,
      GMIptr, locElem, locIEN, pNode, locnbc, locebc, nz_estimate );  

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_ptr, locIEN, locnbc );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();
  
  // ===== Initialize the dot_sol vector by solving mass matrix =====
  if( is_restart == false )
  {
    PLinear_Solver_PETSc * lsolver_acce = new PLinear_Solver_PETSc(
        1.0e-14, 1.0e-85, 1.0e30, 1000, "mass_", "mass_" );

    KSPSetType(lsolver_acce->ksp, KSPGMRES);
    KSPGMRESSetOrthogonalization(lsolver_acce->ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);
    KSPGMRESSetRestart(lsolver_acce->ksp, 500);

    PC preproc; lsolver_acce->GetPC(&preproc);
    PCSetType( preproc, PCHYPRE );
    PCHYPRESetType( preproc, "boomeramg" );
   
    SYS_T::commPrint("===> Assembly mass matrix and residual vector.\n"); 
    gloAssem_ptr->Assem_mass_residual( sol, locElem, locAssem_ptr, elementv,
        elements, quadv, quads, locIEN, fNode, locnbc, locebc );

    lsolver_acce->Solve( gloAssem_ptr->K, gloAssem_ptr->G, dot_sol );
  
    dot_sol -> ScaleValue( -1.0 );

    SYS_T::commPrint("\n===> Consistent initial acceleration is obtained. \n");
    lsolver_acce -> print_info();
    delete lsolver_acce;
    SYS_T::commPrint(" The mass matrix lsolver is destroyed.\n");
  }

  // ===== Linear solver context =====
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();
  
  // ===== Nonlinear solver context =====
  PNonlinear_LinearPDE_Solver * nsolver = new PNonlinear_LinearPDE_Solver(
      nl_rtol, nl_atol, nl_dtol, nl_maxits, nl_refreq, nl_threshold );

  nsolver->print_info();

  // ===== Temporal solver context =====
  PTime_LinearPDE_Solver * tsolver = new PTime_LinearPDE_Solver( sol_bName,
      sol_record_freq, ttan_renew_freq, final_time );

  tsolver->print_info();

  
  // ===== FEM analysis =====
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");

  tsolver->TM_GenAlpha_Transport(is_restart, dot_sol, sol,
      tm_galpha_ptr, timeinfo, locElem, locIEN, pNode, fNode,
      locnbc, locebc, pmat, elementv, elements, quadv, quads,
      locAssem_ptr, gloAssem_ptr, lsolver, nsolver);

  // ===== Print complete solver info =====
  lsolver -> print_info();

  delete tsolver; delete nsolver;
  delete lsolver; delete gloAssem_ptr; delete dot_sol; delete timeinfo;
  delete fNode; delete locIEN; delete GMIptr; delete locElem; delete pNode; delete PartBasic;
  delete locnbc; delete locebc; delete quadv; delete quads; delete elementv; delete elements;
  delete pmat; delete tm_galpha_ptr; delete locAssem_ptr; delete sol;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
