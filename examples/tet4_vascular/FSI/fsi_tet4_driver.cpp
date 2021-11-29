// ==================================================================
// fsi_tet4_driver.cpp
//
// Author: Ju Liu
// Date: Apr. 4 2019
// ==================================================================
#include "HDF5_Writer.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "APart_Node_FSI.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Tet4.hpp"
#include "CVFlowRate_Unsteady.hpp"
#include "CVFlowRate_Linear2Steady.hpp"
#include "CVFlowRate_Steady.hpp"
#include "GenBC_Resistance.hpp"
#include "GenBC_RCR.hpp"
#include "GenBC_Inductance.hpp"
#include "GenBC_Coronary.hpp"
#include "GenBC_Pressure.hpp"
#include "MaterialModel_NeoHookean_M94_Mixed.hpp"
#include "MaterialModel_NeoHookean_Incompressible_Mixed.hpp"
#include "PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha.hpp"
#include "PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha.hpp"
#include "PLocAssem_Tet4_VMS_Seg_Incompressible.hpp"
#include "PLocAssem_Tet4_FSI_Mesh_Elastostatic.hpp"
#include "PGAssem_FSI_FEM.hpp"
#include "PGAssem_Seg_FEM.hpp"
#include "PDNSolution_Mixed_UPV_3D.hpp"
#include "PTime_Seg_Solver.hpp"

int main(int argc, char *argv[])
{
  int nqp_tet = 5, nqp_tri = 4;

  // Estimate the nonzero per row for the sparse matrix
  int nz_estimate = 300;

  // fluid properties
  double fluid_density = 1.06;
  double fluid_mu = 4.0e-2;

  // solid properties
  double solid_density = 1.0;
  double solid_E = 2.0e6;
  double solid_nu = 0.5;

  // mesh motion elasticity solver parameters
  double mesh_E  = 1.0;
  double mesh_nu = 0.3;

  // flag for determining inflow type 0 pulsatile flow;
  // 1 linear-to-steady; 2 steady
  int inflow_type = 0;

  // inflow file
  std::string inflow_file("inflow_fourier_series.txt");

  double inflow_thd_time = 1.0;      // time for linearly increasing inflow to reach steady state

  // LPN file
  std::string lpn_file("lpn_rcr_input.txt");

  // back flow stabilization
  double bs_beta = 0.2;

  // Generalized-alpha method
  double genA_rho_inf = 0.5;
  bool is_backward_Euler = false;

  // part file location
  const std::string part_file("./apart/part");

  // Nonlinear solver parameters
  double nl_rtol = 1.0e-3;
  double nl_atol = 1.0e-6;
  double nl_dtol = 10.0;
  int nl_maxits = 20;
  int nl_refreq = 4;

  // Time stepping parameters
  double initial_time = 0.0;
  double initial_step = 0.1;
  int initial_index = 0;
  double final_time = 1.0;
  std::string sol_bName("SOL_");
  int ttan_renew_freq = 1;
  int sol_record_freq = 1;

  // Restart options
  bool is_restart = false;
  int restart_index = 0;
  double restart_time = 0.0;
  double restart_step = 1.0e-3;
  std::string restart_name = "SOL_";

  // ===== Initialization of PETSc =====
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();
 
  SYS_T::commPrint("Job starts at %s %s \n", SYS_T::get_time().c_str(), SYS_T::get_date().c_str());

  SYS_T::print_perigee_art();

  // ===== Command Line Argument =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");

  SYS_T::GetOptionInt(   "-nqp_tet",           nqp_tet);
  SYS_T::GetOptionInt(   "-nqp_tri",           nqp_tri);
  SYS_T::GetOptionInt(   "-nz_estimate",       nz_estimate);
  SYS_T::GetOptionReal(  "-bs_beta",           bs_beta);
  SYS_T::GetOptionReal(  "-rho_inf",           genA_rho_inf);
  SYS_T::GetOptionBool(  "-is_backward_Euler", is_backward_Euler);
  SYS_T::GetOptionReal(  "-fl_density",        fluid_density);
  SYS_T::GetOptionReal(  "-fl_mu",             fluid_mu);
  SYS_T::GetOptionReal(  "-sl_density",        solid_density);
  SYS_T::GetOptionReal(  "-sl_E",              solid_E);
  SYS_T::GetOptionReal(  "-sl_nu",             solid_nu);
  SYS_T::GetOptionReal(  "-mesh_E",            mesh_E);
  SYS_T::GetOptionReal(  "-mesh_nu",           mesh_nu);
  SYS_T::GetOptionInt(   "-inflow_type",       inflow_type);
  SYS_T::GetOptionString("-inflow_file",       inflow_file);
  SYS_T::GetOptionReal(  "-inflow_thd_time",   inflow_thd_time);
  SYS_T::GetOptionString("-lpn_file",          lpn_file);
  SYS_T::GetOptionReal(  "-nl_rtol",           nl_rtol);
  SYS_T::GetOptionReal(  "-nl_atol",           nl_atol);
  SYS_T::GetOptionReal(  "-nl_dtol",           nl_dtol);
  SYS_T::GetOptionInt(   "-nl_maxits",         nl_maxits);
  SYS_T::GetOptionInt(   "-nl_refreq",         nl_refreq);
  SYS_T::GetOptionReal(  "-init_time",         initial_time);
  SYS_T::GetOptionReal(  "-fina_time",         final_time);
  SYS_T::GetOptionReal(  "-init_step",         initial_step);
  SYS_T::GetOptionInt(   "-init_index",        initial_index);
  SYS_T::GetOptionInt(   "-ttan_freq",         ttan_renew_freq);
  SYS_T::GetOptionInt(   "-sol_rec_freq",      sol_record_freq);
  SYS_T::GetOptionString("-sol_name",          sol_bName);
  SYS_T::GetOptionBool(  "-is_restart",        is_restart);
  SYS_T::GetOptionInt(   "-restart_index",     restart_index);
  SYS_T::GetOptionReal(  "-restart_time",      restart_time);
  SYS_T::GetOptionReal(  "-restart_step",      restart_step);
  SYS_T::GetOptionString("-restart_name",      restart_name);

  // ===== Print the command line argumetn on screen =====
  SYS_T::cmdPrint("-nqp_tet:", nqp_tet);
  SYS_T::cmdPrint("-nqp_tri:", nqp_tri);
  SYS_T::cmdPrint("-nz_estimate:", nz_estimate);
  SYS_T::cmdPrint("-bs_beta:", bs_beta);
  SYS_T::cmdPrint("-fl_density:", fluid_density);
  SYS_T::cmdPrint("-fl_mu:", fluid_mu);
  SYS_T::cmdPrint("-sl_density:", solid_density);
  SYS_T::cmdPrint("-sl_E:", solid_E);
  SYS_T::cmdPrint("-sl_nu:", solid_nu);
  SYS_T::cmdPrint("-mesh_E:", mesh_E);
  SYS_T::cmdPrint("-mesh_nu:", mesh_nu);
  
  if( is_backward_Euler )
    SYS_T::commPrint(     "-is_backward_Euler: true \n");
  else
    SYS_T::cmdPrint(      "-rho_inf:",         genA_rho_inf);
   
  if( inflow_type == 0 )
  {
    SYS_T::commPrint(   "-inflow_type: 0 (pulsatile flow) \n");
    SYS_T::cmdPrint(    "-inflow_file:",     inflow_file);
  }
  else if( inflow_type == 1 )
  {
    SYS_T::commPrint(   "-inflow_type: 1 (linear-to-steady flow) \n");
    SYS_T::cmdPrint(    "-inflow_file:",     inflow_file);
    SYS_T::cmdPrint(    "-inflow_thd_time:", inflow_thd_time);
  }
  else if( inflow_type == 2 )
  {
    SYS_T::commPrint(   "-inflow_type: 2 (steady flow) \n");
    SYS_T::cmdPrint(    "-inflow_file:",     inflow_file);
  }
  else
    SYS_T::print_fatal("Error: unrecognized inflow_type = %d. \n", inflow_type);  

  SYS_T::cmdPrint("-lpn_file:", lpn_file);
  SYS_T::cmdPrint("-nl_rtol:", nl_rtol);
  SYS_T::cmdPrint("-nl_atol:", nl_atol);
  SYS_T::cmdPrint("-nl_dtol:", nl_dtol);
  SYS_T::cmdPrint("-nl_maxits:", nl_maxits);
  SYS_T::cmdPrint("-nl_refreq:", nl_refreq);
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

  // ===== Record important parameters =====
  if(rank == 0)
  {
    hid_t cmd_file_id = H5Fcreate("solver_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

    cmdh5w->write_doubleScalar(  "fl_density",      fluid_density);
    cmdh5w->write_doubleScalar(  "fl_mu",           fluid_mu);
    cmdh5w->write_doubleScalar(  "sl_density",      solid_density);
    cmdh5w->write_doubleScalar(  "sl_E",            solid_E);
    cmdh5w->write_doubleScalar(  "sl_nu",           solid_nu);
    cmdh5w->write_doubleScalar(  "mesh_E",          mesh_E);
    cmdh5w->write_doubleScalar(  "mesh_nu",         mesh_nu);
    cmdh5w->write_doubleScalar(  "init_step",       initial_step);
    
    cmdh5w->write_string(        "lpn_file",        lpn_file);

    cmdh5w->write_intScalar(     "inflow_type",     inflow_type);
    cmdh5w->write_string(        "inflow_file",     inflow_file);
    if( inflow_type == 1 )
      cmdh5w->write_doubleScalar("inflow_thd_time", inflow_thd_time );

    cmdh5w->write_string("date",              SYS_T::get_date() );
    cmdh5w->write_string("time",              SYS_T::get_time() );

    delete cmdh5w; H5Fclose(cmd_file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== Main Data Strucutre =====
  FEANode * fNode = new FEANode(part_file, rank);
  
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);
  
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);
  
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
  
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);
  
  ALocal_NodalBC * locnbc = new ALocal_NodalBC(part_file, rank);
 
  ALocal_Inflow_NodalBC * locinfnbc = new ALocal_Inflow_NodalBC(part_file, rank);
   
  ALocal_NodalBC * mesh_locnbc = new ALocal_NodalBC(part_file, rank, "mesh_nbc");
  
  ALocal_EBC * locebc = new ALocal_EBC_outflow(part_file, rank);
  
  ALocal_EBC * mesh_locebc = new ALocal_EBC(part_file, rank, "mesh_ebc");
  
  APart_Node * pNode = new APart_Node_FSI(part_file, rank);
  
  SYS_T::commPrint("===> Mesh HDF5 files are read from disk.\n");

  // ===== Basic Checking =====
  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for FEM analysis. \n", size);

  // ===== Inflow rate function =====
  SYS_T::commPrint("===> Setup inflow flow rate. \n");
  
  ICVFlowRate * inflow_rate_ptr = nullptr;
 
  if( inflow_type == 0 )
    inflow_rate_ptr = new CVFlowRate_Unsteady( inflow_file );
  else if( inflow_type == 1 )
    inflow_rate_ptr = new CVFlowRate_Linear2Steady( inflow_thd_time, inflow_file );
  else if( inflow_type == 2 )
    inflow_rate_ptr = new CVFlowRate_Steady( inflow_file );
  else
    SYS_T::print_fatal("Error: unrecognized inflow_type = %d. \n", inflow_type);
  
  inflow_rate_ptr->print_info();

  SYS_T::print_fatal_if(locinfnbc->get_num_nbc() != inflow_rate_ptr->get_num_nbc(),
      "Error: ALocal_Inflow_NodalBC number of faces does not match with that in ICVFlowRate.\n");

  // ===== Quadrature rules and FEM container =====
  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );

  SYS_T::commPrint("===> Setup element container. \n");
  FEAElement * elementv = new FEAElement_Tet4( nqp_tet );
  FEAElement * elements = new FEAElement_Triangle3_3D_der0( nqp_tri );

  // ===== Generate a sparse matrix for strong enforcement of essential BC
  Matrix_PETSc * pmat = new Matrix_PETSc(pNode, locnbc);
  pmat->gen_perm_bc(pNode, locnbc);

  Matrix_PETSc * mmat = new Matrix_PETSc(pNode, mesh_locnbc);
  mmat->gen_perm_bc(pNode, mesh_locnbc);

  // ===== Generate the generalized-alpha method
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");
  
  TimeMethod_GenAlpha * tm_galpha_ptr = nullptr;

  if( is_backward_Euler )
    tm_galpha_ptr = new TimeMethod_GenAlpha( 1.0, 1.0, 1.0 );
  else
    tm_galpha_ptr = new TimeMethod_GenAlpha( genA_rho_inf, false );
  
  tm_galpha_ptr->print_info();
  
  // ===== Local assembly =====
  IPLocAssem * locAssem_fluid_ptr = new PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha(
      tm_galpha_ptr, quadv->get_num_quadPts(), fluid_density, fluid_mu, bs_beta );

  IMaterialModel * matmodel = nullptr;
  IPLocAssem * locAssem_solid_ptr = nullptr;

  if( solid_nu == 0.5 )
  {
    matmodel = new MaterialModel_NeoHookean_Incompressible_Mixed( solid_density, solid_E );

    locAssem_solid_ptr = new PLocAssem_Tet4_VMS_Seg_Incompressible(
        matmodel, tm_galpha_ptr, quadv->get_num_quadPts() );
  }
  else
  {
    matmodel = new MaterialModel_NeoHookean_M94_Mixed( solid_density, solid_E, solid_nu );

    locAssem_solid_ptr = new PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha(
        matmodel, tm_galpha_ptr, quadv->get_num_quadPts() );
  }

  // Pseudo elastic mesh motion
  IPLocAssem * locAssem_mesh_ptr = new PLocAssem_Tet4_FSI_Mesh_Elastostatic( mesh_E, mesh_nu );

  // ===== Initial condition =====
  PDNSolution * base = new PDNSolution_Mixed_UPV_3D( pNode, fNode, locinfnbc, 1 );

  PDNSolution * sol = new PDNSolution_Mixed_UPV_3D( pNode, fNode, 0 );

  PDNSolution * dot_sol = new PDNSolution_Mixed_UPV_3D( pNode, fNode, 0 );

  if( is_restart )
  {
    initial_index = restart_index;
    initial_time  = restart_time;
    initial_step  = restart_step;

    // Read sol file
    SYS_T::file_check(restart_name.c_str());
    sol->ReadBinary(restart_name.c_str());

    // Read dot_sol file
    std::string restart_dot_name = "dot_";
    restart_dot_name.append(restart_name);
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

  // ===== GenBC =====
  IGenBC * gbc = nullptr;

  if( GENBC_T::get_genbc_file_type( lpn_file ) == 1  )
    gbc = new GenBC_Resistance( lpn_file );
  else if( GENBC_T::get_genbc_file_type( lpn_file ) == 2  )
    gbc = new GenBC_RCR( lpn_file, 1000, initial_step );
  else if( GENBC_T::get_genbc_file_type( lpn_file ) == 3  )
    gbc = new GenBC_Inductance( lpn_file );
  else if( GENBC_T::get_genbc_file_type( lpn_file ) == 4  )
    gbc = new GenBC_Coronary( lpn_file, 1000, initial_step, initial_index );
  else if( GENBC_T::get_genbc_file_type( lpn_file ) == 5  )
    gbc = new GenBC_Pressure( lpn_file, initial_time );
  else
    SYS_T::print_fatal( "Error: GenBC input file %s format cannot be recongnized.\n", lpn_file.c_str() );

  gbc -> print_info(); 

  SYS_T::print_fatal_if(gbc->get_num_ebc() != locebc->get_num_ebc(),
      "Error: GenBC number of faces does not match with that in ALocal_EBC.\n");

  // ===== Global assembly routine =====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  IPGAssem * gloAssem_ptr = new PGAssem_FSI_FEM( locAssem_fluid_ptr,
      locAssem_solid_ptr, elements, quads, GMIptr, locElem, locIEN, 
      pNode, locnbc, locebc, gbc, nz_estimate );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_fluid_ptr,
      locAssem_solid_ptr, elements, quads, locIEN, pNode, locnbc, locebc, gbc );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();

  // ===== Global assembly for mesh motion =====
  SYS_T::commPrint("===> Initializing Mat K_mesh and Vec G_mesh ... \n");
  IPGAssem * gloAssem_mesh_ptr = new PGAssem_Seg_FEM( locAssem_mesh_ptr,
      GMIptr, locElem, locIEN, pNode, mesh_locnbc, mesh_locebc );

  SYS_T::commPrint("===> Assembly nonzero estimate for K_mesh ... \n");
  gloAssem_mesh_ptr->Assem_nonzero_estimate( locElem, locAssem_mesh_ptr, locIEN,
      pNode, mesh_locnbc );

  SYS_T::commPrint("===> Matrix K_mesh nonzero structure fixed. \n");
  gloAssem_mesh_ptr->Fix_nonzero_err_str();
  gloAssem_mesh_ptr->Clear_KG();

  // ===== Initialize dot sol =====
  if( is_restart == false )
  {
    SYS_T::commPrint("===> Assembly mass matrix and residual vector.\n");
    PLinear_Solver_PETSc * lsolver_acce = new PLinear_Solver_PETSc(
        1.0e-14, 1.0e-85, 1.0e30, 1000, "mass_", "mass_" );

    KSPSetType(lsolver_acce->ksp, KSPGMRES);
    KSPGMRESSetOrthogonalization(lsolver_acce->ksp,
        KSPGMRESModifiedGramSchmidtOrthogonalization);
    KSPGMRESSetRestart(lsolver_acce->ksp, 500);

    lsolver_acce->Monitor();

    PC preproc; lsolver_acce->GetPC(&preproc);
    PCSetType( preproc, PCHYPRE );
    PCHYPRESetType( preproc, "boomeramg" );

    gloAssem_ptr->Assem_mass_residual( sol, locElem, locAssem_fluid_ptr,
        locAssem_solid_ptr, elementv,
        elements, quadv, quads, locIEN, pNode, fNode, locnbc, locebc );

    PDNSolution * dot_pres_velo = new PDNSolution_P_V_Mixed_3D( pNode, fNode, 0, false );

    lsolver_acce->Solve( gloAssem_ptr->K, gloAssem_ptr->G, dot_pres_velo);

    SEG_SOL_T::PlusAiPV(0.0, -1.0, -1.0, dot_pres_velo, dot_sol);

    SEG_SOL_T::PlusAiVPV(1.0, 0.0, 0.0, sol, dot_sol);

    delete lsolver_acce; delete dot_pres_velo;
    SYS_T::commPrint("\n===> Consistent initial acceleration is obtained.");
    SYS_T::commPrint("\n===> The mass matrix lsolver is destroyed.\n");
  }

  // ===== Linear and nonlinear solver context =====
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();

  PC upc; lsolver->GetPC(&upc);
  const PetscInt pfield[1] = {0}, vfields[] = {1,2,3};
  PCFieldSplitSetBlockSize(upc,4);
  PCFieldSplitSetFields(upc,"u",3,vfields,vfields);
  PCFieldSplitSetFields(upc,"p",1,pfield,pfield);

  PLinear_Solver_PETSc * mesh_lsolver = new PLinear_Solver_PETSc(
      1.0e-12, 1.0e-55, 1.0e30, 500, "mesh_", "mesh_" );

  gloAssem_mesh_ptr->Assem_tangent_residual( sol, sol, 0.0,
      timeinfo->get_step(), locElem, locAssem_mesh_ptr, elementv,
      elements, quadv, quads, locIEN, pNode, fNode, mesh_locnbc,
      mesh_locebc );

  mesh_lsolver -> SetOperator( gloAssem_mesh_ptr->K );
  PC mesh_pc; mesh_lsolver->GetPC(&mesh_pc);
  PCFieldSplitSetBlockSize( mesh_pc, 3 );

  SYS_T::commPrint("===> mesh solver LHS setted up.\n");

  // ===== Nonlinear solver context =====
  PNonlinear_Seg_Solver * nsolver = new PNonlinear_Seg_Solver(
      pNode, fNode, nl_rtol, nl_atol, nl_dtol, nl_maxits, nl_refreq);
  SYS_T::commPrint("===> Nonlinear solver setted up:\n");
  nsolver->print_info();

  // ===== Temporal solver context =====
  PTime_Seg_Solver * tsolver = new PTime_Seg_Solver( sol_bName,
      sol_record_freq, ttan_renew_freq, final_time );
  SYS_T::commPrint("===> Time marching solver setted up:\n");
  tsolver->print_info();

  // ===== Outlet flowrate recording files =====
  for(int ff=0; ff<locebc->get_num_ebc(); ++ff)
  {
    const double dot_face_flrate = gloAssem_ptr -> Assem_surface_flowrate(
        dot_sol, locAssem_fluid_ptr, elements, quads, locebc, ff );

    const double face_flrate = gloAssem_ptr -> Assem_surface_flowrate(
        sol, locAssem_fluid_ptr, elements, quads, locebc, ff );

    const double face_avepre = gloAssem_ptr -> Assem_surface_ave_pressure(
        sol, locAssem_fluid_ptr, elements, quads, locebc, ff );

    // set the gbc initial conditions using the 3D data
    gbc -> reset_initial_sol( ff, face_flrate, face_avepre, timeinfo->get_time(), is_restart );

    const double dot_lpn_flowrate = dot_face_flrate;
    const double lpn_flowrate = face_flrate;
    const double lpn_pressure = gbc -> get_P( ff, dot_lpn_flowrate, lpn_flowrate, timeinfo->get_time() );

    if(rank == 0)
    {
      std::ofstream ofile;

      // If this is NOT a restart run, generate a new file, otherwise append to
      // a existing file
      if( !is_restart )
        ofile.open( locebc->gen_flowfile_name(ff).c_str(), std::ofstream::out | std::ofstream::trunc );
      else
        ofile.open( locebc->gen_flowfile_name(ff).c_str(), std::ofstream::out | std::ofstream::app );

      // if this is NOT a restart run, record the initial values
      if( !is_restart )
      {
        ofile<<"Time-index"<<'\t'<<"Time"<<'\t'<<"Flow-rate"<<'\t'<<"Face-averaged-pressure"<<'\t'<<"Reduced-model-pressure"<<'\n';
        ofile<<timeinfo->get_index()<<'\t'<<timeinfo->get_time()<<'\t'<<face_flrate<<'\t'<<face_avepre<<'\t'<<lpn_pressure<<'\n';
      }

      ofile.close();
    }
  }

  // Write all 0D solutions into a file
  if( rank == 0 ) gbc -> write_0D_sol ( initial_index, initial_time );

  // ===== Inlet data recording files =====
  for(int ff=0; ff<locinfnbc->get_num_nbc(); ++ff)
  {
    const double inlet_face_flrate = gloAssem_ptr -> Assem_surface_flowrate(
        sol, locAssem_fluid_ptr, elements, quads, locinfnbc, ff );

    const double inlet_face_avepre = gloAssem_ptr -> Assem_surface_ave_pressure(
        sol, locAssem_fluid_ptr, elements, quads, locinfnbc, ff );

    if( rank == 0 )
    {
      std::ofstream ofile;
      if( !is_restart )
        ofile.open( locinfnbc->gen_flowfile_name(ff).c_str(), std::ofstream::out | std::ofstream::trunc );
      else
        ofile.open( locinfnbc->gen_flowfile_name(ff).c_str(), std::ofstream::out | std::ofstream::app );

      if( !is_restart )
      {
        ofile<<"Time-index"<<'\t'<<"Time"<<'\t'<<"Flow-rate"<<'\t'<<"Face-averaged-pressure"<<'\n';
        ofile<<timeinfo->get_index()<<'\t'<<timeinfo->get_time()<<'\t'<<inlet_face_flrate<<'\t'<<inlet_face_avepre<<'\n';
      }

      ofile.close();
    }
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== FEM analysis =====
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
  tsolver->TM_FSI_GenAlpha(is_restart, base, dot_sol, sol, tm_galpha_ptr,
      timeinfo, inflow_rate_ptr, locElem, locIEN, pNode, fNode, locnbc,
      locinfnbc, mesh_locnbc, locebc, mesh_locebc, gbc, pmat, mmat, 
      elementv, elements, quadv, quads, 
      locAssem_fluid_ptr, locAssem_solid_ptr, locAssem_mesh_ptr,
      gloAssem_ptr, gloAssem_mesh_ptr, lsolver, mesh_lsolver, nsolver);

  // ===== PETSc Finalize =====
  delete tsolver; delete nsolver; delete lsolver; delete mesh_lsolver;
  delete gloAssem_ptr; delete gloAssem_mesh_ptr;
  delete timeinfo; delete sol; delete dot_sol; delete base;
  delete locAssem_solid_ptr; delete locAssem_fluid_ptr; delete locAssem_mesh_ptr;
  delete matmodel; delete inflow_rate_ptr;
  delete elementv; delete elements; delete quadv; delete quads;
  delete pmat; delete mmat; delete tm_galpha_ptr;
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic;
  delete locElem; delete locnbc; delete locinfnbc; delete mesh_locnbc;
  delete locebc; delete pNode; delete mesh_locebc; delete gbc;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
