// ==================================================================
// ale_ns_tet4_driver.cpp
// 
// 4-node tetrahedral finite element analysis code for 3D Navier-Stokes
// equations with variational multiscale analysis.
//
// Ref: CMAME 197 (2007) 173-201
// Date: June 10 2017
// Author: Ju Liu
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_Inflow_NodalBC.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Tet4.hpp"
#include "CVFlowRate_Unsteady.hpp"
#include "CVFlowRate_Linear2Steady.hpp"
#include "GenBC_Resistance.hpp"
#include "GenBC_RCR.hpp"
#include "GenBC_Inductance.hpp"
#include "GenBC_Coronary.hpp"
#include "GenBC_Pressure.hpp"
#include "PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha.hpp"
#include "PLocAssem_Tet4_FSI_Mesh_Laplacian.hpp"
#include "PGAssem_ALE_NS_FEM.hpp"
#include "PDNSolution_Tet4_ALE_NS_3D.hpp"
#include "PTime_Seg_Solver.hpp"

int main(int argc, char *argv[])
{
  // Number of quadrature points for tets and triangles
  int nqp_tet = 5, nqp_tri = 4;
  
  // Estimate of the nonzero per row for the sparse matrix 
  int nz_estimate = 60;

  // fluid properties
  double fluid_density = 1.065;
  double fluid_mu = 3.5e-2;
  
  // inflow file
  std::string inflow_file("inflow_fourier_series.txt");

  // LPN file
  std::string lpn_file("lpn_rcr_input.txt");

  // back flow stabilization
  double bs_beta = 0.2;

  // part file location
  std::string part_file("part");

  // determine if the mesh moving
  bool is_ale = true;

  // determine if we want to print solver info
  bool is_LS_info = false;

  // nonlinear solver parameters
  double nl_rtol = 1.0e-3;
  double nl_atol = 1.0e-6;
  double nl_dtol = 10.0;
  int nl_maxits = 20;
  int nl_refreq = 4;

  // time stepping parameters
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

  PetscMPIInt rank, size;
  // ===== Initialization of PETSc =====
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // Print art letters
  SYS_T::print_perigee_art();

  // ===== Command Line Argument =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");

  SYS_T::GetOptionInt("-nqp_tet", nqp_tet);
  SYS_T::GetOptionInt("-nqp_tri", nqp_tri);
  SYS_T::GetOptionInt("-nz_estimate", nz_estimate);
  SYS_T::GetOptionReal("-bs_beta", bs_beta);
  SYS_T::GetOptionReal("-fl_density", fluid_density);
  SYS_T::GetOptionReal("-fl_mu", fluid_mu);
  SYS_T::GetOptionString("-inflow_file", inflow_file);
  SYS_T::GetOptionString("-lpn_file", lpn_file);
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionBool("-is_ale", is_ale);
  SYS_T::GetOptionBool("-is_ls_info", is_LS_info);
  SYS_T::GetOptionReal("-nl_rtol", nl_rtol);
  SYS_T::GetOptionReal("-nl_atol", nl_atol);
  SYS_T::GetOptionReal("-nl_dtol", nl_dtol);
  SYS_T::GetOptionInt("-nl_maxits", nl_maxits);
  SYS_T::GetOptionInt("-nl_refreq", nl_refreq);
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

  // ===== Print Command Line Arguement =====
  SYS_T::cmdPrint("-nqp_tet:", nqp_tet);
  SYS_T::cmdPrint("-nqp_tri:", nqp_tri);
  SYS_T::cmdPrint("-nz_estimate:", nz_estimate);
  SYS_T::cmdPrint("-bs_beta:", bs_beta);
  SYS_T::cmdPrint("-fl_density:", fluid_density);
  SYS_T::cmdPrint("-fl_mu:", fluid_mu);
  SYS_T::cmdPrint("-inflow_file:", inflow_file);
  SYS_T::cmdPrint("-lpn_file:", lpn_file);
  SYS_T::cmdPrint("-part_file:", part_file);
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
    PetscPrintf(PETSC_COMM_WORLD, "-is_restart: true \n");
    SYS_T::cmdPrint("-restart_index:", restart_index);
    SYS_T::cmdPrint("-restart_time:", restart_time);
    SYS_T::cmdPrint("-restart_step:", restart_step);
    SYS_T::cmdPrint("-restart_name:", restart_name);
  }
  else PetscPrintf(PETSC_COMM_WORLD, "-is_restart: false \n");
  
  if(is_ale) PetscPrintf(PETSC_COMM_WORLD, "-is_ale: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-is_ale: false \n");

  if(is_LS_info) PetscPrintf(PETSC_COMM_WORLD, "-is_ls_info: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-is_ls_info: false \n");

  // ===== Main Data Structure =====
  FEANode * fNode = new FEANode(part_file, rank);

  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);

  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  ALocal_NBC * locnbc = new ALocal_NBC(part_file, rank);

  ALocal_Inflow_NodalBC * locinfnbc = new ALocal_Inflow_NodalBC(part_file, rank);

  ALocal_NBC * mesh_locnbc = new ALocal_NBC(part_file, rank, "mesh_nbc");

  ALocal_EBC * locebc = new ALocal_EBC_outflow(part_file, rank);

  ALocal_EBC * mesh_locebc = new ALocal_EBC(part_file, rank, "mesh_ebc");

  APart_Node * pNode = new APart_Node(part_file, rank);

  SYS_T::commPrint("===> Mesh HDF5 files are read from disk.\n");

  // ===== Do Basic checking =====
  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  // ===== Generate Quadrature rules and FEM container =====
  PetscPrintf(PETSC_COMM_WORLD,
      "===> %d processor(s) are assigned for FEM analysis. \n", size);

  // ===== Inflow flow rate function =====
  SYS_T::commPrint("===> Setup inflow flow rate. \n");

  // Unsteady flow rate
  ICVFlowRate * inflow_rate_ptr = new CVFlowRate_Unsteady( inflow_file.c_str() );

  inflow_rate_ptr->print_info();

  // ===== Quadrature points and element container =====
  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );

  SYS_T::commPrint("===> Setup element container. \n");
  FEAElement * elementv = new FEAElement_Tet4( quadv-> get_num_quadPts() );
  FEAElement * elements = new FEAElement_Triangle3_3D_der0(
      quads-> get_num_quadPts() );

  // ===== Generate a sparse matrix for strong enforcement of essential BC
  Matrix_PETSc * pmat = new Matrix_PETSc(pNode, locnbc);
  pmat->gen_perm_bc(pNode, locnbc);

  Matrix_PETSc * mmat = new Matrix_PETSc(pNode, mesh_locnbc);
  mmat->gen_perm_bc(pNode, mesh_locnbc);

  // ===== Generate Generalized-alpha method
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");
  const double genA_spectrium = 0.5;
  const bool genA_is2ndSystem = false;
  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha(
      genA_spectrium, genA_is2ndSystem);
  tm_galpha_ptr->print_info();

  // ===== Local assembly initialization =====
  IPLocAssem * locAssem_ptr = new PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha(
      tm_galpha_ptr, quadv->get_num_quadPts(), fluid_density, fluid_mu, bs_beta );

  IPLocAssem * locAssem_mesh_ptr = new PLocAssem_Tet4_FSI_Mesh_Laplacian();

  // ===== Initial condition =====
  PDNSolution * base = new PDNSolution_Tet4_ALE_NS_3D( pNode, fNode, locinfnbc, 1 );

  PDNSolution * sol = new PDNSolution_Tet4_ALE_NS_3D( pNode, 0 );

  PDNSolution * dot_sol = new PDNSolution_Tet4_ALE_NS_3D( pNode, 0 );

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
    
    PetscPrintf(PETSC_COMM_WORLD, "===> Read sol from disk as a restart run... \n");
    PetscPrintf(PETSC_COMM_WORLD, "     restart_name: %s \n", restart_name.c_str());
    PetscPrintf(PETSC_COMM_WORLD, "     restart_dot_name: %s \n", restart_dot_name.c_str());
    PetscPrintf(PETSC_COMM_WORLD, "     restart_time: %e \n", restart_time);
    PetscPrintf(PETSC_COMM_WORLD, "     restart_index: %d \n", restart_index);
    PetscPrintf(PETSC_COMM_WORLD, "     restart_step: %e \n", restart_step);
  }

  // Initialize the time info class
  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step); 

  // ===== GenBC RCR ====
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

  // Make sure gbc number of faces matches that of ALocal_EBC
  SYS_T::print_fatal_if(gbc->get_num_ebc() != locebc->get_num_ebc(),
      "Error: GenBC number of faces does not match with that in ALocal_EBC.\n");

  // ===== Global Assembly Routine ====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  IPGAssem * gloAssem_ptr = new PGAssem_ALE_NS_FEM( locAssem_ptr, elements, quads, 
      GMIptr, locElem, locIEN, pNode, locnbc, locebc, gbc, nz_estimate );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_ptr, 
      elements, quads, locIEN, pNode, locnbc, locebc, gbc );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();

  // ===== Global Assemlby for Mesh motion =====
  IPGAssem * gloAssem_mesh_ptr = NULL;
  if( is_ale == true )
  { 
    SYS_T::commPrint("===> Initializing Mat K_mesh and Vec G_mesh ... \n");
    gloAssem_mesh_ptr = new PGAssem_ALE_NS_FEM( locAssem_mesh_ptr,
        GMIptr, locElem, locIEN, pNode, mesh_locnbc, mesh_locebc );

    SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
    gloAssem_mesh_ptr->Assem_nonzero_estimate( locElem, locAssem_mesh_ptr, locIEN,
        pNode, mesh_locnbc );

    SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
    gloAssem_mesh_ptr->Fix_nonzero_err_str();
    gloAssem_mesh_ptr->Clear_KG();
  }

  // ===== Initialize the dot_sol vector by solving mass matrix =====
  if( is_restart == false )
  {
    SYS_T::commPrint("===> Assembly mass matrix and residual vector.\n");
    PLinear_Solver_PETSc * lsolver_acce = new PLinear_Solver_PETSc(
        1.0e-14, 1.0e-85, 1.0e30, 1000, "ls_mass_", "pc_mass_" );

    KSPSetType(lsolver_acce->ksp, KSPGMRES);
    KSPGMRESSetOrthogonalization(lsolver_acce->ksp,
        KSPGMRESModifiedGramSchmidtOrthogonalization);
    KSPGMRESSetRestart(lsolver_acce->ksp, 500);

    PC preproc; lsolver_acce->GetPC(&preproc);
    PCSetType( preproc, PCHYPRE );
    PCHYPRESetType( preproc, "boomeramg" );

    PDNSolution * dot_pres_velo = new PDNSolution_P_V_Mixed_3D( pNode, fNode, 0 );

    gloAssem_ptr->Assem_mass_residual( sol, locElem, locAssem_ptr, elementv,
        elements, quadv, quads, locIEN, pNode, fNode, locnbc, locebc );

    lsolver_acce->Solve( gloAssem_ptr->K, gloAssem_ptr->G, dot_pres_velo );

    SEG_SOL_T::PlusAiPV(0.0, -1.0, -1.0, dot_pres_velo, dot_sol);

    delete lsolver_acce; delete dot_pres_velo;
    SYS_T::commPrint("\n===> Consistent initial acceleration is obtained.");
    SYS_T::commPrint(" The mass matrix lsolver is destroyed. \n\n");
  }

  // ===== Linear & Nonlinear & Time solvers =====
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();

  PC upc; lsolver->GetPC(&upc);
  const PetscInt pfield[1] = {0}, vfields[] = {1,2,3};
  PCFieldSplitSetBlockSize(upc,4);
  PCFieldSplitSetFields(upc,"u",3,vfields,vfields);
  PCFieldSplitSetFields(upc,"p",1,pfield,pfield);

  // ===== Mesh motion solver =====
  PLinear_Solver_PETSc * mesh_lsolver = NULL;
  if( is_ale == true )
  { 
    mesh_lsolver = new PLinear_Solver_PETSc(
        1.0e-12, 1.0e-55, 1.0e30, 500, "ls_mesh_", "pc_mesh_" );
    KSPSetType( mesh_lsolver->ksp, KSPGMRES );

    gloAssem_mesh_ptr->Assem_tangent_residual( dot_sol, sol, 
        locElem, locAssem_mesh_ptr, elementv, elements, quadv, quads, 
        locIEN, pNode, fNode, mesh_locnbc, mesh_locebc );

    mesh_lsolver -> SetOperator( gloAssem_mesh_ptr->K );

    PC mesh_pc; mesh_lsolver->GetPC(&mesh_pc);
    PCSetType( mesh_pc, PCHYPRE );
    SYS_T::commPrint("===> mesh solver LHS setted up.\n");
  }

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
        dot_sol, locAssem_ptr, elements, quads, locebc, ff );

    const double face_flrate = gloAssem_ptr -> Assem_surface_flowrate( 
        sol, locAssem_ptr, elements, quads, locebc, ff );

    const double face_avepre = gloAssem_ptr -> Assem_surface_ave_pressure( 
        sol, locAssem_ptr, elements, quads, locebc, ff );

    // set the gbc initial conditions using the 3D data
    gbc -> reset_initial_sol( ff, face_flrate, face_avepre, timeinfo->get_time(), is_restart );

    const double dot_lpn_flowrate = dot_face_flrate;
    const double lpn_flowrate = face_flrate;
    const double lpn_pressure = gbc -> get_P( ff, dot_lpn_flowrate, lpn_flowrate, timeinfo->get_time() );

    // Create the txt files and write the initial flow rates
    if(rank == 0)
    {
      std::ofstream ofile;
      
      // If this is NOT a restart run, generate a new file, otherwise append to
      // existing file
      if( !is_restart )
        ofile.open( locebc->gen_flowfile_name(ff).c_str(), std::ofstream::out | std::ofstream::trunc );
      else
        ofile.open( locebc->gen_flowfile_name(ff).c_str(), std::ofstream::out | std::ofstream::app );

      // If this is NOT a restart, then record the initial values
      if( !is_restart )
        ofile<<timeinfo->get_index()<<'\t'<<timeinfo->get_time()<<'\t'<<face_flrate<<'\t'<<face_avepre<<'\t'<<lpn_pressure<<'\n';
      
      ofile.close();
    }
  }

  // ===== FEM analysis =====
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");

  tsolver->TM_ALE_NS_GenAlpha(is_restart, is_ale, base, dot_sol, sol, 
      tm_galpha_ptr, timeinfo, inflow_rate_ptr, locElem, locIEN, pNode, fNode,
      locnbc, locinfnbc, mesh_locnbc, locebc, mesh_locebc, gbc, pmat, mmat,
      elementv, elements, quadv, quads,
      locAssem_ptr, locAssem_mesh_ptr,
      gloAssem_ptr, gloAssem_mesh_ptr, lsolver, mesh_lsolver, nsolver);

  // ===== PETSc Finalize =====
  delete tsolver; delete nsolver; delete lsolver; delete timeinfo;
  delete gloAssem_ptr; delete base;  delete sol; delete dot_sol;
  delete locAssem_ptr; delete locAssem_mesh_ptr; delete tm_galpha_ptr;
  delete pmat; delete mmat; delete locinfnbc; delete inflow_rate_ptr;
  delete elementv; delete elements; delete quadv; delete quads;
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic; 
  delete locElem; delete locnbc; delete locebc; delete pNode;
  delete mesh_locnbc; delete mesh_locebc; delete gbc;
  
  if( is_ale == true )
  {
    delete mesh_lsolver;
    delete gloAssem_mesh_ptr;
  }

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF 
