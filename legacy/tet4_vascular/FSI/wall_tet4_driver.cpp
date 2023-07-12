// ============================================================================
// wall_tet4_driver.cpp
//
// wall mechanics solver for generating the prestress.
//
// Date: Sep. 13 2021
// ============================================================================
#include "HDF5_Reader.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "APart_Node_FSI.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Tet4.hpp"
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
#include "PGAssem_Wall_Prestress.hpp"
#include "PGAssem_Seg_FEM.hpp"
#include "PDNSolution_Mixed_UPV_3D.hpp"
#include "PTime_Seg_Solver.hpp"

int main( int argc, char *argv[] )
{
  // solution file name to be loaded for prestressing
  std::string restart_name = "SOL_re";

  // (Pseudo-) time integration parameters
  double genA_rho_inf = 0.0;
  bool is_backward_Euler = true;
  const bool is_load_ps = false;

  // Estimate of num nonzeros per row for the sparse tangent matrix
  int nz_estimate = 300;

  // Prestress tolerance
  double prestress_disp_tol = 1.0e-6;

  // Nonlinear solver parameters
  double nl_rtol = 1.0e-3;           // convergence criterion relative tolerance
  double nl_atol = 1.0e-6;           // convergence criterion absolute tolerance
  double nl_dtol = 1.0e3;            // divergence criterion
  int    nl_maxits = 20;             // maximum number if nonlinear iterations
  int    nl_refreq = 4;              // frequency of tangent matrix renewal
  int    nl_threshold = 4;           // threshold of tangent matrix renewal

  // Time stepping parameters
  double initial_time = 0.0;         // time of initial condition
  double initial_step = 0.1;         // time step size
  int    initial_index = 0;          // index of initial condition
  double final_time = 1.0;           // end time of simulation
  bool   is_record_sol = false;      // bool flag to decide if one wants to record the solution
  std::string sol_bName("PS_");      // base name of the solution file
  int    ttan_renew_freq = 1;        // frequency of tangent matrix renewal
  int    sol_record_freq = 1;        // frequency for recording the solution

  // We assume that a 3D solver has been called (to generate the wall traction)
  // and a suite of command line arguments has been saved to disk
  hid_t solver_cmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( solver_cmd_file );

  const int nqp_tet       = cmd_h5r -> read_intScalar(    "/", "nqp_tet");
  const int nqp_tri       = cmd_h5r -> read_intScalar(    "/", "nqp_tri");
  const double sl_nu      = cmd_h5r -> read_doubleScalar( "/", "sl_nu");

  delete cmd_h5r; H5Fclose(solver_cmd_file);

  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  HDF5_Reader * pcmd_h5r = new HDF5_Reader( prepcmd_file );

  const std::string part_file = pcmd_h5r -> read_string(    "/", "part_file" );
  const int fsiBC_type        = pcmd_h5r -> read_intScalar( "/", "fsiBC_type" );

  delete pcmd_h5r; H5Fclose(prepcmd_file);

  // Initialize PETSc
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::print_fatal_if( fsiBC_type != 2, "Error: fsiBC_type should be 2.\n" );

  // Clean potentially pre-existing hdf5 files of prestress saved in the folder
  // named as prestress
  if(rank == 0 )
  {
    SYS_T::execute("rm -rf prestress");
    SYS_T::execute("mkdir prestress");
  }

  SYS_T::GetOptionInt(   "-nz_estimate",         nz_estimate);
  SYS_T::GetOptionReal(  "-prestress_disp_tol",  prestress_disp_tol);
  SYS_T::GetOptionReal(  "-rho_inf",             genA_rho_inf);
  SYS_T::GetOptionReal(  "-nl_rtol",             nl_rtol);
  SYS_T::GetOptionReal(  "-nl_atol",             nl_atol);
  SYS_T::GetOptionReal(  "-nl_dtol",             nl_dtol);
  SYS_T::GetOptionInt(   "-nl_maxits",           nl_maxits);
  SYS_T::GetOptionInt(   "-nl_refreq",           nl_refreq);
  SYS_T::GetOptionInt(   "-nl_threshold",        nl_threshold);
  SYS_T::GetOptionBool(  "-is_backward_Euler",   is_backward_Euler);
  SYS_T::GetOptionReal(  "-init_time",           initial_time);
  SYS_T::GetOptionReal(  "-fina_time",           final_time);
  SYS_T::GetOptionReal(  "-init_step",           initial_step);
  SYS_T::GetOptionInt(   "-init_index",          initial_index);
  SYS_T::GetOptionInt(   "-ttan_freq",           ttan_renew_freq);
  SYS_T::GetOptionBool(  "-is_record_sol",       is_record_sol);
  SYS_T::GetOptionInt(   "-sol_rec_freq",        sol_record_freq);
  SYS_T::GetOptionString("-restart_name",        restart_name);

  // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint(      "part_file:",          part_file);
  SYS_T::cmdPrint(       "-prestress_disp_tol:", prestress_disp_tol);
  SYS_T::cmdPrint(       "-nl_rtol:",            nl_rtol);
  SYS_T::cmdPrint(       "-nl_atol:",            nl_atol);
  SYS_T::cmdPrint(       "-nl_dtol:",            nl_dtol);
  SYS_T::cmdPrint(       "-nl_maxits:",          nl_maxits);
  SYS_T::cmdPrint(       "-nl_refreq:",          nl_refreq);
  SYS_T::cmdPrint(       "-nl_threshold:",       nl_threshold);

  if( is_backward_Euler )
    SYS_T::commPrint(    "-is_backward_Euler: true \n");
  else
    SYS_T::cmdPrint(     "-rho_inf:",            genA_rho_inf);

  SYS_T::cmdPrint(       "-restart_name:",       restart_name);
  SYS_T::cmdPrint(       "-init_time:",          initial_time);
  SYS_T::cmdPrint(       "-init_step:",          initial_step);
  SYS_T::cmdPrint(       "-init_index:",         initial_index);
  SYS_T::cmdPrint(       "-fina_time:",          final_time);
  SYS_T::cmdPrint(       "-ttan_freq:",          ttan_renew_freq);

  if( is_record_sol )
    SYS_T::cmdPrint(     "-sol_rec_freq:",       sol_record_freq);
  else
    SYS_T::commPrint(    "-is_record_sol: false \n");

  // ====== Data for Analysis ======
  FEANode * fNode = new FEANode(part_file, rank);

  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);

  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  ALocal_NBC * locnbc = new ALocal_NBC(part_file, rank);

  ALocal_EBC * locebc = new ALocal_EBC(part_file, rank);

  APart_Node * pNode = new APart_Node_FSI(part_file, rank);

  Tissue_prestress * ps_data = new Tissue_prestress( locElem, nqp_tet, rank, is_load_ps, "prestress" );

  SYS_T::commPrint("===> Mesh HDF5 files are read from disk.\n");

  // ===== Basic Checking =====
  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for FEM analysis. \n", size);

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

  // ===== Generate the generalized-alpha method
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");

  TimeMethod_GenAlpha * tm_galpha_ptr = nullptr;

  if( is_backward_Euler )
    tm_galpha_ptr = new TimeMethod_GenAlpha( 1.0, 1.0, 1.0 );
  else
    tm_galpha_ptr = new TimeMethod_GenAlpha( genA_rho_inf, false );

  tm_galpha_ptr->print_info();

  // ===== Local assembly =====
  IMaterialModel * matmodel       = nullptr;
  IPLocAssem * locAssem_solid_ptr = nullptr;

  if( sl_nu == 0.5 )
  {
    matmodel = new MaterialModel_NeoHookean_Incompressible_Mixed( "material_model.h5" );

    locAssem_solid_ptr = new PLocAssem_Tet4_VMS_Seg_Incompressible(
        matmodel, tm_galpha_ptr, quadv->get_num_quadPts() );
  }
  else
  {
    matmodel = new MaterialModel_NeoHookean_M94_Mixed( "material_model.h5" );

    locAssem_solid_ptr = new PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha(
        matmodel, tm_galpha_ptr, quadv->get_num_quadPts() );
  }

  // ===== Initial condition =====
  PDNSolution * sol = new PDNSolution_Mixed_UPV_3D( pNode, fNode, 0 );

  PDNSolution * dot_sol = new PDNSolution_Mixed_UPV_3D( pNode, fNode, 0 );

  // Read sol file
  SYS_T::file_check(restart_name.c_str());
  sol->ReadBinary(restart_name.c_str());

  // Read dot_sol file
  std::string restart_dot_name = "dot_";
  restart_dot_name.append(restart_name);
  SYS_T::file_check(restart_dot_name.c_str());
  dot_sol->ReadBinary(restart_dot_name.c_str());

  SYS_T::commPrint("===> Read sol from disk as a restart run: \n");
  SYS_T::commPrint("     restart_name:     %s \n",     restart_name.c_str());
  SYS_T::commPrint("     restart_dot_name: %s \n", restart_dot_name.c_str());

  // ===== Time step info =====
  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);

  // ===== Global assembly routine =====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  IPGAssem * gloAssem_ptr = new PGAssem_Wall_Prestress( locAssem_solid_ptr, 
      elements, quads, GMIptr, locElem, locIEN,
      pNode, locnbc, locebc, nz_estimate );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_solid_ptr, elements, 
      quads, locIEN, pNode, locnbc, locebc );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();
  
  // ===== Linear and nonlinear solver context =====
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();

  PC upc; lsolver->GetPC(&upc);
  const PetscInt pfield[1] = {0}, vfields[] = {1,2,3};
  PCFieldSplitSetBlockSize(upc,4);
  PCFieldSplitSetFields(upc,"u",3,vfields,vfields);
  PCFieldSplitSetFields(upc,"p",1,pfield,pfield);

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

  // ===== FEM analysis =====
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
  tsolver->TM_FSI_Prestress( is_record_sol, prestress_disp_tol, dot_sol, sol, tm_galpha_ptr,
      timeinfo, locElem, locIEN, pNode, fNode, locnbc,
      locebc, pmat, elementv, elements, quadv, quads, ps_data,
      locAssem_solid_ptr, gloAssem_ptr, lsolver, nsolver);

  // ====== Finalization ======
  delete fNode; delete locIEN;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
