// ============================================================================
// wall_ps_tet4_driver.cpp
//
// Wall mechanics solver for generating the prestress field.
//
// Date: Jan 28 2022
// ============================================================================
#include "HDF5_Tools.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "APart_Node_FSI.hpp"
#include "PGAssem_Wall_Prestress.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Tet4.hpp"
#include "MaterialModel_NeoHookean_M94_Mixed.hpp"
#include "MaterialModel_NeoHookean_Incompressible_Mixed.hpp"
#include "PLocAssem_2x2Block_VMS_Incompressible.hpp"
#include "PLocAssem_2x2Block_VMS_Hyperelasticity.hpp"
#include "PTime_FSI_Solver.hpp"

int main( int argc, char *argv[] )
{
  // solution file name to be loaded for prestressing
  std::string restart_velo_name = "SOL_velo_re";
  std::string restart_pres_name = "SOL_pres_re";

  // (Pseudo-) time integration parameters
  double genA_rho_inf = 0.0;
  bool is_backward_Euler = true;
  const bool is_load_ps = false;

  // Estimate of num nonzeros per row for the sparse tangent matrix
  int nz_estimate = 300;

  // Prestress tolerance
  double prestress_disp_tol = 1.0e-6;

  // Nonlinear solver parameters
  double nl_rtol    = 1.0e-3;        // convergence criterion relative tolerance
  double nl_atol    = 1.0e-6;        // convergence criterion absolute tolerance
  double nl_dtol    = 1.0e3;         // divergence criterion
  int    nl_maxits  = 20;            // maximum number if nonlinear iterations
  int    nl_refreq  = 4;             // frequency of tangent matrix renewal
  int    nl_rethred = 4;             // threshold of tangent matrix renewal

  // Time stepping parameters
  double initial_time = 0.0;         // time of initial condition
  double initial_step = 0.1;         // time step size
  int    initial_index = 0;          // index of initial condition
  double final_time = 1.0;           // end time of simulation
  bool   is_record_sol = false;      // bool flag to decide if one wants to record the solution
  std::string sol_bName("PS_");      // base name of the solution file
  int    ttan_renew_freq = 1;        // frequency of tangent matrix renewal
  int    sol_record_freq = 1;        // frequency for recording the solution

  // Solid properties
  bool   is_read_material = true;    // bool flag to decide if one wants to read material model from h5 file
  double solid_density = 1.0;
  double solid_E = 2.0e6;
  double solid_nu = 0.5;

  // We assume that a 3D solver has been called (to generate the wall traction)
  // and a suite of command line arguments has been saved to disk
  hid_t solver_cmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( solver_cmd_file );

  const int nqp_tet  = cmd_h5r -> read_intScalar(    "/", "nqp_tet");
  const int nqp_tri  = cmd_h5r -> read_intScalar(    "/", "nqp_tri");
  const double sl_nu = cmd_h5r -> read_doubleScalar( "/", "sl_nu");

  delete cmd_h5r; H5Fclose(solver_cmd_file);

  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  HDF5_Reader * pcmd_h5r = new HDF5_Reader( prepcmd_file );

  const std::string part_v_file = pcmd_h5r -> read_string(    "/", "part_file_v" );
  const std::string part_p_file = pcmd_h5r -> read_string(    "/", "part_file_p" );
  const int fsiBC_type          = pcmd_h5r -> read_intScalar( "/", "fsiBC_type" );

  delete pcmd_h5r; H5Fclose(prepcmd_file);

  // Initialize PETSc
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  // Assert that the fsiBC type is 2, which clamped the lumen nodes
  SYS_T::print_fatal_if( fsiBC_type != 2, "Error: fsiBC_type should be 2.\n" );

  // Clean potentially pre-existing hdf5 files of prestress saved in the folder
  // named as ps_data
  if(rank == 0 )
  {
    if( SYS_T::directory_exist("ps_data") )
    {
      std::cout<<"Clean the folder ps_data.\n";
      SYS_T::execute("rm -rf ps_data");
    }

    SYS_T::execute("mkdir ps_data");
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  const std::string ps_file_name("./ps_data/prestress");

  SYS_T::GetOptionString("-restart_velo_name",   restart_velo_name);
  SYS_T::GetOptionString("-restart_pres_name",   restart_pres_name);
  SYS_T::GetOptionReal(  "-rho_inf",             genA_rho_inf);
  SYS_T::GetOptionBool(  "-is_backward_Euler",   is_backward_Euler);
  SYS_T::GetOptionInt(   "-nz_estimate",         nz_estimate);
  SYS_T::GetOptionReal(  "-prestress_disp_tol",  prestress_disp_tol);
  SYS_T::GetOptionReal(  "-nl_rtol",             nl_rtol);
  SYS_T::GetOptionReal(  "-nl_atol",             nl_atol);
  SYS_T::GetOptionReal(  "-nl_dtol",             nl_dtol);
  SYS_T::GetOptionInt(   "-nl_maxits",           nl_maxits);
  SYS_T::GetOptionInt(   "-nl_refreq",           nl_refreq);
  SYS_T::GetOptionInt(   "-nl_rethred",          nl_rethred);
  SYS_T::GetOptionBool(  "-is_backward_Euler",   is_backward_Euler);
  SYS_T::GetOptionReal(  "-init_time",           initial_time);
  SYS_T::GetOptionReal(  "-fina_time",           final_time);
  SYS_T::GetOptionReal(  "-init_step",           initial_step);
  SYS_T::GetOptionInt(   "-init_index",          initial_index);
  SYS_T::GetOptionInt(   "-ttan_freq",           ttan_renew_freq);
  SYS_T::GetOptionBool(  "-is_record_sol",       is_record_sol);
  SYS_T::GetOptionInt(   "-sol_rec_freq",        sol_record_freq);
  SYS_T::GetOptionBool(  "-is_read_material",    is_read_material);
  SYS_T::GetOptionReal(  "-sl_density",          solid_density);
  SYS_T::GetOptionReal(  "-sl_E",                solid_E);
  SYS_T::GetOptionReal(  "-sl_nu",               solid_nu);

  // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint(      "part_v_file:",          part_v_file);
  SYS_T::cmdPrint(      "part_p_file:",          part_p_file);
  SYS_T::cmdPrint(       "-prestress_disp_tol:", prestress_disp_tol);
  SYS_T::cmdPrint(       "-nl_rtol:",            nl_rtol);
  SYS_T::cmdPrint(       "-nl_atol:",            nl_atol);
  SYS_T::cmdPrint(       "-nl_dtol:",            nl_dtol);
  SYS_T::cmdPrint(       "-nl_maxits:",          nl_maxits);
  SYS_T::cmdPrint(       "-nl_refreq:",          nl_refreq);
  SYS_T::cmdPrint(       "-nl_rethred",          nl_rethred);

  if( is_backward_Euler )
    SYS_T::commPrint(    "-is_backward_Euler: true \n");
  else
    SYS_T::cmdPrint(     "-rho_inf:",            genA_rho_inf);

  SYS_T::cmdPrint(       "-init_time:",          initial_time);
  SYS_T::cmdPrint(       "-init_step:",          initial_step);
  SYS_T::cmdPrint(       "-init_index:",         initial_index);
  SYS_T::cmdPrint(       "-fina_time:",          final_time);
  SYS_T::cmdPrint(       "-ttan_freq:",          ttan_renew_freq);

  if( is_record_sol )
    SYS_T::cmdPrint(     "-sol_rec_freq:",       sol_record_freq);
  else
    SYS_T::commPrint(    "-is_record_sol: false \n");

  if( is_read_material )
  {
    SYS_T::commPrint(    "-is_read_material: true \n");
    SYS_T::file_check(   "material_model.h5" );
    SYS_T::commPrint(    "material_model.h5 found. \n");
  }
  else
  {
    SYS_T::commPrint(    "-is_read_material: false \n");
    SYS_T::cmdPrint(     "-sl_density:",         solid_density);
    SYS_T::cmdPrint(     "-sl_E:",               solid_E);
    SYS_T::cmdPrint(     "-sl_nu:",              solid_nu);
  }

  // ====== Record important parameters ======
  if(rank == 0)
  {
    hid_t cmd_file_id = H5Fcreate("wall_ps_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);
    
    cmdh5w -> write_string(      "ps_file_name",       ps_file_name);
    cmdh5w -> write_doubleScalar("prestress_disp_tol", prestress_disp_tol );

    delete cmdh5w; H5Fclose(cmd_file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ====== Data for Analysis ======
  FEANode * fNode = new FEANode(part_v_file, rank);

  ALocal_IEN * locIEN_v = new ALocal_IEN(part_v_file, rank);

  ALocal_IEN * locIEN_p = new ALocal_IEN(part_p_file, rank);

  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_v_file, rank);

  ALocal_Elem * locElem = new ALocal_Elem(part_v_file, rank);

  APart_Node * pNode_v = new APart_Node_FSI(part_v_file, rank);

  APart_Node * pNode_p = new APart_Node_FSI(part_p_file, rank);

  ALocal_EBC * locebc_v = new ALocal_EBC(part_v_file, rank);

  ALocal_EBC * locebc_p = new ALocal_EBC( part_p_file, rank );

  ALocal_NBC * locnbc_v = new ALocal_NBC(part_v_file, rank, "/nbc/MF");

  ALocal_NBC * locnbc_p = new ALocal_NBC(part_p_file, rank, "/nbc/MF");

  Tissue_prestress * ps_data = new Tissue_prestress(locElem, nqp_tet, rank, is_load_ps, ps_file_name);  
  SYS_T::commPrint("===> Mesh HDF5 files are read from disk.\n");

  // Group APart_Node and ALocal_NBC into a vector
  std::vector<APart_Node *> pNode_list { pNode_v, pNode_p };

  std::vector<ALocal_NBC *> locnbc_list { locnbc_v, locnbc_p };

  std::vector<APart_Node *> pNode_m_list { pNode_v };

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

  // ===== Generate the IS for pres and velo =====
  const int idx_v_start = HDF5_T::read_intScalar( SYS_T::gen_partfile_name(part_v_file, rank).c_str(), "/DOF_mapper", "start_idx" );
  const int idx_p_start = HDF5_T::read_intScalar( SYS_T::gen_partfile_name(part_p_file, rank).c_str(), "/DOF_mapper", "start_idx" );

  const int idx_v_len = pNode_v->get_dof() * pNode_v -> get_nlocalnode();
  const int idx_p_len = pNode_p->get_dof() * pNode_p -> get_nlocalnode();

  PetscInt * is_array_velo = new PetscInt[ idx_v_len ];
  for(int ii=0; ii<idx_v_len; ++ii) is_array_velo[ii] = idx_v_start + ii;

  PetscInt * is_array_pres = new PetscInt[ idx_p_len ];
  for(int ii=0; ii<idx_p_len; ++ii) is_array_pres[ii] = idx_p_start + ii;

  IS is_velo, is_pres;
  ISCreateGeneral(PETSC_COMM_WORLD, idx_v_len, is_array_velo, PETSC_COPY_VALUES, &is_velo);
  ISCreateGeneral(PETSC_COMM_WORLD, idx_p_len, is_array_pres, PETSC_COPY_VALUES, &is_pres);

  delete [] is_array_velo; is_array_velo = nullptr;
  delete [] is_array_pres; is_array_pres = nullptr;
  // ================================================================

  // ===== Generate a sparse matrix for strong enforcement of essential BC
  std::vector<int> start_idx{ idx_v_start, idx_p_start };

  Matrix_PETSc * pmat = new Matrix_PETSc( idx_v_len + idx_p_len );
  pmat -> gen_perm_bc( pNode_list, locnbc_list, start_idx );

  // ===== Generate the generalized-alpha method
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");

  TimeMethod_GenAlpha * tm_galpha_ptr = nullptr;

  if( is_backward_Euler )
    tm_galpha_ptr = new TimeMethod_GenAlpha( 1.0, 1.0, 1.0 );
  else
    tm_galpha_ptr = new TimeMethod_GenAlpha( genA_rho_inf, false );

  tm_galpha_ptr->print_info();

  // ===== Local assembly =====
  IMaterialModel * matmodel = nullptr;
  IPLocAssem_2x2Block * locAssem_solid_ptr = nullptr;

  if( is_read_material )
  {
    if( sl_nu == 0.5 )
    {
      matmodel = new MaterialModel_NeoHookean_Incompressible_Mixed( "material_model.h5" );

      locAssem_solid_ptr = new PLocAssem_2x2Block_VMS_Incompressible(
          matmodel, tm_galpha_ptr, elementv -> get_nLocBas(), elements->get_nLocBas() );
    }
    else
    {
      matmodel = new MaterialModel_NeoHookean_M94_Mixed( "material_model.h5" );

      locAssem_solid_ptr = new PLocAssem_2x2Block_VMS_Hyperelasticity(
          matmodel, tm_galpha_ptr, elementv -> get_nLocBas(), elements->get_nLocBas() );
    }
  }
  else
  {
    if( solid_nu == 0.5 )
    {
      matmodel = new MaterialModel_NeoHookean_Incompressible_Mixed( solid_density, solid_E );

      locAssem_solid_ptr = new PLocAssem_2x2Block_VMS_Incompressible(
          matmodel, tm_galpha_ptr, elementv -> get_nLocBas(), elements->get_nLocBas() );
    }
    else
    {
      matmodel = new MaterialModel_NeoHookean_M94_Mixed( solid_density, solid_E, solid_nu );

      locAssem_solid_ptr = new PLocAssem_2x2Block_VMS_Hyperelasticity(
          matmodel, tm_galpha_ptr, elementv -> get_nLocBas(), elements->get_nLocBas() );
    }
  }

  // ===== Initial conditions =====
  PDNSolution * velo = new PDNSolution_V(pNode_v, 0, true, "velo");
  PDNSolution * disp = new PDNSolution_V(pNode_v, 0, true, "disp");
  PDNSolution * pres = new PDNSolution_P(pNode_p, 0, true, "pres");

  PDNSolution * dot_velo = new PDNSolution_V(pNode_v, 0, true, "dot_velo");
  PDNSolution * dot_disp = new PDNSolution_V(pNode_v, 0, true, "dot_disp");
  PDNSolution * dot_pres = new PDNSolution_P(pNode_p, 0, true, "dot_pres");

  // Read sol file
  SYS_T::file_check(restart_velo_name.c_str());
  velo->ReadBinary(restart_velo_name.c_str());

  SYS_T::file_check(restart_pres_name.c_str());
  pres->ReadBinary(restart_pres_name.c_str());

  // Read dot_sol file
  std::string restart_dot_velo_name = "dot_";
  restart_dot_velo_name.append(restart_velo_name);
  SYS_T::file_check(restart_dot_velo_name.c_str());
  dot_velo->ReadBinary(restart_dot_velo_name.c_str());
 
  std::string restart_dot_pres_name = "dot_";
  restart_dot_pres_name.append(restart_pres_name);
  SYS_T::file_check(restart_dot_pres_name.c_str());
  dot_pres->ReadBinary(restart_dot_pres_name.c_str());

  SYS_T::commPrint("===> Read sol from disk as a restart run: \n");
  SYS_T::commPrint("     restart_velo_name:     %s \n", restart_velo_name.c_str());
  SYS_T::commPrint("     restart_dot_velo_name: %s \n", restart_dot_velo_name.c_str());
  SYS_T::commPrint("     restart_pres_name:     %s \n", restart_pres_name.c_str());
  SYS_T::commPrint("     restart_dot_pres_name: %s \n", restart_dot_pres_name.c_str());

  // ===== Time step info =====
  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);

  // ===== Global assembly routine =====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  IPGAssem * gloAssem_ptr = new PGAssem_Wall_Prestress( locAssem_solid_ptr, locElem, 
      locIEN_v, locIEN_p, pNode_v, pNode_p, locnbc_v, locnbc_p, locebc_v, nz_estimate );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_solid_ptr, locIEN_v, locIEN_p, 
      locnbc_v, locnbc_p );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();

  // ===== Linear and nonlinear solver context =====
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();

  PC upc; lsolver->GetPC(&upc);
  PCFieldSplitSetIS(upc, "u", is_velo);
  PCFieldSplitSetIS(upc, "p", is_pres);

  // ===== Nonlinear solver context =====
  PNonlinear_FSI_Solver * nsolver = new PNonlinear_FSI_Solver(
      nl_rtol, nl_atol, nl_dtol, nl_maxits, nl_refreq, nl_rethred);
  SYS_T::commPrint("===> Nonlinear solver setted up:\n");
  nsolver->print_info();

  // ===== Temporal solver context =====
  PTime_FSI_Solver * tsolver = new PTime_FSI_Solver( sol_bName,
      sol_record_freq, ttan_renew_freq, final_time );
  SYS_T::commPrint("===> Time marching solver setted up:\n");
  tsolver->print_info();

  // ===== FEM analysis =====
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
  tsolver -> TM_FSI_Prestress( is_record_sol, prestress_disp_tol, is_velo, is_pres,
      dot_disp, dot_velo, dot_pres, disp, velo, pres, tm_galpha_ptr,
      timeinfo, locElem, locIEN_v, locIEN_p, pNode_v, pNode_p, fNode,
      locnbc_v, locnbc_p, locebc_v, locebc_p, pmat, elementv, elements,
      quadv, quads, ps_data, locAssem_solid_ptr, gloAssem_ptr, lsolver, nsolver );

  // ===== Record the wall prestress to h5 file =====
  ps_data -> write_prestress_hdf5();

  // ==========================================================================
  // Clean the memory
  delete fNode; delete locIEN_v; delete locIEN_p; delete PartBasic; delete locElem;
  delete pNode_v; delete pNode_p; delete locebc_v; delete locebc_p; 
  delete locnbc_v; delete locnbc_p;
  delete ps_data; delete quadv; delete quads; delete elementv; delete elements;
  delete pmat; delete tm_galpha_ptr; delete matmodel; delete locAssem_solid_ptr;
  delete velo; delete disp; delete pres; delete dot_velo; delete dot_disp; delete dot_pres;
  delete timeinfo; delete gloAssem_ptr; delete lsolver; delete nsolver; delete tsolver;
  ISDestroy(&is_velo); ISDestroy(&is_pres);
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF 
