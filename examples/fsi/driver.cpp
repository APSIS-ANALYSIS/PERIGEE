// ============================================================================
// fsi_driver.cpp
// 
// Author: Ju Liu
// Date: Dec. 26 2021
// ============================================================================
#include "HDF5_Tools.hpp"
#include "AGlobal_Mesh_Info.hpp"
#include "ANL_Tools.hpp"
#include "APart_Node_FSI.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss_Hex.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Hex8.hpp"
#include "FlowRateFactory.hpp"
#include "GenBC_Resistance.hpp"
#include "GenBC_RCR.hpp"
#include "GenBC_Inductance.hpp"
#include "GenBC_Coronary.hpp"
#include "GenBC_Pressure.hpp"
#include "MaterialModel_vol_Incompressible.hpp"
#include "MaterialModel_vol_ST91.hpp"
#include "MaterialModel_vol_M94.hpp"
#include "MaterialModel_ich_NeoHookean.hpp"
#include "MaterialModel_ich_GOH14.hpp"
#include "PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha.hpp"
#include "PLocAssem_2x2Block_VMS_Incompressible.hpp"
#include "PLocAssem_2x2Block_VMS_Hyperelasticity.hpp"
#include "PLocAssem_FSI_Mesh_Elastostatic.hpp"
#include "PLocAssem_FSI_Mesh_Laplacian.hpp"
#include "PGAssem_FSI.hpp"
#include "PGAssem_Mesh.hpp"
#include "PTime_FSI_Solver.hpp"

int main(int argc, char *argv[])
{
  // Number of quadrature points 5/4 for tet and triangle
  // Number of quadrature points 8/4  for hex and quadrangle
  int nqp_vol = 5, nqp_sur = 4;

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

  // inflow file
  std::string inflow_file("inflow_fourier_series.txt");

  // LPN file
  std::string lpn_file("lpn_rcr_input.txt");

  // back flow stabilization
  double bs_beta = 0.2;

  // flag that determines if the prestress data to be loaded
  bool is_load_ps = false;

  // Generalized-alpha method
  double genA_rho_inf = 0.5;
  bool is_backward_Euler = false;

  // part file location
  const std::string part_v_file("./apart/part_v");
  const std::string part_p_file("./apart/part_p");

  // Nonlinear solver parameters
  double nl_rtol = 1.0e-3;
  double nl_atol = 1.0e-6;
  double nl_dtol = 10.0;
  int nl_maxits  = 20;
  int nl_refreq  = 4;
  int nl_rethred = 4;

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
  std::string restart_u_name = "SOL_U_";
  std::string restart_v_name = "SOL_V_";
  std::string restart_p_name = "SOL_P_";

  // Yaml options
  bool is_loadYaml = true;
  std::string yaml_file("./runscript.yml");

  // ===== Initialization of PETSc =====
#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::print_perigee_art();

  SYS_T::commPrint("Job started on %s %s \n", SYS_T::get_time().c_str(), SYS_T::get_date().c_str());
  SYS_T::commPrint("PETSc version: %s \n", PETSc_T::get_version().c_str());

  // ===== Yaml Argument =====
  SYS_T::GetOptionBool(  "-is_loadYaml",       is_loadYaml);
  SYS_T::GetOptionString("-yaml_file",         yaml_file);

  // load the YAML file to pass the argument values
  if(is_loadYaml) SYS_T::InsertFileYAML( yaml_file,  false );

  // ===== Command Line Argument =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");

  SYS_T::GetOptionInt(   "-nqp_vol",           nqp_vol);
  SYS_T::GetOptionInt(   "-nqp_sur",           nqp_sur);;
  SYS_T::GetOptionInt(   "-nz_estimate",       nz_estimate);
  SYS_T::GetOptionReal(  "-bs_beta",           bs_beta);
  SYS_T::GetOptionReal(  "-rho_inf",           genA_rho_inf);
  SYS_T::GetOptionBool(  "-is_backward_Euler", is_backward_Euler);
  SYS_T::GetOptionBool(  "-is_load_ps",        is_load_ps);
  SYS_T::GetOptionReal(  "-fl_density",        fluid_density);
  SYS_T::GetOptionReal(  "-fl_mu",             fluid_mu);
  SYS_T::GetOptionReal(  "-sl_density",        solid_density);
  SYS_T::GetOptionReal(  "-sl_E",              solid_E);
  SYS_T::GetOptionReal(  "-sl_nu",             solid_nu);
  SYS_T::GetOptionReal(  "-mesh_E",            mesh_E);
  SYS_T::GetOptionReal(  "-mesh_nu",           mesh_nu);
  SYS_T::GetOptionString("-inflow_file",       inflow_file);
  SYS_T::GetOptionString("-lpn_file",          lpn_file);
  SYS_T::GetOptionReal(  "-nl_rtol",           nl_rtol);
  SYS_T::GetOptionReal(  "-nl_atol",           nl_atol);
  SYS_T::GetOptionReal(  "-nl_dtol",           nl_dtol);
  SYS_T::GetOptionInt(   "-nl_maxits",         nl_maxits);
  SYS_T::GetOptionInt(   "-nl_refreq",         nl_refreq);
  SYS_T::GetOptionInt(   "-nl_rethred",        nl_rethred);
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
  SYS_T::GetOptionString("-restart_u_name",    restart_u_name);
  SYS_T::GetOptionString("-restart_v_name",    restart_v_name);
  SYS_T::GetOptionString("-restart_p_name",    restart_p_name);

  // ===== Print the command line argumetn on screen =====
  SYS_T::cmdPrint("-nqp_vol:", nqp_vol);
  SYS_T::cmdPrint("-nqp_sur:", nqp_sur);
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

  if( is_load_ps )
    SYS_T::commPrint(     "-is_load_ps: true \n");
  else
    SYS_T::commPrint(     "-is_load_ps: false \n");

  SYS_T::cmdPrint("-inflow_file:", inflow_file);
  SYS_T::cmdPrint("-lpn_file:", lpn_file);
  SYS_T::cmdPrint("-nl_rtol:", nl_rtol);
  SYS_T::cmdPrint("-nl_atol:", nl_atol);
  SYS_T::cmdPrint("-nl_dtol:", nl_dtol);
  SYS_T::cmdPrint("-nl_maxits:", nl_maxits);
  SYS_T::cmdPrint("-nl_refreq:", nl_refreq);
  SYS_T::cmdPrint("-nl_rethred", nl_rethred);
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
    SYS_T::cmdPrint("-restart_u_name:", restart_u_name);
    SYS_T::cmdPrint("-restart_v_name:", restart_v_name);
    SYS_T::cmdPrint("-restart_p_name:", restart_p_name);
  }
  else SYS_T::commPrint("-is_restart: false \n");

  // ===== Record important parameters =====
  if(rank == 0)
  {
    hid_t cmd_file_id = H5Fcreate("solver_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

    cmdh5w->write_doubleScalar(  "fl_density",      fluid_density);
    cmdh5w->write_doubleScalar(  "fl_mu",           fluid_mu);
    cmdh5w->write_doubleScalar(  "fl_bs_beta",      bs_beta);
    cmdh5w->write_doubleScalar(  "sl_density",      solid_density);
    cmdh5w->write_doubleScalar(  "sl_E",            solid_E);
    cmdh5w->write_doubleScalar(  "sl_nu",           solid_nu);
    cmdh5w->write_doubleScalar(  "mesh_E",          mesh_E);
    cmdh5w->write_doubleScalar(  "mesh_nu",         mesh_nu);
    cmdh5w->write_doubleScalar(  "init_step",       initial_step);
    cmdh5w->write_intScalar(     "sol_record_freq", sol_record_freq);
    cmdh5w->write_intScalar(     "nqp_vol",         nqp_vol);
    cmdh5w->write_intScalar(     "nqp_sur",         nqp_sur);
    cmdh5w->write_string(        "lpn_file",        lpn_file);
    cmdh5w->write_string(        "inflow_file",     inflow_file);
    cmdh5w->write_string(        "sol_bName",       sol_bName);

    cmdh5w->write_string("date",              SYS_T::get_date() );
    cmdh5w->write_string("time",              SYS_T::get_time() );
    cmdh5w->write_string("petsc-version",     PETSc_T::get_version() );

    delete cmdh5w; H5Fclose(cmd_file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
    
  // ===== Main Data Strucutre =====
  // Control points are only stored for the geometry-defining field, that is the velo/disp
  // field.
  auto locElem = SYS_T::make_unique<ALocal_Elem>(part_v_file, rank);
  
  auto locIEN_v = SYS_T::make_unique<ALocal_IEN>(part_v_file, rank);

  auto locIEN_p = SYS_T::make_unique<ALocal_IEN>(part_p_file, rank);

  auto fNode = SYS_T::make_unique<FEANode>(part_v_file, rank);

  auto pNode_v = SYS_T::make_unique<APart_Node>(part_v_file, rank);

  auto pNode_p = SYS_T::make_unique<APart_Node>(part_p_file, rank);

  auto locinfnbc = SYS_T::make_unique<ALocal_InflowBC>(part_v_file, rank);

  auto locnbc_v = SYS_T::make_unique<ALocal_NBC>(part_v_file, rank, "/nbc/MF");

  auto locnbc_p = SYS_T::make_unique<ALocal_NBC>(part_v_file, rank, "/nbc/MF");

  auto mesh_locnbc = SYS_T::make_unique<ALocal_NBC>(part_v_file, rank, "/mesh_nbc/MF");

  std::unique_ptr<ALocal_EBC> locebc_v = SYS_T::make_unique<ALocal_EBC_outflow>(part_v_file, rank);

  std::unique_ptr<ALocal_EBC> locebc_p = SYS_T::make_unique<ALocal_EBC>(part_p_file, rank);

  std::unique_ptr<ALocal_EBC> mesh_locebc = SYS_T::make_unique<ALocal_EBC>(part_v_file, rank, "/mesh_ebc");

  auto ps_data = SYS_T::make_unique<Tissue_prestress>(locElem.get(), nqp_vol, rank, is_load_ps, "./ps_data/prestress");

  // Group APart_Node and ALocal_NBC into a vector
  std::vector<std::unique_ptr<APart_Node>> pNode_list { std::move(pNode_v), std::move(pNode_p) };

  std::vector<std::unique_ptr<ALocal_NBC>> locnbc_list { std::move(locnbc_v), std::move(locnbc_p) };

  std::vector<std::unique_ptr<APart_Node>> pNode_m_list { std::move(pNode_v) };

  std::vector<std::unique_ptr<ALocal_NBC>> locnbc_m_list { std::move(mesh_locnbc) };

  // ===== Basic Checking =====
  SYS_T::print_fatal_if( size!= ANL_T::get_cpu_size(part_v_file, rank),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for FEM analysis. \n", size);

  // ===== Inflow rate function =====
  SYS_T::commPrint("===> Setup inflow flow rate. \n");

  auto inflow_rate = FlowRateFactory::createFlowRate(inflow_file);

  inflow_rate->print_info();

  SYS_T::print_fatal_if(locinfnbc->get_num_nbc() != inflow_rate->get_num_nbc(),
      "Error: ALocal_InflowBC number of faces does not match with that in IFlowRate.\n");

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

  const int idx_m_start = pNode_v->get_node_loc(0) * locnbc_v->get_dof_LID();
  std::vector<int> start_m_idx{ idx_m_start };
  
  Matrix_PETSc * mmat = new Matrix_PETSc( pNode_v -> get_nlocalnode() * pNode_v -> get_dof() );
  mmat -> gen_perm_bc( pNode_m_list, locnbc_m_list, start_m_idx );

  // ===== Generate the generalized-alpha method
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");

  TimeMethod_GenAlpha * tm_galpha_ptr = nullptr;

  if( is_backward_Euler )
    tm_galpha_ptr = new TimeMethod_GenAlpha( 1.0, 1.0, 1.0 );
  else
    tm_galpha_ptr = new TimeMethod_GenAlpha( genA_rho_inf, false );

  tm_galpha_ptr->print_info();

  // ===== Local assembly =====
  IPLocAssem_2x2Block * locAssem_fluid_ptr = new PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha(
      ANL_T::get_elemType(part_file, rank), nqp_vol, nqp_sur,
      tm_galpha.get(), fluid_density, fluid_mu, bs_beta );  

  IPLocAssem_2x2Block * locAssem_solid_ptr = nullptr;

  const double solid_mu = solid_E/(2.0+2.0*solid_nu);
  std::unique_ptr<IMaterialModel_ich> imodel = SYS_T::make_unique<MaterialModel_ich_NeoHookean>(solid_mu);

  if( solid_nu == 0.5 )
  {
    std::unique_ptr<IMaterialModel_vol> vmodel = SYS_T::make_unique<MaterialModel_vol_Incompressible>(solid_density);
    std::unique_ptr<MaterialModel_Mixed_Elasticity> matmodel = SYS_T::make_unique<MaterialModel_Mixed_Elasticity>(std::move(vmodel), std::move(imodel));
    locAssem_solid_ptr = new PLocAssem_2x2Block_VMS_Incompressible(
        ANL_T::get_elemType(part_file, rank), nqp_vol, nqp_sur, tm_galpha_ptr, std::move(matmodel) );
  }
  else
  {
    const double solid_lambda = solid_nu * solid_E / ((1+solid_nu) * (1-2.0*solid_nu));
    const double solid_kappa  = solid_lambda + 2.0 * solid_mu / 3.0;
    std::unique_ptr<IMaterialModel_vol> vmodel = SYS_T::make_unique<MaterialModel_vol_M94>(solid_density, solid_kappa);
    std::unique_ptr<MaterialModel_Mixed_Elasticity> matmodel = SYS_T::make_unique<MaterialModel_Mixed_Elasticity>(std::move(vmodel), std::move(imodel));
    locAssem_solid_ptr = new PLocAssem_2x2Block_VMS_Hyperelasticity(
        ANL_T::get_elemType(part_file, rank), nqp_vol, nqp_sur, tm_galpha_ptr, std::move(matmodel) );
  }

  // Pseudo elastic mesh motion
  IPLocAssem * locAssem_mesh_ptr = new PLocAssem_FSI_Mesh_Laplacian( ANL_T::get_elemType(part_file, rank), nqp_vol, nqp_sur );
  
  // ===== Initial condition =====
  PDNSolution * base = new PDNSolution_V( pNode_v, fNode, locinfnbc, 1, true, "base" ); 
  
  PDNSolution * velo = new PDNSolution_V(pNode_v, 0, true, "velo");
  PDNSolution * disp = new PDNSolution_V(pNode_v, 0, true, "disp");
  PDNSolution * pres = new PDNSolution_P(pNode_p, 0, true, "pres");

  PDNSolution * dot_velo = new PDNSolution_V(pNode_v, 0, true, "dot_velo");
  PDNSolution * dot_disp = new PDNSolution_V(pNode_v, 0, true, "dot_disp");
  PDNSolution * dot_pres = new PDNSolution_P(pNode_p, 0, true, "dot_pres");
  
  if( is_restart )
  {
    initial_index = restart_index;
    initial_time  = restart_time;
    initial_step  = restart_step;

    // Read sol file
    SYS_T::file_check(restart_u_name);
    SYS_T::file_check(restart_v_name);
    SYS_T::file_check(restart_p_name);
    
    disp->ReadBinary(restart_u_name);
    velo->ReadBinary(restart_v_name);
    pres->ReadBinary(restart_p_name);
    
    // Read dot_sol file
    std::string restart_dot_u_name = "dot_";
    restart_dot_u_name.append(restart_u_name);
    SYS_T::file_check(restart_dot_u_name);
    dot_disp->ReadBinary(restart_dot_u_name);

    std::string restart_dot_v_name = "dot_";
    restart_dot_v_name.append(restart_v_name);
    SYS_T::file_check(restart_dot_v_name);
    dot_velo->ReadBinary(restart_dot_v_name);

    std::string restart_dot_p_name = "dot_";
    restart_dot_p_name.append(restart_p_name);
    SYS_T::file_check(restart_dot_p_name);
    dot_pres->ReadBinary(restart_dot_p_name);

    SYS_T::commPrint("===> Read sol from disk as a restart run... \n");
    SYS_T::commPrint("     restart_u_name: %s \n", restart_u_name.c_str());
    SYS_T::commPrint("     restart_v_name: %s \n", restart_v_name.c_str());
    SYS_T::commPrint("     restart_p_name: %s \n", restart_p_name.c_str());
    SYS_T::commPrint("     restart_dot_u_name: %s \n", restart_dot_u_name.c_str());
    SYS_T::commPrint("     restart_dot_v_name: %s \n", restart_dot_v_name.c_str());
    SYS_T::commPrint("     restart_dot_p_name: %s \n", restart_dot_p_name.c_str());
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

  SYS_T::print_fatal_if(gbc->get_num_ebc() != locebc_v->get_num_ebc(),
      "Error: GenBC number of faces does not match with that in ALocal_EBC.\n");

  // ===== Global assembly routine =====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  IPGAssem * gloAssem_ptr = new PGAssem_FSI( 
      gbc.get(), std::move(locIEN_v), std::move(locIEN_p), std::move(locElem), 
      std::move(fNode), std::move(pNode_v), std::move(pNode_p), std::move(locnbc_v), 
      std::move(locnbc_p), std::move(locebc_v), std::move(locebc_p), 
      std::move(locAssem_fluid_ptr), std::move(locAssem_solid_ptr), nz_estimate );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem_ptr->Assem_nonzero_estimate( gbc );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();

  // ===== Global assembly for mesh motion =====
  SYS_T::commPrint("===> Initializing Mat K_mesh and Vec G_mesh ... \n");
  IPGAssem * gloAssem_mesh_ptr = new PGAssem_Mesh( locAssem_mesh_ptr,
      locElem, locIEN_v, pNode_v, mesh_locnbc, mesh_locebc, nz_estimate );

  SYS_T::commPrint("===> Assembly nonzero estimate for K_mesh ... \n");
  gloAssem_mesh_ptr->Assem_nonzero_estimate( locElem, locAssem_mesh_ptr, locIEN_v, mesh_locnbc );

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
    KSPGMRESSetOrthogonalization(lsolver_acce->ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);
    KSPGMRESSetRestart(lsolver_acce->ksp, 500);

    lsolver_acce->Monitor();

    PC preproc; lsolver_acce->GetPC(&preproc);
    PCSetType( preproc, PCHYPRE );
    PCHYPRESetType( preproc, "boomeramg" );

    gloAssem_ptr->Assem_mass_residual( disp, velo, pres,  ps_data );

    Vec proj_vp, proj_v, proj_p;
    VecDuplicate( gloAssem_ptr->G, &proj_vp );
    lsolver_acce -> Solve( gloAssem_ptr->K, gloAssem_ptr->G, proj_vp );

    SYS_T::commPrint("\n===> Consistent initial acceleration is obtained.\n");
    lsolver_acce -> print_info();

    VecGetSubVector(proj_vp, is_velo, &proj_v);
    VecGetSubVector(proj_vp, is_pres, &proj_p);

    // [dot_v, dot_p] = -M^-1 [Res_v, Res_p]
    dot_velo -> PlusAX(proj_v, -1.0);
    dot_pres -> PlusAX(proj_p, -1.0);

    // dot_u = v
    dot_disp -> Copy( velo ); 

    VecRestoreSubVector(proj_vp, is_velo, &proj_v);
    VecRestoreSubVector(proj_vp, is_pres, &proj_p);
    VecDestroy(&proj_vp);
    delete lsolver_acce;
    SYS_T::commPrint("===> The mass matrix lsolver is destroyed.\n");
  }

  // ===== Linear and nonlinear solver context =====
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();

  PC upc; lsolver->GetPC(&upc);
  PCFieldSplitSetIS(upc, "u", is_velo);
  PCFieldSplitSetIS(upc, "p", is_pres);

  PLinear_Solver_PETSc * mesh_lsolver = new PLinear_Solver_PETSc(
      1.0e-12, 1.0e-55, 1.0e30, 500, "mesh_", "mesh_" );

  gloAssem_mesh_ptr->Assem_tangent_residual( disp, disp, 0.0,
      timeinfo->get_step(), locElem, locAssem_mesh_ptr, elementv,
      elements, quadv, quads, locIEN_v, fNode, mesh_locnbc,
      mesh_locebc );

  mesh_lsolver -> SetOperator( gloAssem_mesh_ptr->K );
  PC mesh_pc; mesh_lsolver->GetPC(&mesh_pc);
  PCFieldSplitSetBlockSize( mesh_pc, 3 );

  mesh_lsolver -> print_info();

  SYS_T::commPrint("===> mesh solver LHS setted up.\n");

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

  // ===== Outlet flowrate recording files =====
  for(int ff=0; ff<locebc_v->get_num_ebc(); ++ff)
  {
    const double dot_face_flrate = gloAssem_ptr -> Assem_surface_flowrate( disp, dot_velo, ff );

    const double face_flrate = gloAssem_ptr -> Assem_surface_flowrate( disp, velo, ff );

    const double face_avepre = gloAssem_ptr -> Assem_surface_ave_pressure( disp, pres, ff );

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
        ofile.open( tsolver->gen_flowfile_name("Outlet_", ff).c_str(), std::ofstream::out | std::ofstream::trunc );
      else
        ofile.open( tsolver->gen_flowfile_name("Outlet_", ff).c_str(), std::ofstream::out | std::ofstream::app );

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
    const double inlet_face_flrate = gloAssem_ptr -> Assem_surface_flowrate( disp, velo, locinfnbc, ff );

    const double inlet_face_avepre = gloAssem_ptr -> Assem_surface_ave_pressure( disp, pres, locinfnbc, ff );

    if( rank == 0 )
    {
      std::ofstream ofile;
      if( !is_restart )
        ofile.open( tsolver->gen_flowfile_name("Inlet_", ff).c_str(), std::ofstream::out | std::ofstream::trunc );
      else
        ofile.open( tsolver->gen_flowfile_name("Inlet_", ff).c_str(), std::ofstream::out | std::ofstream::app );

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
#ifdef PETSC_USE_LOG
  PetscLogEvent tsolver_event;
  PetscClassId classid;
  PetscClassIdRegister("user-log-info-time-solver", &classid);
  PetscLogEventRegister("time_solver", classid, &tsolver_event);
#endif

  SYS_T::commPrint("===> Start Finite Element Analysis:\n");

#ifdef PETSC_USE_LOG
  PetscLogEventBegin(tsolver_event, 0,0,0,0);
#endif

  tsolver->TM_FSI_GenAlpha(is_restart, is_velo, is_pres, base, 
      dot_disp, dot_velo, dot_pres, disp, velo, pres, 
      tm_galpha_ptr, timeinfo, inflow_rate.get(), locElem, locIEN_v, locIEN_p, 
      pNode_v, pNode_p, fNode, locnbc_v, locnbc_p, locinfnbc, mesh_locnbc, 
      locebc_v, locebc_p, mesh_locebc, 
      gbc, pmat, mmat, elementv, elements, quadv, quads, ps_data,
      locAssem_fluid_ptr, locAssem_solid_ptr, locAssem_mesh_ptr,
      gloAssem_ptr, gloAssem_mesh_ptr, lsolver, mesh_lsolver, nsolver);

#ifdef PETSC_USE_LOG
  PetscLogEventEnd(tsolver_event,0,0,0,0);
#endif

  // Print complete solver info
  lsolver -> print_info();
  mesh_lsolver -> print_info();

  // ===== PETSc Finalize =====
  delete tsolver; delete nsolver; delete lsolver; delete mesh_lsolver;
  delete gloAssem_ptr; delete gloAssem_mesh_ptr;
  delete timeinfo; delete gbc;
  delete pres; delete dot_pres;
  delete base; delete dot_velo; delete dot_disp; delete velo; delete disp;
  delete locAssem_mesh_ptr; delete locAssem_fluid_ptr;
  delete locAssem_solid_ptr; delete pmat; delete mmat; delete tm_galpha_ptr;
  ISDestroy(&is_velo); ISDestroy(&is_pres);
  delete elements; delete elementv; delete quadv; delete quads;
  delete GMIptr; delete locElem; delete fNode; delete pNode_v; delete pNode_p;
  delete locinfnbc; delete locnbc_v; delete locnbc_p; delete mesh_locnbc; 
  delete locebc_v; delete locebc_p; delete mesh_locebc; 
  delete locIEN_v; delete locIEN_p; delete ps_data;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
