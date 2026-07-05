// ==================================================================
// ns code driver.cpp
//
// Finite element code for 3D Navier-Stokes equations using the 
// Variational Multiscale Formulation and Generalized-alpha time 
// stepping.
//
// Author: Ju Liu, liujuy@gmail.com
// Date: Feb. 6 2020
// ==================================================================
#include "HDF5_Writer.hpp"
#include "ANL_Tools.hpp"
#include "FlowRateFactory.hpp"
#include "GenBCFactory.hpp"
#include "HIModelFactory.hpp"
#include "HITools.hpp"
#include "InitHelpers.hpp"
#include "PLocAssem_Transport_VMS_NS_GenAlpha.hpp"
#include "PLocAssem_VMS_NS_GenAlpha_WeakBC.hpp"
#include "PGAssem_Transport_NS_FEM.hpp"
#include "PGAssem_NS_FEM.hpp"
#include "PTime_NS_Solver.hpp"

int main(int argc, char *argv[])
{
  // Coefficient for weak bc
  double C_bI = 4.0;

  // Number of quadrature points for tets and triangles
  // Suggested values: 5 / 4 for linear, 17 / 13 for quadratic
  // Number of quadrature points for hexs and quadrangles
  // Suggested values: 8 / 4 for linear, 64 / 16 for quadratic
  int nqp_vol = 5, nqp_sur = 4;

  // Estimate of the nonzero per row for the sparse matrix
  int nz_estimate = 300;

  // fluid properties
  double fluid_density = 1.065;
  double fluid_mu = 3.5e-2;
  double c_tauc = 1.0; // scaling factor for tau_c, take 0.0, 0.125, or 1.0
  double c_ct = 4.0;   // C_T parameter for defining tau_M

  // inflow file
  std::string inflow_file("inflow_fourier_series.txt");

  // LPN file
  std::string lpn_file("lpn_rcr_input.txt");

  // back flow stabilization
  double bs_beta = 0.2;

  // generalized-alpha rho_inf
  double genA_rho_inf = 0.5;
  bool is_backward_Euler = false;

  // part file location
  std::string part_file("part");

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

  // Optional weakly coupled scalar transport solve
  bool is_transport = false;
  double transport_ct = 1.0;
  std::string transport_name("TRNS_");
  std::string hi_model_name("none");

  // Restart options
  bool is_restart = false;
  int restart_index = 0;             // restart solution time index
  double restart_time = 0.0;         // restart time
  double restart_step = 1.0e-3;      // restart simulation time step size
  std::string restart_name = "SOL_"; // restart solution base name

  // Yaml options
  bool is_loadYaml = true;
  std::string yaml_file("./runscript.yml");

#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::print_perigee_art();

  // ===== Yaml Arguments =====
  SYS_T::GetOptionBool("-is_loadYaml", is_loadYaml);
  SYS_T::GetOptionString("-yaml_file", yaml_file);

  if (is_loadYaml) SYS_T::InsertFileYAML( yaml_file,  false );

  // ===== Read Command Line Arguments =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");

  SYS_T::GetOptionInt("-nqp_vol", nqp_vol);
  SYS_T::GetOptionInt("-nqp_sur", nqp_sur);
  SYS_T::GetOptionInt("-nz_estimate", nz_estimate);
  SYS_T::GetOptionReal("-bs_beta", bs_beta);
  SYS_T::GetOptionReal("-rho_inf", genA_rho_inf);
  SYS_T::GetOptionBool("-is_backward_Euler", is_backward_Euler);
  SYS_T::GetOptionReal("-fl_density", fluid_density);
  SYS_T::GetOptionReal("-fl_mu", fluid_mu);
  SYS_T::GetOptionReal("-c_tauc", c_tauc);
  SYS_T::GetOptionReal("-c_ct", c_ct);
  SYS_T::GetOptionString("-inflow_file", inflow_file);
  SYS_T::GetOptionString("-lpn_file", lpn_file);
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
  SYS_T::GetOptionBool("-is_transport", is_transport);
  SYS_T::GetOptionReal("-transport_ct", transport_ct);
  SYS_T::GetOptionString("-transport_name", transport_name);
  SYS_T::GetOptionString("-hi_model", hi_model_name);
  SYS_T::GetOptionBool("-is_restart", is_restart);
  SYS_T::GetOptionInt("-restart_index", restart_index);
  SYS_T::GetOptionReal("-restart_time", restart_time);
  SYS_T::GetOptionReal("-restart_step", restart_step);
  SYS_T::GetOptionString("-restart_name", restart_name);
  SYS_T::GetOptionReal("-C_bI", C_bI);

  // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint("-nqp_vol:", nqp_vol);
  SYS_T::cmdPrint("-nqp_sur:", nqp_sur);
  if( is_backward_Euler )
    SYS_T::commPrint(   "-is_backward_Euler: true \n");
  else
    SYS_T::cmdPrint(    "-rho_inf:",         genA_rho_inf);
  SYS_T::cmdPrint("-nz_estimate:", nz_estimate);
  SYS_T::cmdPrint("-bs_beta:", bs_beta);
  SYS_T::cmdPrint("-rho_inf:", genA_rho_inf);
  SYS_T::cmdPrint("-fl_density:", fluid_density);
  SYS_T::cmdPrint("-fl_mu:", fluid_mu);
  SYS_T::cmdPrint("-c_tauc:", c_tauc);
  SYS_T::cmdPrint("-c_ct:", c_ct);
  SYS_T::cmdPrint("-inflow_file:", inflow_file);
  SYS_T::cmdPrint("-lpn_file:", lpn_file);
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
  SYS_T::cmdPrint("-transport_ct:", transport_ct);
  SYS_T::cmdPrint("-transport_name:", transport_name);
  SYS_T::cmdPrint("-hi_model:", hi_model_name);
  if(is_transport) SYS_T::commPrint("-is_transport: true \n");
  else SYS_T::commPrint("-is_transport: false \n");
  if(is_restart)
  {
    SYS_T::commPrint("-is_restart: true \n");
    SYS_T::cmdPrint("-restart_index:", restart_index);
    SYS_T::cmdPrint("-restart_time:", restart_time);
    SYS_T::cmdPrint("-restart_step:", restart_step);
    SYS_T::cmdPrint("-restart_name:", restart_name);
  }
  else SYS_T::commPrint("-is_restart: false \n");

  std::unique_ptr<IHIModel> restart_hi_model = nullptr;
  if( is_transport ) restart_hi_model = HIModelFactory::createModel(hi_model_name);

  // ===== Record important solver options =====
  if(rank == 0)
  {
    auto cmdh5w = SYS_T::make_unique<HDF5_Writer>("solver_cmd.h5");

    cmdh5w->write_doubleScalar("fl_density", fluid_density);
    cmdh5w->write_doubleScalar("fl_mu", fluid_mu);
    cmdh5w->write_doubleScalar("init_step", initial_step);
    cmdh5w->write_intScalar("sol_record_freq", sol_record_freq);
    cmdh5w->write_intScalar("is_transport", is_transport ? 1 : 0);
    cmdh5w->write_string("lpn_file", lpn_file);
    cmdh5w->write_string("inflow_file", inflow_file);
    cmdh5w->write_string("sol_bName", sol_bName);
    cmdh5w->write_string("transport_name", transport_name);
    cmdh5w->write_string("hi_model", hi_model_name);
    cmdh5w->write_doubleScalar("hi_beta",
        restart_hi_model == nullptr ? 1.0 : restart_hi_model->get_beta());
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== Data from Files =====
  // Control points' xyz coordinates
  auto fNode = SYS_T::make_unique<FEANode>(part_file, rank);

  // Local sub-domain's IEN array
  auto locIEN = SYS_T::make_unique<ALocal_IEN>(part_file, rank);
  
  // Local sub-domain's element indices
  auto locElem = SYS_T::make_unique<ALocal_Elem>(part_file, rank);

  // Local sub-domain's nodal bc
  auto locnbc = SYS_T::make_unique<ALocal_NBC>(part_file, rank);

  // Local sub-domain's inflow bc
  auto locinfnbc = SYS_T::make_unique<ALocal_InflowBC>(part_file, rank);

  // Local sub-domain's elemental bc
  std::unique_ptr<ALocal_EBC> locebc = SYS_T::make_unique<ALocal_EBC_outflow>(part_file, rank);

  // Local sub_domain's weak bc
  auto locwbc = SYS_T::make_unique<ALocal_WeakBC>(part_file, rank);
  locwbc -> print_info();

  // Local sub-domain's nodal indices
  auto pNode = SYS_T::make_unique<APart_Node>(part_file, rank);

  SYS_T::commPrint("===> Data from HDF5 files are read from disk.\n");

  SYS_T::print_fatal_if( size!= ANL_T::get_cpu_size(part_file, rank),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for FEM analysis. \n", size);

  // ===== Inflow flow rate =====
  SYS_T::commPrint("===> Setup inflow flow rate. \n");

  auto inflow_rate = FlowRateFactory::createFlowRate(inflow_file);

  inflow_rate->print_info();

  // ===== LPN models =====
  auto gbc = GenBCFactory::createGenBC(lpn_file, initial_time, initial_step, 
      initial_index, 1000);

  gbc -> print_info();

  // Make sure the gbc number of faces matches that of ALocal_EBC
  SYS_T::print_fatal_if(gbc->get_num_ebc() != locebc->get_num_ebc(),
      "Error: GenBC number of faces does not match with that in ALocal_EBC.\n");

  // ===== Generate a sparse matrix for the enforcement of essential BCs
  auto pmat = SYS_T::make_unique<Matrix_PETSc>(pNode.get(), locnbc.get());

  pmat->gen_perm_bc(pNode.get(), locnbc.get());

  // ===== Generalized-alpha =====
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");

  auto tm_galpha = is_backward_Euler
    ? SYS_T::make_unique<TimeMethod_GenAlpha>(1.0, 1.0, 1.0)
    : SYS_T::make_unique<TimeMethod_GenAlpha>(genA_rho_inf, false);

  tm_galpha->print_info();

  // ===== Local Assembly routine =====
  std::unique_ptr<IPLocAssem> locAssem_ptr = nullptr;

  if( locwbc->get_wall_model_type() == 0 )
  {
    locAssem_ptr = SYS_T::make_unique<PLocAssem_VMS_NS_GenAlpha>(
      ANL_T::get_elemType(part_file, rank), nqp_vol, nqp_sur,
      tm_galpha.get(), fluid_density, fluid_mu, bs_beta, c_ct, c_tauc );    
  }
  else if( locwbc->get_wall_model_type() == 1 )
  {
    locAssem_ptr = SYS_T::make_unique<PLocAssem_VMS_NS_GenAlpha_WeakBC>(
      ANL_T::get_elemType(part_file, rank), nqp_vol, nqp_sur,
      tm_galpha.get(), fluid_density, fluid_mu, bs_beta, c_ct, c_tauc, C_bI );    
  }
  else SYS_T::print_fatal("Error: Unknown wall model type.\n");

  // ===== Initial condition =====
  std::unique_ptr<PDNSolution> base = SYS_T::make_unique<PDNSolution_NS>( 
      pNode.get(), fNode.get(), locinfnbc.get(), 1 );

  std::unique_ptr<PDNSolution> sol = nullptr;
  std::unique_ptr<PDNSolution> dot_sol = nullptr;
  NS_INIT::initialize_solution_state(pNode.get(), is_restart, 
      restart_index, restart_time, restart_step, restart_name, 
      sol, dot_sol, initial_index, initial_time, initial_step);

  std::unique_ptr<PDNSolution> transport_sol = nullptr;
  std::unique_ptr<PDNSolution> transport_dot_sol = nullptr;
  std::unique_ptr<PGAssem_Transport_NS_FEM> transport_gassem = nullptr;
  std::unique_ptr<PLinear_Solver_PETSc> transport_lsolver = nullptr;

  if( is_transport )
  {
    transport_sol = SYS_T::make_unique<PDNSolution>( pNode.get(), 1 );
    transport_dot_sol = SYS_T::make_unique<PDNSolution>( pNode.get(), 1 );
    VecSet(transport_sol->solution, 0.0);
    transport_sol->Assembly_GhostUpdate();
    VecSet(transport_dot_sol->solution, 0.0);
    transport_dot_sol->Assembly_GhostUpdate();

    if( is_restart )
    {
      const std::string transport_restart_name =
          transport_name + SYS_T::fixed_length_index(restart_index);
      const std::string dot_transport_restart_name =
          "dot_" + transport_name + SYS_T::fixed_length_index(restart_index);

      if( SYS_T::file_exist(transport_restart_name) )
      {
        transport_sol->ReadBinary(transport_restart_name);
        if( restart_hi_model != nullptr )
          HI_T::convert_HI_to_D(transport_sol.get(), restart_hi_model->get_beta());
      }

      if( SYS_T::file_exist(dot_transport_restart_name) )
      {
        transport_dot_sol->ReadBinary(dot_transport_restart_name);
      }
    }
  }

  // ===== Global assembly =====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  std::unique_ptr<IPGAssem> gloAssem = SYS_T::make_unique<PGAssem_NS_FEM>( 
      gbc.get(), std::move(locIEN), std::move(locElem), std::move(fNode), 
      std::move(pNode), std::move(locnbc), std::move(locebc), 
      std::move(locwbc), std::move(locAssem_ptr), nz_estimate );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem->Assem_nonzero_estimate( gbc.get() );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem->Fix_nonzero_err_str();
  gloAssem->Clear_KG();

  if( is_transport )
  {
    auto transport_locIEN = SYS_T::make_unique<ALocal_IEN>(part_file, rank);
    auto transport_locElem = SYS_T::make_unique<ALocal_Elem>(part_file, rank);
    auto transport_fNode = SYS_T::make_unique<FEANode>(part_file, rank);
    auto transport_pNode = SYS_T::make_unique<APart_Node>(part_file, rank);
    auto transport_locebc = SYS_T::make_unique<ALocal_EBC_outflow>(part_file, rank);
    auto transport_locAssem =
      SYS_T::make_unique<PLocAssem_Transport_VMS_NS_GenAlpha>(
          ANL_T::get_elemType(part_file, rank), nqp_vol, transport_ct,
          fluid_mu, tm_galpha.get(), std::move(restart_hi_model) );

    transport_gassem = SYS_T::make_unique<PGAssem_Transport_NS_FEM>(
        locinfnbc.get(), std::move(transport_locebc),
        std::move(transport_locIEN), std::move(transport_locElem),
        std::move(transport_fNode), std::move(transport_pNode),
        ANL_T::get_elemType(part_file, rank), nqp_sur,
        std::move(transport_locAssem), nz_estimate );

    transport_gassem->Fix_nonzero_err_str();
    transport_gassem->Clear_KG();

    transport_lsolver = SYS_T::make_unique<PLinear_Solver_PETSc>(
        1.0e-5, 1.0e-50, 1.0e50, 10000,
        "transport_", "transport_" );
  }

  // ===== Initialize the dot_sol vector by solving mass matrix =====
  NS_INIT::initialize_dot_solution(gloAssem.get(), sol.get(), dot_sol.get(), is_restart);

  // ===== Linear solver context =====
  auto lsolver = SYS_T::make_unique<PLinear_Solver_PETSc>();

  PC upc; lsolver->GetPC(&upc);
  const PetscInt pfield[1] = {0}, vfields[] = {1,2,3};
  PCFieldSplitSetBlockSize(upc,4);
  PCFieldSplitSetFields(upc,"u",3,vfields,vfields);
  PCFieldSplitSetFields(upc,"p",1,pfield,pfield);

  // ===== Nonlinear solver context =====
  auto nsolver = SYS_T::make_unique<PNonlinear_NS_Solver>(
      std::move(lsolver), std::move(pmat), std::move(tm_galpha), 
      std::move(inflow_rate), std::move(base), nl_rtol, nl_atol, 
      nl_dtol, nl_maxits, nl_refreq, nl_threshold );

  nsolver->print_info();

  // ===== Time step info =====
  auto timeinfo = SYS_T::make_unique<PDNTimeStep>(initial_index, initial_time, 
      initial_step);

  // ===== Temporal solver context =====
  auto tsolver = SYS_T::make_unique<PTime_NS_Solver>(
      std::move(nsolver), sol_bName, sol_record_freq, 
      ttan_renew_freq, final_time );

  tsolver->print_info();

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== Outlet data recording files =====
  tsolver->record_outlet_data(sol.get(), dot_sol.get(), timeinfo.get(), gbc.get(),
      gloAssem.get(), true, is_restart);

  // ===== Inlet data recording files =====
  tsolver->record_inlet_data(sol.get(), timeinfo.get(), locinfnbc.get(), 
      gloAssem.get(), true, is_restart);

  if( is_transport && transport_gassem->has_hi_model() )
  {
    tsolver->record_outlet_HI_data(sol.get(), transport_sol.get(), timeinfo.get(),
        gbc.get(), gloAssem.get(), transport_gassem.get(), fluid_density,
        true, is_restart);
  }

  // ===== FEM analysis =====
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
  if( is_transport )
  {
    tsolver->TM_NS_Transport_GenAlpha(is_restart, std::move(dot_sol), std::move(sol),
        std::move(transport_dot_sol), std::move(transport_sol), std::move(timeinfo),
        locinfnbc.get(), gbc.get(), gloAssem.get(), transport_gassem.get(),
        transport_lsolver.get(), transport_name, fluid_density );
  }
  else
  {
    tsolver->TM_NS_GenAlpha(is_restart, std::move(dot_sol), std::move(sol),
        std::move(timeinfo), locinfnbc.get(), gbc.get(), gloAssem.get() );
  }

  // ===== Print complete solver info =====
  tsolver -> print_lsolver_info();

  tsolver.reset(); locinfnbc.reset(); gbc.reset(); gloAssem.reset();

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
