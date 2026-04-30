// ============================================================================
// driver.cpp
//
// Stablized finite element code for 3D hyperelastic solid problems.
//
// Date: Jan. 29 2026
// ==========================================================================
#include "HDF5_Writer.hpp"
#include "ANL_Tools.hpp"
#include "MaterialModel_ich_NeoHookean.hpp"
#include "MaterialModel_vol_Incompressible.hpp"
#include "MaterialModel_Mixed_Elasticity.hpp"
#include "PLocAssem_2x2Block_VMS_Incompressible.hpp"
#include "InitHelpers.hpp"
#include "ALocal_NBC.hpp"
#include "PGAssem_Solid_FEM.hpp"
#include "PNonlinear_Solid_Solver.hpp"
#include "PTime_Solid_Solver.hpp"

int main(int argc, char *argv[])
{
  // number of quadrature points
  int nqp_vol = 5, nqp_sur = 4;

  // generalized-alpha rho_inf
  double genA_rho_inf = 0.5;
  bool is_backward_Euler = false;

  // Estimate of the nonzero per row for the sparse matrix
  int nz_estimate = 300;

  // part file location
  std::string part_file("./apart/part");

  // nonlinear solver parameters
  double nl_rtol = 1.0e-10; // convergence criterion relative tolerance
  double nl_atol = 1.0e-10; // convergence criterion absolute tolerance
  double nl_dtol = 10.0;   // divergence criterion
  int nl_maxits = 20;      // maximum number if nonlinear iterations
  int nl_refreq = 4;       // frequency of tangent matrix renewal
  int nl_threshold = 4;    // threshold of tangent matrix renewal

  // time stepping parameters
  double initial_time = 0.0; // time of the initial condition
  double initial_step = 0.01; // time step size
  int initial_index = 0;     // indiex of the initial condition
  double final_time = 1.0;   // final time
  std::string sol_bName("SOL_"); // base name of the solution file
  int ttan_renew_freq = 1;   // frequency of tangent matrix renewal
  int sol_record_freq = 1;   // frequency of recording the solution

  // solid material parameters
  const double solid_density = 1.0e3;
  const double solid_mu = 6.666666666e4;

  // displacement-driven BC parameters (edit here)

  // Restart options
  bool is_restart = false;
  int restart_index = 0;     // restart solution time index
  double restart_time = 0.0; // restart time
  double restart_step = 1.0e-3; // restart simulation time step size
  const std::string restart_u_name = "SOL_disp_";
  const std::string restart_v_name = "SOL_velo_";
  const std::string restart_p_name = "SOL_pres_";

  // Yaml options
  bool is_loadYaml = true;
  std::string yaml_file("./driver.yml");

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

  // ===== Yaml Arguments =====
  SYS_T::GetOptionBool("-is_loadYaml", is_loadYaml);
  SYS_T::GetOptionString("-yaml_file", yaml_file);

  if(is_loadYaml) SYS_T::InsertFileYAML( yaml_file,  false );

  // ===== Read Command Line Arguments =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionReal("-rho_inf", genA_rho_inf);
  SYS_T::GetOptionBool("-is_backward_Euler", is_backward_Euler);
  SYS_T::GetOptionInt("-nqp_vol", nqp_vol);
  SYS_T::GetOptionInt("-nqp_sur", nqp_sur);
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

  // ===== Print Command Line Arguments =====
  if( is_backward_Euler )
    SYS_T::commPrint(   "-is_backward_Euler: true \n");
  else
    SYS_T::cmdPrint(    "-rho_inf:", genA_rho_inf);

  SYS_T::cmdPrint("-nz_estimate:", nz_estimate);
  SYS_T::cmdPrint("-nqp_vol:", nqp_vol);
  SYS_T::cmdPrint("-nqp_sur:", nqp_sur);
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
  }
  else SYS_T::commPrint("-is_restart: false \n");

  MPI_Barrier(PETSC_COMM_WORLD);

  SYS_T::print_fatal_if( size != ANL_T::get_cpu_size(part_file, rank),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for FEM analysis.\n", size);

  const FEType elemType = ANL_T::get_elemType(part_file, rank);

  // ===== Record important solver options =====
  if(rank == 0)
  {
    auto cmdh5w = SYS_T::make_unique<HDF5_Writer>("solver_cmd.h5");

    const double nstep_real = (final_time - initial_time) / initial_step;
    int nstep = static_cast<int>(nstep_real);
    if( initial_time + nstep * initial_step < final_time - 1.0e-12 ) nstep += 1;
    const int final_index = initial_index + nstep;

    cmdh5w->write_doubleScalar("init_time", initial_time);
    cmdh5w->write_intScalar("init_index", initial_index);
    cmdh5w->write_doubleScalar("init_step", initial_step);
    cmdh5w->write_doubleScalar("final_time", final_time);
    cmdh5w->write_intScalar("final_index", final_index);
    cmdh5w->write_intScalar("sol_record_freq", sol_record_freq);
    cmdh5w->write_intScalar("nqp_vol", nqp_vol);
    cmdh5w->write_intScalar("nqp_sur", nqp_sur);
  }

  auto locIEN = SYS_T::make_unique<ALocal_IEN>(part_file, rank);
  auto locElem = SYS_T::make_unique<ALocal_Elem>(part_file, rank);
  auto fNode = SYS_T::make_unique<FEANode>(part_file, rank);
  auto locnbc = SYS_T::make_unique<ALocal_NBC>(part_file, rank);
  auto locnbc_disp =
    SYS_T::make_unique<ALocal_NBC>(part_file, rank, "/nbc_disp_driven");
  auto locebc = SYS_T::make_unique<ALocal_EBC>(part_file, rank);
  auto pNode = SYS_T::make_unique<APart_Node>(part_file, rank);

  // ===== Generate a sparse matrix for the enforcement of essential BCs
  auto pmat = SYS_T::make_unique<Matrix_PETSc>(pNode.get(), locnbc.get());
  pmat->gen_perm_bc(pNode.get(), locnbc.get());

  const int dof_mat = locnbc->get_dof_LID();
  const int nlocal = pNode->get_nlocalnode();

  std::vector<PetscInt> idx_v(3 * nlocal);
  std::vector<PetscInt> idx_p(nlocal);
  for(int ii=0; ii<nlocal; ++ii)
  {
    const PetscInt gid = pNode->get_node_loc(ii);
    idx_p[ii] = dof_mat * gid;
    idx_v[3*ii  ] = dof_mat * gid + 1;
    idx_v[3*ii+1] = dof_mat * gid + 2;
    idx_v[3*ii+2] = dof_mat * gid + 3;
  }

  IS is_velo, is_pres;
  ISCreateGeneral(PETSC_COMM_WORLD, static_cast<PetscInt>(idx_v.size()), idx_v.data(), PETSC_COPY_VALUES, &is_velo);
  ISCreateGeneral(PETSC_COMM_WORLD, static_cast<PetscInt>(idx_p.size()), idx_p.data(), PETSC_COPY_VALUES, &is_pres);

  // ===== Generalized-alpha =====
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");

  auto tm_galpha = is_backward_Euler
    ? SYS_T::make_unique<TimeMethod_GenAlpha>(1.0, 1.0, 1.0)
    : SYS_T::make_unique<TimeMethod_GenAlpha>(genA_rho_inf, false);

  tm_galpha->print_info();

  // ===== Local Assembly Routine =====
  std::unique_ptr<IMaterialModel_ich> imodel =
    SYS_T::make_unique<MaterialModel_ich_NeoHookean>(solid_mu);

  std::unique_ptr<IMaterialModel_vol> vmodel =
    SYS_T::make_unique<MaterialModel_vol_Incompressible>(solid_density);

  std::unique_ptr<MaterialModel_Mixed_Elasticity> matmodel =
    SYS_T::make_unique<MaterialModel_Mixed_Elasticity>(std::move(vmodel), std::move(imodel));

  std::unique_ptr<IPLocAssem_2x2Block> locAssem_ptr =
    SYS_T::make_unique<PLocAssem_2x2Block_VMS_Incompressible>(
        elemType, nqp_vol, nqp_sur,
        tm_galpha.get(), std::move(matmodel));

  // ===== Initial condition =====
  auto disp = PDNSolution::Gen_zero_ptr( pNode.get(), 3 );
  auto velo = PDNSolution::Gen_zero_ptr( pNode.get(), 3 );
  auto pres = PDNSolution::Gen_zero_ptr( pNode.get(), 1 );

  auto dot_disp = PDNSolution::Gen_zero_ptr( pNode.get(), 3 );
  auto dot_velo = PDNSolution::Gen_zero_ptr( pNode.get(), 3 );
  auto dot_pres = PDNSolution::Gen_zero_ptr( pNode.get(), 1 );

  // ===== Global Assembly Routine =====
  std::unique_ptr<PGAssem_Solid_FEM> gloAssem_ptr =
    SYS_T::make_unique<PGAssem_Solid_FEM>(
        std::move(locIEN), std::move(locElem), std::move(fNode),
        std::move(pNode), std::move(locnbc), std::move(locebc),
        std::move(locAssem_ptr), nz_estimate);

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();

  SOLID_INIT::initialize_solution_state( is_restart,
      restart_index, restart_time, restart_step,
      restart_u_name, restart_v_name, restart_p_name,
      disp, velo, pres, dot_disp, dot_velo, dot_pres,
      initial_index, initial_time, initial_step );

  // ===== Initialize the dot_sol vectors by solving mass matrix =====
  SOLID_INIT::initialize_dot_solution( gloAssem_ptr.get(), is_velo, is_pres,
      dot_disp.get(), dot_velo.get(), dot_pres.get(),
      disp.get(), velo.get(), pres.get(), is_restart );

  // ===== Linear and nonlinear solver context =====
  auto lsolver = SYS_T::make_unique<PLinear_Solver_PETSc>();

  auto nsolver = SYS_T::make_unique<PNonlinear_Solver>(
      std::move(gloAssem_ptr), std::move(lsolver), std::move(pmat),
      std::move(tm_galpha),
      nl_rtol, nl_atol, nl_dtol, nl_maxits, nl_refreq, nl_threshold );

  nsolver->print_info();

  // ===== Time step info =====
  auto timeinfo = SYS_T::make_unique<PDNTimeStep>(initial_index, initial_time, initial_step);

  // ===== Temporal solver context =====
  auto tsolver = SYS_T::make_unique<PTime_Solver>(
      std::move(nsolver), sol_bName, sol_record_freq, ttan_renew_freq, final_time );

  tsolver->print_info();

  SYS_T::commPrint("===> Start Finite Element Analysis:\n");

  tsolver->TM_Solid_GenAlpha( is_restart, is_velo, is_pres,
      locnbc_disp.get(),
      std::move(dot_disp), std::move(dot_velo), std::move(dot_pres),
      std::move(disp), std::move(velo), std::move(pres),
      std::move(timeinfo) );

  ISDestroy(&is_velo);
  ISDestroy(&is_pres);

  // Ensure PETSc objects are destroyed before PetscFinalize
  tsolver.reset();

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
