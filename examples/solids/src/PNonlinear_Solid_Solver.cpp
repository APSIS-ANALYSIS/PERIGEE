#include "PNonlinear_Solid_Solver.hpp"

PNonlinear_Solid_Solver::PNonlinear_Solid_Solver(
    std::unique_ptr<PGAssem_Solid_FEM> in_gassem,
    std::unique_ptr<PLinear_Solver_PETSc> in_lsolver,
    std::unique_ptr<Matrix_PETSc> in_bc_mat,
    std::unique_ptr<TimeMethod_GenAlpha> in_tmga,
    std::unique_ptr<APart_Node> in_pnode,
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol, const int &input_max_iteration,
    const int &input_renew_freq, const int &input_renew_threshold )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq),
  nrenew_threshold(input_renew_threshold),
  gassem(std::move(in_gassem)),
  lsolver(std::move(in_lsolver)),
  bc_mat(std::move(in_bc_mat)),
  tmga(std::move(in_tmga)),
  pnode(std::move(in_pnode))
{}

void PNonlinear_Solid_Solver::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint("relative tolerance: %e \n", nr_tol);
  SYS_T::commPrint("absolute tolerance: %e \n", na_tol);
  SYS_T::commPrint("divergence tolerance: %e \n", nd_tol);
  SYS_T::commPrint("maximum iteration: %d \n", nmaxits);
  SYS_T::commPrint("tangent matrix renew frequency: %d \n", nrenew_freq);
  SYS_T::commPrint("tangent matrix renew threshold: %d \n", nrenew_threshold);
  SYS_T::print_sep_line();
}

void PNonlinear_Solid_Solver::update_solid_kinematics( const double &val,
    const Vec &input,
    PDNSolution * const &output ) const
{
  const int nlocal = pnode->get_nlocalnode();

  Vec local_output;
  VecGhostGetLocalForm(output->solution, &local_output);

  double * array_input, * array_output;
  VecGetArray(input,        &array_input);
  VecGetArray(local_output, &array_output);

  for(int ii=0; ii<nlocal; ++ii)
  {
    const int ii3 = ii * 3;
    array_output[ii3  ] += val * array_input[ii3  ];
    array_output[ii3+1] += val * array_input[ii3+1];
    array_output[ii3+2] += val * array_input[ii3+2];
  }

  VecRestoreArray(input,        &array_input);
  VecRestoreArray(local_output, &array_output);
  VecGhostRestoreLocalForm(output->solution, &local_output);

  output->GhostUpdate();
}

void PNonlinear_Solid_Solver::GenAlpha_Seg_solve_Solid(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const IS &is_v,
    const IS &is_p,
    const PDNSolution * const &pre_dot_disp,
    const PDNSolution * const &pre_dot_velo,
    const PDNSolution * const &pre_dot_pres,
    const PDNSolution * const &pre_disp,
    const PDNSolution * const &pre_velo,
    const PDNSolution * const &pre_pres,
    PDNSolution * const &dot_disp,
    PDNSolution * const &dot_velo,
    PDNSolution * const &dot_pres,
    PDNSolution * const &disp,
    PDNSolution * const &velo,
    PDNSolution * const &pres,
    bool &conv_flag, int &nl_counter ) const
{
  // Initialization
  nl_counter = 0;
  double residual_norm = 0.0, initial_norm = 0.0, relative_error = 0.0;

  const double gamma   = tmga->get_gamma();
  const double alpha_m = tmga->get_alpha_m();
  const double alpha_f = tmga->get_alpha_f();

  const double val_1 = alpha_f * gamma * dt / alpha_m;

  // Same-Y predictor
  dot_disp -> Copy( pre_dot_disp ); dot_disp -> ScaleValue( (gamma-1.0)/gamma );
  dot_velo -> Copy( pre_dot_velo ); dot_velo -> ScaleValue( (gamma-1.0)/gamma );
  dot_pres -> Copy( pre_dot_pres ); dot_pres -> ScaleValue( (gamma-1.0)/gamma );

  disp -> Copy( pre_disp );
  velo -> Copy( pre_velo );
  pres -> Copy( pre_pres );

  gassem->Apply_Dirichlet_BC( curr_time + dt, dot_disp, dot_velo, disp, velo );

  // Define intermediate solutions
  auto dot_disp_alpha = SYS_T::make_unique<PDNSolution>(pre_dot_disp);
  dot_disp_alpha -> ScaleValue( 1.0 - alpha_m );
  dot_disp_alpha -> PlusAX( dot_disp, alpha_m );

  auto dot_velo_alpha = SYS_T::make_unique<PDNSolution>(pre_dot_velo);
  dot_velo_alpha -> ScaleValue( 1.0 - alpha_m );
  dot_velo_alpha -> PlusAX( dot_velo, alpha_m );

  auto dot_pres_alpha = SYS_T::make_unique<PDNSolution>(pre_dot_pres);
  dot_pres_alpha -> ScaleValue( 1.0 - alpha_m );
  dot_pres_alpha -> PlusAX( dot_pres, alpha_m );

  auto disp_alpha = SYS_T::make_unique<PDNSolution>(pre_disp);
  disp_alpha -> ScaleValue( 1.0 - alpha_f );
  disp_alpha -> PlusAX( disp, alpha_f );

  auto velo_alpha = SYS_T::make_unique<PDNSolution>(pre_velo);
  velo_alpha -> ScaleValue( 1.0 - alpha_f );
  velo_alpha -> PlusAX( velo, alpha_f );

  auto pres_alpha = SYS_T::make_unique<PDNSolution>(pre_pres);
  pres_alpha -> ScaleValue( 1.0 - alpha_f );
  pres_alpha -> PlusAX( pres, alpha_f );

  // Get Delta_dot_disp by assuming Delta_v is zero
  auto Delta_dot_disp = SYS_T::make_unique<PDNSolution>(pre_disp);
  Delta_dot_disp -> ScaleValue( 0.0 );

  update_solid_kinematics( -1.0 / alpha_m, dot_disp_alpha->solution, Delta_dot_disp.get() );
  update_solid_kinematics(  1.0 / alpha_m,     velo_alpha->solution, Delta_dot_disp.get() );

  // Update displacement solutions
  dot_disp       -> PlusAX( Delta_dot_disp.get(), 1.0 );
  disp           -> PlusAX( Delta_dot_disp.get(), gamma * dt );
  dot_disp_alpha -> PlusAX( Delta_dot_disp.get(), alpha_m );
  disp_alpha     -> PlusAX( Delta_dot_disp.get(), alpha_f * gamma * dt );

  // Assemble residual (& tangent if needed)
  if( new_tangent_flag )
  {
    gassem -> Clear_KG();
    gassem -> Assem_Tangent_Residual( curr_time, dt,
        dot_disp_alpha.get(), dot_velo_alpha.get(), dot_pres_alpha.get(),
        disp_alpha.get(), velo_alpha.get(), pres_alpha.get() );

    SYS_T::commPrint("  --- M updated");
    lsolver->SetOperator(gassem->K);
  }
  else
  {
    gassem -> Clear_G();
    gassem -> Assem_Residual( curr_time, dt,
        dot_disp_alpha.get(), dot_velo_alpha.get(), dot_pres_alpha.get(),
        disp_alpha.get(), velo_alpha.get(), pres_alpha.get() );
  }

  VecNorm(gassem->G, NORM_2, &initial_norm);
  SYS_T::commPrint("  Init res 2-norm: %e \n", initial_norm);

  Vec sol_vp, sol_v, sol_p;
  VecDuplicate(gassem->G, &sol_vp);

  // Newton-Raphson iteration
  do
  {
    lsolver->Solve( gassem->G, sol_vp );
    bc_mat -> MatMultSol( sol_vp );

    nl_counter += 1;

    VecGetSubVector(sol_vp, is_v, &sol_v);
    VecGetSubVector(sol_vp, is_p, &sol_p);

    dot_velo       -> PlusAX( sol_v, -1.0 );
    dot_velo_alpha -> PlusAX( sol_v, -1.0 * alpha_m );
    velo           -> PlusAX( sol_v, -1.0 * gamma * dt );
    velo_alpha     -> PlusAX( sol_v, -1.0 * gamma * alpha_f * dt );

    dot_pres       -> PlusAX( sol_p, -1.0 );
    dot_pres_alpha -> PlusAX( sol_p, -1.0 * alpha_m );
    pres           -> PlusAX( sol_p, -1.0 * gamma * dt );
    pres_alpha     -> PlusAX( sol_p, -1.0 * gamma * alpha_f * dt );

    update_solid_kinematics( -1.0 * val_1, sol_v, dot_disp );
    update_solid_kinematics( -1.0 * val_1 * gamma * dt, sol_v, disp );
    update_solid_kinematics( -1.0 * val_1 * alpha_m, sol_v, dot_disp_alpha.get() );
    update_solid_kinematics( -1.0 * val_1 * alpha_f * gamma * dt, sol_v, disp_alpha.get() );

    VecRestoreSubVector(sol_vp, is_v, &sol_v);
    VecRestoreSubVector(sol_vp, is_p, &sol_p);

    if( nl_counter % nrenew_freq == 0 || nl_counter >= nrenew_threshold )
    {
      gassem -> Clear_KG();
      gassem -> Assem_Tangent_Residual( curr_time, dt,
          dot_disp_alpha.get(), dot_velo_alpha.get(), dot_pres_alpha.get(),
          disp_alpha.get(), velo_alpha.get(), pres_alpha.get() );

      SYS_T::commPrint("  --- M updated");
      lsolver->SetOperator(gassem->K);
    }
    else
    {
      gassem -> Clear_G();
      gassem -> Assem_Residual( curr_time, dt,
          dot_disp_alpha.get(), dot_velo_alpha.get(), dot_pres_alpha.get(),
          disp_alpha.get(), velo_alpha.get(), pres_alpha.get() );
    }

    VecNorm(gassem->G, NORM_2, &residual_norm);
    SYS_T::commPrint("  --- nl_res: %e \n", residual_norm);

    relative_error = residual_norm / initial_norm;

    if( relative_error >= nd_tol )
    {
      SYS_T::commPrint("Warning: nonlinear solver is diverging with error %e. \n", relative_error);
      break;
    }

  }while( nl_counter < nmaxits && relative_error > nr_tol && residual_norm > na_tol );

  gassem->Apply_Dirichlet_BC( curr_time + dt, dot_disp, dot_velo, disp, velo );

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if( relative_error <= nr_tol || residual_norm <= na_tol ) conv_flag = true;
  else conv_flag = false;

  VecDestroy(&sol_vp);
}

// EOF
