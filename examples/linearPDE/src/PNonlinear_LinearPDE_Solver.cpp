#include "PNonlinear_LinearPDE_Solver.hpp"

PNonlinear_LinearPDE_Solver::PNonlinear_LinearPDE_Solver(
    std::unique_ptr<IPGAssem> in_gassem,
    std::unique_ptr<PLinear_Solver_PETSc> in_lsolver,
    std::unique_ptr<Matrix_PETSc> in_bc_mat,
    std::unique_ptr<TimeMethod_GenAlpha> in_tmga,
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol, const int &input_max_iteration,
    const int &input_renew_freq,
    const int &input_renew_threshold )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq),
  nrenew_threshold(input_renew_threshold),
  gassem(std::move(in_gassem)),
  lsolver(std::move(in_lsolver)),
  bc_mat(std::move(in_bc_mat)),
  tmga(std::move(in_tmga))
{}

void PNonlinear_LinearPDE_Solver::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("Nonlinear solver setted up:\n");
  SYS_T::commPrint("  relative tolerance: %e \n", nr_tol);
  SYS_T::commPrint("  absolute tolerance: %e \n", na_tol);
  SYS_T::commPrint("  divergence tolerance: %e \n", nd_tol);
  SYS_T::commPrint("  maximum iteration: %d \n", nmaxits);
  SYS_T::commPrint("  tangent matrix renew frequency: %d \n", nrenew_freq);
  SYS_T::commPrint("  tangent matrix renew threshold: %d \n", nrenew_threshold);
  SYS_T::commPrint("----------------------------------------------------------- \n");
}

void PNonlinear_LinearPDE_Solver::GenAlpha_Solve_Transport(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &pre_dot_sol,
    const PDNSolution * const &pre_sol,
    PDNSolution * const &dot_sol,
    PDNSolution * const &sol,
    bool &conv_flag, int &nl_counter ) const
{
  // Initialize the counter and error
  nl_counter = 0;
  double residual_norm = 0.0, initial_norm = 0.0, relative_error = 0.0;

  // Gen-alpha parameters
  const double gamma   = tmga->get_gamma();
  const double alpha_m = tmga->get_alpha_m();
  const double alpha_f = tmga->get_alpha_f();

  // Same-Y predictor
  sol     -> Copy(*pre_sol);
  dot_sol -> Copy(*pre_dot_sol);
  dot_sol -> ScaleValue( (gamma-1.0)/gamma );

  // Define the dol_sol at alpha_m: dot_sol_alpha
  PDNSolution dot_sol_alpha(*pre_dot_sol);
  dot_sol_alpha.ScaleValue( 1.0 - alpha_m );
  dot_sol_alpha.PlusAX(*dot_sol, alpha_m);

  // Define the sol at alpha_f: sol_alpha
  PDNSolution sol_alpha(*pre_sol);
  sol_alpha.ScaleValue( 1.0 - alpha_f );
  sol_alpha.PlusAX( *sol, alpha_f );

  // If new_tangent_flag == TRUE, update the tangent matrix;
  // otherwise, use the matrix from the previous time step
  if( new_tangent_flag )
  {
    gassem->Clear_KG();

    gassem->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha,
        curr_time, dt );

    SYS_T::commPrint("  --- M updated");

    // SetOperator will pass the tangent matrix to the linear solver and the
    // linear solver will generate the preconditioner based on the new matrix.
    lsolver->SetOperator( gassem->K );
  }
  else
  {
    gassem->Clear_G();

    gassem->Assem_residual( &dot_sol_alpha, &sol_alpha,
        curr_time, dt );
  }

  VecNorm( gassem->G, NORM_2, &initial_norm );
  SYS_T::commPrint("  Init res 2-norm: %e \n", initial_norm);

  PDNSolution * dot_step = new PDNSolution( pre_sol );

  // Now do the Newton-Raphson iteration (multi-corrector stage)
  do
  {
    // solve the equation : K dot_step = G
    lsolver->Solve( gassem->G, dot_step );

    bc_mat -> MatMultSol( dot_step );

    nl_counter += 1;

    dot_sol->PlusAX( dot_step, -1.0 );
    sol->PlusAX( dot_step, (-1.0) * gamma * dt );

    dot_sol_alpha.PlusAX( dot_step, (-1.0) * alpha_m );
    sol_alpha.PlusAX( dot_step, (-1.0) * alpha_f * gamma * dt );

    // Assembly residual (& tangent if condition satisfied)
    if( nl_counter % nrenew_freq == 0 || nl_counter >= nrenew_threshold )
    {
      gassem->Clear_KG();

      gassem->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha,
          curr_time, dt );

      SYS_T::commPrint("  --- M updated");
      lsolver->SetOperator(gassem->K);
    }
    else
    {
      gassem->Clear_G();

      gassem->Assem_residual( &dot_sol_alpha, &sol_alpha,
          curr_time, dt );
    }

    VecNorm(gassem->G, NORM_2, &residual_norm);

    SYS_T::print_fatal_if( residual_norm != residual_norm, "Error: nonlinear solver residual norm is NaN. Job killed.\n" );

    SYS_T::commPrint("  --- nl_res: %e \n", residual_norm);

    relative_error = residual_norm / initial_norm;

    SYS_T::print_fatal_if( relative_error >= nd_tol, "Error: nonlinear solver is diverging with error %e. Job killed.\n", relative_error);

  }while(nl_counter<nmaxits && relative_error > nr_tol && residual_norm > na_tol);

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol) conv_flag = true;
  else conv_flag = false;

  delete dot_step;
}

void PNonlinear_LinearPDE_Solver::GenAlpha_Solve_Elastodynamics(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &pre_dot_disp,
    const PDNSolution * const &pre_dot_velo,
    const PDNSolution * const &pre_disp,
    const PDNSolution * const &pre_velo,
    PDNSolution * const &dot_disp,
    PDNSolution * const &dot_velo,
    PDNSolution * const &disp,
    PDNSolution * const &velo,
    bool &conv_flag, int &nl_counter ) const
{
  // Initialize the counter and error
  nl_counter = 0;
  double residual_norm = 0.0, initial_norm = 0.0, relative_error = 0.0;

  // Gen-alpha parameters
  const double gamma   = tmga->get_gamma();
  const double alpha_m = tmga->get_alpha_m();
  const double alpha_f = tmga->get_alpha_f();

  // Same-Velocity predictor
  velo -> Copy( pre_velo );
  dot_velo -> Copy( pre_dot_velo ); dot_velo -> ScaleValue( (gamma-1.0)/gamma );

  // get dot_disp by the kinematic equations
  dot_disp -> Copy( pre_dot_disp );
  dot_disp -> ScaleValue( alpha_m-1.0 );
  dot_disp -> PlusAX( pre_velo, 1.0 );
  dot_disp -> ScaleValue( 1.0/alpha_m );

  disp -> Copy( pre_disp );
  disp -> PlusAX( pre_dot_disp, dt - gamma * dt / alpha_m );
  disp -> PlusAX( pre_velo, gamma * dt / alpha_m );

  // Define intermediate solutions
  PDNSolution * dot_velo_alpha = new PDNSolution( pre_dot_velo );
  dot_velo_alpha -> ScaleValue( 1.0 - alpha_m );
  dot_velo_alpha -> PlusAX( dot_velo, alpha_m );

  PDNSolution * disp_alpha = new PDNSolution( pre_disp );
  disp_alpha -> ScaleValue( 1.0 - alpha_f );
  disp_alpha -> PlusAX( disp, alpha_f );

  // If new_tangent_flag == TRUE, update the tangent matrix;
  // otherwise, use the matrix from the previous time step
  if( new_tangent_flag )
  {
    gassem->Clear_KG();

    gassem->Assem_tangent_residual( dot_velo_alpha, disp_alpha,
        curr_time, dt );

    SYS_T::commPrint("  --- M updated");

    // SetOperator will pass the tangent matrix to the linear solver and the
    // linear solver will generate the preconditioner based on the new matrix.
    lsolver->SetOperator( gassem->K );
  }
  else
  {
    gassem->Clear_G();

    gassem->Assem_residual( dot_velo_alpha, disp_alpha,
        curr_time, dt );
  }

  VecNorm( gassem->G, NORM_2, &initial_norm );
  SYS_T::commPrint("  Init res 2-norm: %e \n", initial_norm);

  PDNSolution * dot_step = new PDNSolution( pre_velo );

  // Now do the Newton-Raphson iteration (multi-corrector stage)
  do
  {
    // solve the equation : K dot_step = G
    lsolver->Solve( gassem->G, dot_step );

    bc_mat -> MatMultSol( dot_step );

    nl_counter += 1;

    dot_velo->PlusAX( dot_step, -1.0 );
    dot_disp->PlusAX( dot_step, - alpha_f * gamma * dt / alpha_m );
    velo->PlusAX( dot_step, - gamma * dt );
    disp->PlusAX( dot_step, - alpha_f * gamma * gamma * dt * dt / alpha_m );

    dot_velo_alpha->PlusAX( dot_step, - alpha_m );
    disp_alpha->PlusAX( dot_step, - alpha_f * alpha_f * gamma * gamma * dt * dt / alpha_m);

    // Assembly residual (& tangent if condition satisfied)
    if( nl_counter % nrenew_freq == 0 || nl_counter >= nrenew_threshold )
    {
      gassem->Clear_KG();

      gassem->Assem_tangent_residual( dot_velo_alpha, disp_alpha,
          curr_time, dt );

      SYS_T::commPrint("  --- M updated");
      lsolver->SetOperator(gassem->K);
    }
    else
    {
      gassem->Clear_G();

      gassem->Assem_residual( dot_velo_alpha, disp_alpha,
          curr_time, dt );
    }

    VecNorm(gassem->G, NORM_2, &residual_norm);

    SYS_T::print_fatal_if( residual_norm != residual_norm, "Error: nonlinear solver residual norm is NaN. Job killed.\n" );

    SYS_T::commPrint("  --- nl_res: %e \n", residual_norm);

    relative_error = residual_norm / initial_norm;

    SYS_T::print_fatal_if( relative_error >= nd_tol, "Error: nonlinear solver is diverging with error %e. Job killed.\n", relative_error);

  }while(nl_counter<nmaxits && relative_error > nr_tol && residual_norm > na_tol);

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol) conv_flag = true;
  else conv_flag = false;

  delete dot_step; delete dot_velo_alpha; delete disp_alpha;
}

// EOF
