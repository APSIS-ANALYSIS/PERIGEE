#include "PNonlinear_NS_Solver.hpp"

PNonlinear_NS_Solver::PNonlinear_NS_Solver(
    std::unique_ptr<PLinear_Solver_PETSc> in_lsolver,
    std::unique_ptr<Matrix_PETSc> in_bc_mat,
    std::unique_ptr<TimeMethod_GenAlpha> in_tmga,
    std::unique_ptr<IFlowRate> in_flrate,
    std::unique_ptr<PDNSolution> in_sol_base,
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol,
    const int &input_max_iteration, 
    const int &input_renew_freq,
    const int &input_renew_threshold )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq),
  nrenew_threshold(input_renew_threshold),
  lsolver(std::move(in_lsolver)),
  bc_mat(std::move(in_bc_mat)),
  tmga(std::move(in_tmga)),
  flrate(std::move(in_flrate)),
  sol_base(std::move(in_sol_base))
{}

void PNonlinear_NS_Solver::print_info() const
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


void PNonlinear_NS_Solver::GenAlpha_Solve_NS(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &pre_dot_sol,
    const PDNSolution * const &pre_sol,
    PDNSolution * const &dot_sol,
    PDNSolution * const &sol,
    const ALocal_InflowBC * const &infnbc_part,
    const IGenBC * const &gbc,
    IPGAssem * const &gassem_ptr,
    bool &conv_flag, int &nl_counter ) const
{
#ifdef PETSC_USE_LOG
  PetscLogEvent mat_assem_0_event, mat_assem_1_event;
  PetscLogEvent vec_assem_0_event, vec_assem_1_event;
  PetscLogEvent lin_solve_event;
  PetscClassId classid_assembly;
  PetscClassIdRegister("mat_vec_assembly", &classid_assembly);
  PetscLogEventRegister("assembly mat 0", classid_assembly, &mat_assem_0_event);
  PetscLogEventRegister("assembly mat 1", classid_assembly, &mat_assem_1_event);
  PetscLogEventRegister("assembly vec 0", classid_assembly, &vec_assem_0_event);
  PetscLogEventRegister("assembly vec 1", classid_assembly, &vec_assem_1_event);
  PetscLogEventRegister("lin_solve", classid_assembly, &lin_solve_event);
#endif

  // Initialize the counter and error
  nl_counter = 0;
  double residual_norm = 0.0, initial_norm = 0.0, relative_error = 0.0;

  // Gen-alpha parameters
  const double gamma   = tmga->get_gamma();
  const double alpha_m = tmga->get_alpha_m();
  const double alpha_f = tmga->get_alpha_f();

  // Same-Y predictor
  sol->Copy(*pre_sol);
  dot_sol->Copy(*pre_dot_sol);
  dot_sol->ScaleValue( (gamma-1.0)/gamma );

  // Define the dol_sol at alpha_m: dot_sol_alpha
  PDNSolution dot_sol_alpha(*pre_dot_sol);
  dot_sol_alpha.ScaleValue( 1.0 - alpha_m );
  dot_sol_alpha.PlusAX(*dot_sol, alpha_m);

  // Define the sol at alpha_f: sol_alpha
  PDNSolution sol_alpha(*pre_sol);
  sol_alpha.ScaleValue( 1.0 - alpha_f );
  sol_alpha.PlusAX( *sol, alpha_f );

  // ------------------------------------------------- 
  // Update the inflow boundary values
  rescale_inflow_value(curr_time+dt, infnbc_part, sol);
  rescale_inflow_value(curr_time+alpha_f*dt, infnbc_part, &sol_alpha);
  // ------------------------------------------------- 

  // If new_tangent_flag == TRUE, update the tangent matrix;
  // otherwise, use the matrix from the previous time step
  if( new_tangent_flag )
  {
    gassem_ptr->Clear_KG();

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(mat_assem_0_event, 0,0,0,0);
#endif

    gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol, 
        curr_time, dt, gbc );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(mat_assem_0_event,0,0,0,0);
#endif

    SYS_T::commPrint("  --- M updated");
    
    // SetOperator will pass the tangent matrix to the linear solver and the
    // linear solver will generate the preconditioner based on the new matrix.
    lsolver->SetOperator( gassem_ptr->K );
  }
  else
  {
    gassem_ptr->Clear_G();

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(vec_assem_0_event, 0,0,0,0);
#endif

    gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
        curr_time, dt, gbc );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(vec_assem_0_event,0,0,0,0);
#endif
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  SYS_T::commPrint("  Init res 2-norm: %e \n", initial_norm);

  auto dot_step = SYS_T::make_unique<PDNSolution>( pre_sol );

  // Now do consistent Newton-Raphson iteration
  do
  {
#ifdef PETSC_USE_LOG
    PetscLogEventBegin(lin_solve_event, 0,0,0,0);
#endif
    
    // solve the equation K dot_step = G
    lsolver->Solve( gassem_ptr->G, dot_step.get() );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(lin_solve_event,0,0,0,0);
#endif

    bc_mat->MatMultSol( dot_step.get() );

    nl_counter += 1;

    dot_sol->PlusAX( dot_step.get(), -1.0 );
    sol->PlusAX( dot_step.get(), (-1.0) * gamma * dt );

    dot_sol_alpha.PlusAX( dot_step.get(), (-1.0) * alpha_m );
    sol_alpha.PlusAX( dot_step.get(), (-1.0) * alpha_f * gamma * dt );

    // Assembly residual (& tangent if condition satisfied) 
    if( nl_counter % nrenew_freq == 0 || nl_counter >= nrenew_threshold )
    {
      gassem_ptr->Clear_KG();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(mat_assem_1_event, 0,0,0,0);
#endif

      gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
          curr_time, dt, gbc );

#ifdef PETSC_USE_LOG
      PetscLogEventEnd(mat_assem_1_event,0,0,0,0);
#endif

      SYS_T::commPrint("  --- M updated");
      lsolver->SetOperator(gassem_ptr->K);
    }
    else
    {
      gassem_ptr->Clear_G();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(vec_assem_1_event, 0,0,0,0);
#endif

      gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
          curr_time, dt, gbc );

#ifdef PETSC_USE_LOG
      PetscLogEventEnd(vec_assem_1_event,0,0,0,0);
#endif
    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    
    SYS_T::print_fatal_if( residual_norm != residual_norm, "Error: nonlinear solver residual norm is NaN. Job killed.\n" );
    
    SYS_T::commPrint("  --- nl_res: %e \n", residual_norm);

    relative_error = residual_norm / initial_norm;

    SYS_T::print_fatal_if( relative_error >= nd_tol, "Error: nonlinear solver is diverging with error %e. Job killed.\n", relative_error);

  }while(nl_counter<nmaxits && relative_error > nr_tol && residual_norm > na_tol);

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol) conv_flag = true;
  else conv_flag = false;
}

void PNonlinear_NS_Solver::rescale_inflow_value( const double &stime,
    const ALocal_InflowBC * const &infbc,
    PDNSolution * const &sol ) const
{
  const int num_nbc = infbc -> get_num_nbc();

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    const int numnode = infbc -> get_Num_LD( nbc_id );

    const double factor  = flrate -> get_flow_rate( nbc_id, stime );
    const double std_dev = flrate -> get_flow_TI_std_dev( nbc_id );

    for(int ii=0; ii<numnode; ++ii)
    {
      const int node_index = infbc -> get_LDN( nbc_id, ii );
      
      const int base_idx[3] = { node_index*4+1, node_index*4+2, node_index*4+3 };

      double base_vals[3];

      VecGetValues(sol_base->solution, 3, base_idx, base_vals);

      const double perturb_x = MATH_T::gen_double_rand_normal(0, std_dev);
      const double perturb_y = MATH_T::gen_double_rand_normal(0, std_dev);
      const double perturb_z = MATH_T::gen_double_rand_normal(0, std_dev);

      const double vals[3] = { base_vals[0] * factor * (1.0 + perturb_x), 
        base_vals[1] * factor * (1.0 + perturb_y),
        base_vals[2] * factor * (1.0 + perturb_z) };

      VecSetValues(sol->solution, 3, base_idx, vals, INSERT_VALUES);
    }
  }

  sol->Assembly_GhostUpdate();
}

void PNonlinear_NS_Solver::rescale_dot_inflow_value( const double &stime,
    const ALocal_InflowBC * const &infbc,
    PDNSolution * const &sol ) const
{
  const int num_nbc = infbc -> get_num_nbc();

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    const int numnode = infbc -> get_Num_LD( nbc_id );

    const double factor  = flrate -> get_flow_rate( nbc_id, stime );
    const double std_dev = flrate -> get_flow_TI_std_dev( nbc_id );

    for(int ii=0; ii<numnode; ++ii)
    {
      const int node_index = infbc -> get_LDN( nbc_id, ii );
      
      const int base_idx[3] = { node_index*4+1, node_index*4+2, node_index*4+3 };

      double base_vals[3];

      VecGetValues(sol_base->solution, 3, base_idx, base_vals);

      const double perturb_x = MATH_T::gen_double_rand_normal(0, std_dev);
      const double perturb_y = MATH_T::gen_double_rand_normal(0, std_dev);
      const double perturb_z = MATH_T::gen_double_rand_normal(0, std_dev);

      const double vals[3] = { base_vals[0] * factor * (1.0 + perturb_x), 
        base_vals[1] * factor * (1.0 + perturb_y),
        base_vals[2] * factor * (1.0 + perturb_z) };

      VecSetValues(sol->solution, 3, base_idx, vals, INSERT_VALUES);
    }
  }

  sol->Assembly_GhostUpdate();
}

// EOF
