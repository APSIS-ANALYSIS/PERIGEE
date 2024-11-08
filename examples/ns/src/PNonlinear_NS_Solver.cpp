#include "PNonlinear_NS_Solver.hpp"

PNonlinear_NS_Solver::PNonlinear_NS_Solver(
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol,
    const int &input_max_iteration, 
    const int &input_renew_freq,
    const int &input_renew_threshold )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq),
  nrenew_threshold(input_renew_threshold)
{
  // Generate the incremental solution vector used for update 
  // the solution of the nonlinear algebraic system 
  dot_step = new PDNSolution_NS( anode_ptr, 0, false );
}


PNonlinear_NS_Solver::~PNonlinear_NS_Solver()
{
  delete dot_step; dot_step = nullptr;
}


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
    const PDNSolution * const &sol_base,
    const PDNSolution * const &pre_dot_sol,
    const PDNSolution * const &pre_sol,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ICVFlowRate * const flr_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_InflowBC * const &infnbc_part,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc,
    const ALocal_WeakBC * const &wbc_part,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    FEAElement * const &elementvs,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PDNSolution * const &dot_sol,
    PDNSolution * const &sol,
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
  const double gamma   = tmga_ptr->get_gamma();
  const double alpha_m = tmga_ptr->get_alpha_m();
  const double alpha_f = tmga_ptr->get_alpha_f();

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
  rescale_inflow_value(curr_time+dt, infnbc_part, flr_ptr, sol_base, sol);
  rescale_inflow_value(curr_time+alpha_f*dt, infnbc_part, flr_ptr, sol_base, &sol_alpha);
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
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements, elementvs,
        quad_v, quad_s, lien_ptr, feanode_ptr, nbc_part, ebc_part, gbc, wbc_part );
   
#ifdef PETSC_USE_LOG
    PetscLogEventEnd(mat_assem_0_event,0,0,0,0);
#endif

    SYS_T::commPrint("  --- M updated");
    
    // SetOperator will pass the tangent matrix to the linear solver and the
    // linear solver will generate the preconditioner based on the new matrix.
    lsolver_ptr->SetOperator( gassem_ptr->K );
  }
  else
  {
    gassem_ptr->Clear_G();

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(vec_assem_0_event, 0,0,0,0);
#endif

    gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements, elementvs,
        quad_v, quad_s, lien_ptr, feanode_ptr, nbc_part, ebc_part, gbc, wbc_part );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(vec_assem_0_event,0,0,0,0);
#endif
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  SYS_T::commPrint("  Init res 2-norm: %e \n", initial_norm);

  // Now do consistent Newton-Raphson iteration
  do
  {
#ifdef PETSC_USE_LOG
    PetscLogEventBegin(lin_solve_event, 0,0,0,0);
#endif
    
    // solve the equation K dot_step = G
    lsolver_ptr->Solve( gassem_ptr->G, dot_step );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(lin_solve_event,0,0,0,0);
#endif

    bc_mat->MatMultSol( dot_step );

    nl_counter += 1;

    dot_sol->PlusAX( dot_step, -1.0 );
    sol->PlusAX( dot_step, (-1.0) * gamma * dt );

    dot_sol_alpha.PlusAX( dot_step, (-1.0) * alpha_m );
    sol_alpha.PlusAX( dot_step, (-1.0) * alpha_f * gamma * dt );

    // Assembly residual (& tangent if condition satisfied) 
    if( nl_counter % nrenew_freq == 0 || nl_counter >= nrenew_threshold )
    {
      gassem_ptr->Clear_KG();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(mat_assem_1_event, 0,0,0,0);
#endif

      gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements, elementvs,
          quad_v, quad_s, lien_ptr, feanode_ptr, nbc_part, ebc_part, gbc, wbc_part );

#ifdef PETSC_USE_LOG
      PetscLogEventEnd(mat_assem_1_event,0,0,0,0);
#endif

      SYS_T::commPrint("  --- M updated");
      lsolver_ptr->SetOperator(gassem_ptr->K);
    }
    else
    {
      gassem_ptr->Clear_G();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(vec_assem_1_event, 0,0,0,0);
#endif

      gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements, elementvs,
          quad_v, quad_s, lien_ptr, feanode_ptr, nbc_part, ebc_part, gbc, wbc_part );

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

void PNonlinear_NS_Solver::HERK_Solve_NS(
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &sol_base,
    const PDNSolution * const &dot_sol_base,
    PDNSolution ** const &cur_velo_sols,
    PDNSolution * const &cur_velo,
    PDNSolution * const &cur_dot_velo,
    PDNSolution ** const &cur_pres_sols,
    PDNSolution * const &cur_pres,
    PDNSolution ** const &pre_velo_sols,
    PDNSolution * const &pre_velo,
    PDNSolution ** const &pre_pres_sols,
    PDNSolution * const &pre_pres,
    PDNSolution * const &pre_velo_before,
    const Runge_Kutta_Butcher * const &tm_RK_ptr,
    const ICVFlowRate * const flr_ptr,
    const ICVFlowRate * const dot_flr_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_InflowBC * const &infnbc_part,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc,
    const ALocal_WeakBC * const &wbc_part,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    FEAElement * const &elementvs,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PDNSolution * const &cur_sol) const
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

  // HERK steps
  const int ss = tm_RK_ptr->get_RK_step();

  // 第一个子步, u_1 = u_n
  cur_velo_sols[0] -> Copy(*pre_velo);
  
  SYS_T::commPrint(" --- substep = 1 is solved. \n");

  // 子步（从第二个子步开始）
  for(int ii = 1; ii < ss; ++ii)
  {
    // 使得每个子步的速度满足Dirchlet边界
    rescale_inflow_velo(curr_time + tm_RK_ptr->get_RK_c(ii) * dt, infnbc_part, flr_ptr, sol_base, cur_velo_sols[ii]);

    gassem_ptr->Clear_KG();
    
    gassem_ptr->Assem_tangent_residual_substep( ii, cur_velo_sols, cur_pres_sols,
      pre_velo_sols, pre_velo, pre_pres_sols, pre_velo_before, tm_RK_ptr, 
      curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements, elementvs,
      quad_v, quad_s, lien_ptr, feanode_ptr, nbc_part, ebc_part, gbc, wbc_part );     
    
    lsolver_ptr->SetOperator(gassem_ptr->K);
    lsolver_ptr->Solve(gassem_ptr->G, dot_step);

    SYS_T::commPrint(" --- substep = %d is solved. \n", ii);

    Update_pressure_velocity(anode_ptr, cur_velo_sols[ii], cur_pres_sols[ii-1], dot_step);
  }
  // 终步
    // 使得终步的速度满足Dirchlet边界
    rescale_inflow_velo(curr_time + dt, infnbc_part, flr_ptr, sol_base, cur_velo);

    gassem_ptr->Clear_KG();

    gassem_ptr->Assem_tangent_residual_laststep( cur_velo_sols, cur_velo, 
      cur_pres_sols, pre_velo_sols, pre_velo, pre_pres_sols, pre_velo_before,
      tm_RK_ptr, curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements, elementvs,
      quad_v, quad_s, lien_ptr, feanode_ptr, nbc_part, ebc_part, gbc, wbc_part );

    lsolver_ptr->SetOperator(gassem_ptr->K);
    lsolver_ptr->Solve(gassem_ptr->G, dot_step);

    SYS_T::commPrint(" --- laststep is solved. \n");

    Update_pressure_velocity(anode_ptr, cur_velo, cur_pres_sols[ss-1], dot_step);

  // 最终步
    // 使得终步的dot速度满足Dirchlet边界
    rescale_inflow_velo(curr_time + dt, infnbc_part, dot_flr_ptr, dot_sol_base, cur_dot_velo);

    gassem_ptr->Clear_KG();

    gassem_ptr->Assem_tangent_residual_finalstep( cur_dot_velo, cur_velo_sols, 
      cur_velo, cur_pres_sols, pre_velo, cur_pres, tm_RK_ptr, curr_time, dt, 
      alelem_ptr, lassem_ptr, elementv, elements, elementvs, quad_v, quad_s, 
      lien_ptr, feanode_ptr, nbc_part, ebc_part, gbc, wbc_part );

    lsolver_ptr->SetOperator(gassem_ptr->K);
    lsolver_ptr->Solve(gassem_ptr->G, dot_step);

    SYS_T::commPrint(" --- finalstep is solved. \n");

    Update_pressure_velocity(anode_ptr, cur_dot_velo, cur_pres, dot_step);

  // 将n+1步速度和压强组装为解向量
    Update_solutions(anode_ptr, cur_velo, cur_pres, cur_sol);
}

void PNonlinear_NS_Solver::rescale_inflow_value( const double &stime,
    const ALocal_InflowBC * const &infbc,
    const ICVFlowRate * const &flrate,
    const PDNSolution * const &sol_base,
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

void PNonlinear_NS_Solver::rescale_inflow_velo( const double &stime,
    const ALocal_InflowBC * const &infbc,
    const ICVFlowRate * const &flrate,
    const PDNSolution * const &base,
    PDNSolution * const &velo ) const
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

      VecGetValues(base->solution, 3, base_idx, base_vals);

      const double perturb_x = MATH_T::gen_double_rand_normal(0, std_dev);
      const double perturb_y = MATH_T::gen_double_rand_normal(0, std_dev);
      const double perturb_z = MATH_T::gen_double_rand_normal(0, std_dev);

      const double vals[3] = { base_vals[0] * factor * (1.0 + perturb_x), 
        base_vals[1] * factor * (1.0 + perturb_y),
        base_vals[2] * factor * (1.0 + perturb_z) };

      const int velo_idx[3] = { node_index*3, node_index*3+1, node_index*3+2 };

      VecSetValues(velo->solution, 3, velo_idx, vals, INSERT_VALUES);
    }
  }

  velo->Assembly_GhostUpdate();
}

void PNonlinear_NS_Solver::Update_pressure_velocity(     
    const APart_Node * const &anode_ptr, 
    PDNSolution * const &velo,
    PDNSolution * const &pres,
    const PDNSolution * const &step) const
  {
    Vec lvelo, lpres, lstep;
    double * array_velo, * array_pres, * array_step;

    VecGhostGetLocalForm(velo->solution, &lvelo);    
    VecGhostGetLocalForm(pres->solution, &lpres);    
    VecGhostGetLocalForm(step->solution, &lstep);

    VecGetArray(lvelo, &array_velo);
    VecGetArray(lpres, &array_pres); 
    VecGetArray(lstep, &array_step);

    for(int ii=0; ii<anode_ptr->get_nlocalnode(); ++ii)
    {
      array_pres[ii       ] = array_step[ii*4 + 0];
      array_velo[ii*3 + 0 ] = array_step[ii*4 + 1];
      array_velo[ii*3 + 1 ] = array_step[ii*4 + 2];
      array_velo[ii*3 + 2 ] = array_step[ii*4 + 3];
    }

    VecRestoreArray(lvelo, &array_velo);    
    VecRestoreArray(lpres, &array_pres);

    VecGhostRestoreLocalForm(velo->solution, &lvelo);
    VecGhostRestoreLocalForm(pres->solution, &lpres);
    
    velo->GhostUpdate();
    pres->GhostUpdate();  
  }

void PNonlinear_NS_Solver::Update_solutions(     
    const APart_Node * const &anode_ptr, 
    const PDNSolution * const &velo,
    const PDNSolution * const &pres,
    PDNSolution * const &sol) const
  {
    Vec lvelo, lpres, lsol;
    double * array_velo, * array_pres, * array_sol;

    VecGhostGetLocalForm(velo->solution, &lvelo);    
    VecGhostGetLocalForm(pres->solution, &lpres);    
    VecGhostGetLocalForm(sol->solution, &lsol);

    VecGetArray(lvelo, &array_velo);
    VecGetArray(lpres, &array_pres); 
    VecGetArray(lsol, &array_sol);

    for(int ii=0; ii<anode_ptr->get_nlocalnode(); ++ii)
    {
      array_sol[ii*4 + 0] = array_pres[ii       ];
      array_sol[ii*4 + 1] = array_velo[ii*3 + 0 ];
      array_sol[ii*4 + 2] = array_velo[ii*3 + 1 ];
      array_sol[ii*4 + 3] = array_velo[ii*3 + 2 ];
    }

    VecRestoreArray(lsol, &array_sol);    

    VecGhostRestoreLocalForm(sol->solution, &lsol);
    
    sol->GhostUpdate();
  }
  
// EOF
