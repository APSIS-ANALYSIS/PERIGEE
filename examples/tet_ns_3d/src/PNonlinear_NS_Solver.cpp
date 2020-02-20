#include "PNonlinear_NS_Solver.hpp"

PNonlinear_NS_Solver::PNonlinear_NS_Solver(
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol,
    const int &input_max_iteration, const int &input_renew_freq )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq)
{
  // Generate the incremental solution vector 
  dot_step = new PDNSolution_NS( anode_ptr, 0, false );
}


PNonlinear_NS_Solver::~PNonlinear_NS_Solver()
{
  delete dot_step; dot_step = nullptr;
}


void PNonlinear_NS_Solver::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "relative tolerance: %e \n", nr_tol);
  PetscPrintf(PETSC_COMM_WORLD, "absolute tolerance: %e \n", na_tol);
  PetscPrintf(PETSC_COMM_WORLD, "divergence tolerance: %e \n", nd_tol);
  PetscPrintf(PETSC_COMM_WORLD, "maximum iteration: %d \n", nmaxits);
  PetscPrintf(PETSC_COMM_WORLD, "tangent matrix renew frequency: %d \n", nrenew_freq);
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
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
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Inflow_NodalBC * const &infnbc_part,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
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

    gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, sol,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part, gbc );
   
#ifdef PETSC_USE_LOG
    PetscLogEventEnd(mat_assem_0_event,0,0,0,0);
#endif

    PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
    
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
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part, gbc );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(vec_assem_0_event,0,0,0,0);
#endif
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  PetscPrintf(PETSC_COMM_WORLD, "  Init res 2-norm: %e \n", initial_norm);

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

    // Assembly ALE-NS residual (& tangent) 
    if( nl_counter % nrenew_freq == 0 || nl_counter >= 4 )
    {
      gassem_ptr->Clear_KG();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(mat_assem_1_event, 0,0,0,0);
#endif

      gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, sol,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part, gbc );

#ifdef PETSC_USE_LOG
      PetscLogEventEnd(mat_assem_1_event,0,0,0,0);
#endif

      PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
      lsolver_ptr->SetOperator(gassem_ptr->K);
    }
    else
    {
      gassem_ptr->Clear_G();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(vec_assem_1_event, 0,0,0,0);
#endif

      gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part, gbc );

#ifdef PETSC_USE_LOG
      PetscLogEventEnd(vec_assem_1_event,0,0,0,0);
#endif
    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    PetscPrintf(PETSC_COMM_WORLD, "  --- nl_res: %e \n", residual_norm);

    relative_error = residual_norm / initial_norm;

    if( relative_error >= nd_tol )
    {
      PetscPrintf(PETSC_COMM_WORLD,
          "Warning: nonlinear solver is diverging with error %e \n",
          relative_error);
      break;
    }
  }while(nl_counter<nmaxits && relative_error > nr_tol && residual_norm > na_tol);

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol) conv_flag = true;
  else conv_flag = false;
}


void PNonlinear_NS_Solver::rescale_inflow_value( const double &stime,
    const ALocal_Inflow_NodalBC * const &infbc,
    const ICVFlowRate * const &flrate,
    const PDNSolution * const &sol_base,
    PDNSolution * const &sol ) const
{
  const int numnode = infbc -> get_Num_LD();

  const double val = flrate -> get_flow_rate( stime );

  double base_vals[3];
  int base_idx[3];

  for(int ii=0; ii<numnode; ++ii)
  {
    const int node_index = infbc -> get_LDN( ii );

    base_idx[0] = node_index * 4 + 1;
    base_idx[1] = node_index * 4 + 2;
    base_idx[2] = node_index * 4 + 3;

    VecGetValues(sol_base->solution, 3, base_idx, base_vals);

    VecSetValue(sol->solution, node_index*4+1, base_vals[0] * val, INSERT_VALUES);
    VecSetValue(sol->solution, node_index*4+2, base_vals[1] * val, INSERT_VALUES);
    VecSetValue(sol->solution, node_index*4+3, base_vals[2] * val, INSERT_VALUES);
  }

  VecAssemblyBegin(sol->solution); VecAssemblyEnd(sol->solution);
  sol->GhostUpdate();
}

// EOF
