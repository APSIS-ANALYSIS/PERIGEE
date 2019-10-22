#include "PNonlinear_Seg_Solver.hpp"

PNonlinear_Seg_Solver::PNonlinear_Seg_Solver(
    const APart_Node * const &anode_ptr,
    const IPLocAssem * const &lassem_ptr,
    const FEANode * const &feanode_ptr,
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol,
    const int &input_max_iteration, const int &input_renew_freq )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq)
{
  // Generate the solution vector for P-V system as the step updater
  dot_P_V_step = new PDNSolution_P_V_Mixed_Hyperelastic_3D( 
      anode_ptr, lassem_ptr, feanode_ptr, 0 );
}


PNonlinear_Seg_Solver::~PNonlinear_Seg_Solver()
{
  delete dot_P_V_step; dot_P_V_step = NULL;
}


void PNonlinear_Seg_Solver::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "relative tolerance: %e \n", nr_tol);
  PetscPrintf(PETSC_COMM_WORLD, "absolute tolerance: %e \n", na_tol);
  PetscPrintf(PETSC_COMM_WORLD, "divergence tolerance: %e \n", nd_tol);
  PetscPrintf(PETSC_COMM_WORLD, "maximum iteration: %d \n", nmaxits);
  PetscPrintf(PETSC_COMM_WORLD, "tangent matrix renew frequency: %d \n", nrenew_freq);
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PNonlinear_Seg_Solver::GenAlpha_Seg_solve(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &pre_dot_sol,
    const PDNSolution * const &pre_sol,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
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
  // Initialize the counter and error
  nl_counter = 0;
  double residual_norm = 0.0, initial_norm = 0.0, relative_error = 0.0;

  // Gen-alpha parameters
  const double gamma   = tmga_ptr->get_gamma();
  const double alpha_m = tmga_ptr->get_alpha_m();
  const double alpha_f = tmga_ptr->get_alpha_f();

  const double val_1 = alpha_f * gamma * dt / alpha_m;
  const double val_2 = (1.0/gamma) - (1.0/alpha_m);
  const double val_3 = 1.0 / alpha_m;

  // Same-Y predictor
  sol->Copy(*pre_sol);
  dot_sol->Copy(*pre_dot_sol);
  dot_sol->ScaleValue( (gamma-1.0)/gamma );

  // Generate the U-P-V step updater
  PDNSolution dot_step(*dot_sol);

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
    gassem_ptr->Clear_KG();
    gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part );
    PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
    lsolver_ptr->SetOperator(gassem_ptr->K);
  }
  else
  {
    gassem_ptr->Clear_G();
    gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part );
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  PetscPrintf(PETSC_COMM_WORLD, "  Init res 2-norm: %e \n", initial_norm);

  // Solve for the first round
  lsolver_ptr->Solve( gassem_ptr->G, dot_P_V_step );
  PetscPrintf(PETSC_COMM_WORLD, " --> update displacement. \n");
  bc_mat->MatMultSol( dot_P_V_step ); // strictly enforce essential BC

  // Get the dot_step from the dot_P_V_step
  dot_step.ScaleValue(0.0);
  SEG_SOL_T::PlusAiPV( val_1*(-1.0), -1.0, -1.0, dot_P_V_step, &dot_step );
  SEG_SOL_T::PlusAiUPV( val_2, 0.0, 0.0, pre_dot_sol, &dot_step );
  SEG_SOL_T::PlusAiVPV( val_3, 0.0, 0.0, pre_sol, &dot_step );

  // Now use dot_step to udpate solutions
  dot_sol->PlusAX( dot_step, 1.0 );
  sol->PlusAX( dot_step, gamma * dt );

  dot_sol_alpha.PlusAX( dot_step, alpha_m );
  sol_alpha.PlusAX( dot_step, alpha_f * gamma * dt );

  nl_counter += 1;
  
  // Now do consistent Newton-Raphson iteration
  do
  {
    if( nl_counter >= nrenew_freq )
    {
      gassem_ptr->Clear_KG();
      gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part );
      PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
      lsolver_ptr->SetOperator(gassem_ptr->K);
    }
    else
    {
      gassem_ptr->Clear_G();
      gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part );
    }
    
    lsolver_ptr->Solve( gassem_ptr->G, dot_P_V_step );
    bc_mat->MatMultSol( dot_P_V_step );

    nl_counter += 1;

    dot_step.ScaleValue(0.0);
    SEG_SOL_T::PlusAiPV( val_1 * (-1.0), -1.0, -1.0, dot_P_V_step, &dot_step );

    dot_sol->PlusAX( dot_step, 1.0 );
    sol->PlusAX( dot_step, gamma * dt );

    dot_sol_alpha.PlusAX( dot_step, alpha_m );
    sol_alpha.PlusAX( dot_step, alpha_f * gamma * dt );

    // CHECK MY PROPOSITION OF DOT U and V relation
    //SEG_SOL_T::CheckUV(&dot_sol_alpha, &sol_alpha);
    // Finish my check routine

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    PetscPrintf(PETSC_COMM_WORLD, "  --- res norm: %e \n", residual_norm);

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


void PNonlinear_Seg_Solver::GenAlpha_Seg_solve_2(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &pre_dot_sol,
    const PDNSolution * const &pre_sol,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
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
  PetscClassId classid_assembly;
  PetscClassIdRegister("mat_vec_assembly", &classid_assembly);
  PetscLogEventRegister("assembly mat 0", classid_assembly, &mat_assem_0_event);
  PetscLogEventRegister("assembly mat 1", classid_assembly, &mat_assem_1_event);
  PetscLogEventRegister("assembly vec 0", classid_assembly, &vec_assem_0_event);
  PetscLogEventRegister("assembly vec 1", classid_assembly, &vec_assem_1_event);
#endif

  // Initialize the counter and error
  nl_counter = 0;
  double residual_norm = 0.0, initial_norm = 0.0, relative_error = 0.0;

  // Gen-alpha parameters
  const double gamma   = tmga_ptr->get_gamma();
  const double alpha_m = tmga_ptr->get_alpha_m();
  const double alpha_f = tmga_ptr->get_alpha_f();

  const double val_1 = alpha_f * gamma * dt / alpha_m;
  const double val_2 = (1.0/gamma) - (1.0/alpha_m);
  const double val_3 = 1.0 / alpha_m;

  // Same-Y predictor
  sol->Copy(*pre_sol);
  dot_sol->Copy(*pre_dot_sol);
  dot_sol->ScaleValue( (gamma-1.0)/gamma );

  // Generate the U-P-V step updater
  PDNSolution dot_step(*dot_sol);

  // Define the dol_sol at alpha_m: dot_sol_alpha
  PDNSolution dot_sol_alpha(*pre_dot_sol);
  dot_sol_alpha.ScaleValue( 1.0 - alpha_m );
  dot_sol_alpha.PlusAX(*dot_sol, alpha_m);

  // Define the sol at alpha_f: sol_alpha
  PDNSolution sol_alpha(*pre_sol);
  sol_alpha.ScaleValue( 1.0 - alpha_f );
  sol_alpha.PlusAX( *sol, alpha_f );

  // Get the dot_step from the dot_P_V_step
  dot_step.ScaleValue(0.0);
  SEG_SOL_T::PlusAiUPV( val_2, 0.0, 0.0, pre_dot_sol, &dot_step );
  SEG_SOL_T::PlusAiVPV( val_3, 0.0, 0.0, pre_sol, &dot_step );

  // Now use dot_step to udpate solutions
  dot_sol->PlusAX( dot_step, 1.0 );
  sol->PlusAX( dot_step, gamma * dt );

  dot_sol_alpha.PlusAX( dot_step, alpha_m );
  sol_alpha.PlusAX( dot_step, alpha_f * gamma * dt );

  // If new_tangent_flag == TRUE, update the tangent matrix;
  // otherwise, use the matrix from the previous time step
  if( new_tangent_flag )
  {
    gassem_ptr->Clear_KG();

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(mat_assem_0_event, 0,0,0,0);
#endif
    gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part );
#ifdef PETSC_USE_LOG
    PetscLogEventEnd(mat_assem_0_event,0,0,0,0);
#endif
    PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
    lsolver_ptr->SetOperator(gassem_ptr->K);
  }
  else
  {
    gassem_ptr->Clear_G();
    
#ifdef PETSC_USE_LOG
    PetscLogEventBegin(vec_assem_0_event, 0,0,0,0);
#endif
    gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part );
#ifdef PETSC_USE_LOG
    PetscLogEventEnd(vec_assem_0_event,0,0,0,0);
#endif
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  PetscPrintf(PETSC_COMM_WORLD, "  Init res 2-norm: %e \n", initial_norm);

  // Now do consistent Newton-Raphson iteration
  do
  {
    lsolver_ptr->Solve( gassem_ptr->G, dot_P_V_step );
    bc_mat->MatMultSol( dot_P_V_step );

    nl_counter += 1;

    dot_step.ScaleValue(0.0);
    SEG_SOL_T::PlusAiPV( val_1 * (-1.0), -1.0, -1.0, dot_P_V_step, &dot_step );

    dot_sol->PlusAX( dot_step, 1.0 );
    sol->PlusAX( dot_step, gamma * dt );

    dot_sol_alpha.PlusAX( dot_step, alpha_m );
    sol_alpha.PlusAX( dot_step, alpha_f * gamma * dt );

    // CHECK MY PROPOSITION OF DOT U and V relation
    //SEG_SOL_T::CheckUV(&dot_sol_alpha, &sol_alpha);
    // Finish my check routine

    if( nl_counter >= nrenew_freq )
    {
      gassem_ptr->Clear_KG();
      
#ifdef PETSC_USE_LOG
      PetscLogEventBegin(mat_assem_1_event, 0,0,0,0);
#endif
      gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part );
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
      gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part );
#ifdef PETSC_USE_LOG
      PetscLogEventEnd(vec_assem_1_event,0,0,0,0);
#endif
    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    PetscPrintf(PETSC_COMM_WORLD, "  --- res norm: %e \n", residual_norm);

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



void PNonlinear_Seg_Solver::GenAlpha_Seg_solve_DiagScale(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &pre_dot_sol,
    const PDNSolution * const &pre_sol,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_DiagScale * const &lsolver_ptr,
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

  const double val_1 = alpha_f * gamma * dt / alpha_m;
  const double val_2 = (1.0/gamma) - (1.0/alpha_m);
  const double val_3 = 1.0 / alpha_m;

  // Same-Y predictor
  sol->Copy(*pre_sol);
  dot_sol->Copy(*pre_dot_sol);
  dot_sol->ScaleValue( (gamma-1.0)/gamma );

  // Generate the U-P-V step updater
  PDNSolution dot_step(*dot_sol);

  // Define the dol_sol at alpha_m: dot_sol_alpha
  PDNSolution dot_sol_alpha(*pre_dot_sol);
  dot_sol_alpha.ScaleValue( 1.0 - alpha_m );
  dot_sol_alpha.PlusAX(*dot_sol, alpha_m);

  // Define the sol at alpha_f: sol_alpha
  PDNSolution sol_alpha(*pre_sol);
  sol_alpha.ScaleValue( 1.0 - alpha_f );
  sol_alpha.PlusAX( *sol, alpha_f );

  // Get the dot_step from the dot_P_V_step
  dot_step.ScaleValue(0.0);
  SEG_SOL_T::PlusAiUPV( val_2, 0.0, 0.0, pre_dot_sol, &dot_step );
  SEG_SOL_T::PlusAiVPV( val_3, 0.0, 0.0, pre_sol, &dot_step );

  // Now use dot_step to udpate solutions
  dot_sol->PlusAX( dot_step, 1.0 );
  sol->PlusAX( dot_step, gamma * dt );

  dot_sol_alpha.PlusAX( dot_step, alpha_m );
  sol_alpha.PlusAX( dot_step, alpha_f * gamma * dt );

  // If new_tangent_flag == TRUE, update the tangent matrix;
  // otherwise, use the matrix from the previous time step
  if( new_tangent_flag )
  {
    gassem_ptr->Clear_KG();

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(mat_assem_0_event, 0,0,0,0);
#endif
    gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part );
#ifdef PETSC_USE_LOG
    PetscLogEventEnd(mat_assem_0_event,0,0,0,0);
#endif
    
    // Symmetric scaling
    lsolver_ptr->SymmJacobi_MatVec_Scale( gassem_ptr );
    
    lsolver_ptr->SetOperator(gassem_ptr->K);

    PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
  }
  else
  {
    gassem_ptr->Clear_G();
    
#ifdef PETSC_USE_LOG
    PetscLogEventBegin(vec_assem_0_event, 0,0,0,0);
#endif
    gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part );
#ifdef PETSC_USE_LOG
    PetscLogEventEnd(vec_assem_0_event,0,0,0,0);
#endif

    // Scale the G vector
    lsolver_ptr -> SymmJacobi_Vec_Scale( gassem_ptr );
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  PetscPrintf(PETSC_COMM_WORLD, "  Init res 2-norm: %e \n", initial_norm);

  // Now do consistent Newton-Raphson iteration
  do
  {
#ifdef PETSC_USE_LOG
    PetscLogEventBegin(lin_solve_event, 0,0,0,0);
#endif
    lsolver_ptr->Solve( gassem_ptr->G, dot_P_V_step );
#ifdef PETSC_USE_LOG
    PetscLogEventEnd(lin_solve_event,0,0,0,0);
#endif
    
    bc_mat->MatMultSol( dot_P_V_step );

    nl_counter += 1;

    dot_step.ScaleValue(0.0);
    SEG_SOL_T::PlusAiPV( val_1 * (-1.0), -1.0, -1.0, dot_P_V_step, &dot_step );

    dot_sol->PlusAX( dot_step, 1.0 );
    sol->PlusAX( dot_step, gamma * dt );

    dot_sol_alpha.PlusAX( dot_step, alpha_m );
    sol_alpha.PlusAX( dot_step, alpha_f * gamma * dt );

    // CHECK MY PROPOSITION OF DOT U and V relation
    //SEG_SOL_T::CheckUV(&dot_sol_alpha, &sol_alpha);
    // Finish my check routine

    if( nl_counter >= nrenew_freq )
    {
      gassem_ptr->Clear_KG();
      
#ifdef PETSC_USE_LOG
      PetscLogEventBegin(mat_assem_1_event, 0,0,0,0);
#endif
      gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part );
#ifdef PETSC_USE_LOG
      PetscLogEventEnd(mat_assem_1_event,0,0,0,0);
#endif
      // Symmetric scaling
      lsolver_ptr->SymmJacobi_MatVec_Scale( gassem_ptr );

      lsolver_ptr->SetOperator(gassem_ptr->K);
      
      PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
    }
    else
    {
      gassem_ptr->Clear_G();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(vec_assem_1_event, 0,0,0,0);
#endif
      gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part );
#ifdef PETSC_USE_LOG
      PetscLogEventEnd(vec_assem_1_event,0,0,0,0);
#endif

      // Scale the G vector
      lsolver_ptr -> SymmJacobi_Vec_Scale( gassem_ptr );
    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    PetscPrintf(PETSC_COMM_WORLD, "  --- res norm: %e \n", residual_norm);

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

// EOF
