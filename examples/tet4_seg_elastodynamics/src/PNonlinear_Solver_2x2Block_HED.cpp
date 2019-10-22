#include "PNonlinear_Solver_2x2Block_HED.hpp"

PNonlinear_Solver_2x2Block_HED::PNonlinear_Solver_2x2Block_HED( 
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr, const double &input_nrtol,
    const double &input_natol, const double &input_ndtol,
    const int &input_max_iteration, const int &input_renew_freq )
: nr_tol( input_nrtol ), na_tol( input_natol ), nd_tol( input_ndtol ),
  nmaxits( input_max_iteration ), nrenew_freq( input_renew_freq )
{
  step_velo = new PDNSolution_Disp_3D( anode_ptr, feanode_ptr, 0 );
  step_pres = new PDNSolution_Pres_3D( anode_ptr, feanode_ptr, 0 );

  sol_v_alpha = new PDNSolution_Disp_3D( anode_ptr, feanode_ptr, 0 );
  dot_v_alpha = new PDNSolution_Disp_3D( anode_ptr, feanode_ptr, 0 );

  sol_d_alpha = new PDNSolution_Disp_3D( anode_ptr, feanode_ptr, 0 );
  dot_d_alpha = new PDNSolution_Disp_3D( anode_ptr, feanode_ptr, 0 );
  
  sol_p_alpha = new PDNSolution_Pres_3D( anode_ptr, feanode_ptr, 0 );
  dot_p_alpha = new PDNSolution_Pres_3D( anode_ptr, feanode_ptr, 0 );
}


PNonlinear_Solver_2x2Block_HED::~PNonlinear_Solver_2x2Block_HED()
{
  delete step_velo; step_velo = NULL;
  delete step_pres; step_pres = NULL;

  delete sol_v_alpha; sol_v_alpha = NULL;
  delete dot_v_alpha; dot_v_alpha = NULL;
  delete sol_d_alpha; sol_d_alpha = NULL;
  delete dot_d_alpha; dot_d_alpha = NULL;
  delete sol_p_alpha; sol_p_alpha = NULL;
  delete dot_p_alpha; dot_p_alpha = NULL;
}


void PNonlinear_Solver_2x2Block_HED::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "relative tolerance: %e \n", nr_tol);
  PetscPrintf(PETSC_COMM_WORLD, "absolute tolerance: %e \n", na_tol);
  PetscPrintf(PETSC_COMM_WORLD, "divergence tolerance: %e \n", nd_tol);
  PetscPrintf(PETSC_COMM_WORLD, "maximum iteration: %d \n", nmaxits);
  PetscPrintf(PETSC_COMM_WORLD, "tangent matrix renew frequency: %d \n", nrenew_freq);
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PNonlinear_Solver_2x2Block_HED::GenAlpha_solve(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &pre_dot_d,
    const PDNSolution * const &pre_dot_p,
    const PDNSolution * const &pre_dot_v,
    const PDNSolution * const &pre_sol_d,
    const PDNSolution * const &pre_sol_p,
    const PDNSolution * const &pre_sol_v,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem_2x2Block * const &lassem_ptr,
    IPGAssem_2x2Block * const &gassem_ptr,
    IPLinear_Solver_2x2Block * const &lsolver_ptr,
    PDNSolution * const &dot_d,
    PDNSolution * const &dot_p,
    PDNSolution * const &dot_v,
    PDNSolution * const &sol_d,
    PDNSolution * const &sol_p,
    PDNSolution * const &sol_v,
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

  // Gen-alpha parameters
  const double gamma   = tmga_ptr->get_gamma();
  const double alpha_m = tmga_ptr->get_alpha_m();
  const double alpha_f = tmga_ptr->get_alpha_f();

  const double val_1 = alpha_f * gamma * dt / alpha_m;
  const double val_2 = (1.0/gamma) - (1.0/alpha_m);
  const double val_3 = 1.0 / alpha_m;

  // same-Y predictor
  sol_d -> Copy( *pre_sol_d );
  sol_p -> Copy( *pre_sol_p );
  sol_v -> Copy( *pre_sol_v );

  dot_d -> Copy( *pre_dot_d ); dot_d -> ScaleValue( (gamma-1.0)/gamma );
  dot_p -> Copy( *pre_dot_p ); dot_p -> ScaleValue( (gamma-1.0)/gamma );
  dot_v -> Copy( *pre_dot_v ); dot_v -> ScaleValue( (gamma-1.0)/gamma );

  // Define intermediate step solution
  sol_d_alpha -> Copy( *pre_sol_d );
  sol_d_alpha -> ScaleValue( 1.0 - alpha_f );
  sol_d_alpha -> PlusAX( *sol_d, alpha_f );

  sol_p_alpha -> Copy( *pre_sol_p );
  sol_p_alpha -> ScaleValue( 1.0 - alpha_f );
  sol_p_alpha -> PlusAX( *sol_p, alpha_f );

  sol_v_alpha -> Copy( *pre_sol_v );
  sol_v_alpha -> ScaleValue( 1.0 - alpha_f );
  sol_v_alpha -> PlusAX( *sol_v, alpha_f );

  dot_d_alpha -> Copy( *pre_dot_d );
  dot_d_alpha -> ScaleValue( 1.0 - alpha_m );
  dot_d_alpha -> PlusAX( *dot_d, alpha_m );

  dot_p_alpha -> Copy( *pre_dot_p );
  dot_p_alpha -> ScaleValue( 1.0 - alpha_m );
  dot_p_alpha -> PlusAX( *dot_p, alpha_m );

  dot_v_alpha -> Copy( *pre_dot_v );
  dot_v_alpha -> ScaleValue( 1.0 - alpha_m );
  dot_v_alpha -> PlusAX( *dot_v, alpha_m );

  // step_disp = -1/alpha_m R^k_(1), we use step_velo to hold step_disp here
  step_velo -> Copy( pre_dot_d );
  step_velo -> ScaleValue( val_2 );
  step_velo -> PlusAX( pre_sol_v, val_3 );

  dot_d -> PlusAX( step_velo, 1.0 );
  sol_d -> PlusAX( step_velo, gamma * dt );

  dot_d_alpha -> PlusAX( step_velo, alpha_m );
  sol_d_alpha -> PlusAX( step_velo, alpha_f * gamma * dt );

  if( new_tangent_flag )
  {
    gassem_ptr -> Clear_KG();

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(mat_assem_0_event, 0,0,0,0);
#endif
    gassem_ptr -> Assem_tangent_residual( dot_d_alpha, dot_p_alpha,
        dot_v_alpha, sol_d_alpha, sol_p_alpha, sol_v_alpha,
        curr_time, dt, lassem_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, feanode_ptr, nbc_part, ebc_part );
#ifdef PETSC_USE_LOG
    PetscLogEventEnd(mat_assem_0_event,0,0,0,0);
#endif 
    // optional, Symmetric Jacobi Scaling
    lsolver_ptr -> SymmJacobi_MatVec_Scale( gassem_ptr );

    // Destroy previously generated Schur complement object
    //lsolver_ptr -> DestroyOperators();

    // Set Lhs and Rhs
    lsolver_ptr -> SetLHS_RHS( gassem_ptr );

    PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
  }
  else
  {
    gassem_ptr -> Clear_G();

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(vec_assem_0_event, 0,0,0,0);
#endif
    gassem_ptr -> Assem_residual( dot_d_alpha, dot_p_alpha, dot_v_alpha,
        sol_d_alpha, sol_p_alpha, sol_v_alpha, curr_time, dt,
        lassem_ptr, elementv, elements, quad_v, quad_s, lien_ptr,
        feanode_ptr, nbc_part, ebc_part );
#ifdef PETSC_USE_LOG
    PetscLogEventEnd(vec_assem_0_event,0,0,0,0);
#endif 

    // Optional, Symmetric Jacobi Scaling
    lsolver_ptr -> SymmJacobi_Vec_Scale( gassem_ptr );

    // Set Rhs
    lsolver_ptr -> SetRHS( gassem_ptr );
  }

  // Initialization
  nl_counter = 0;
  double residual_norm = 0.0, relative_error = 0.0;
  double res_0 = 0.0, res_1 = 0.0;

  VecNorm( gassem_ptr->G_0, NORM_2, &res_0 );
  VecNorm( gassem_ptr->G_1, NORM_2, &res_1 );
  const double initial_norm = std::sqrt(res_0 * res_0 + res_1 * res_1);
  PetscPrintf(PETSC_COMM_WORLD, "  Init res 2-norm: %e \n", initial_norm);

  // Consistent Newton-Raphson iteration
  do
  {
    // Call solver
    lsolver_ptr -> Solve_RHS( gassem_ptr, step_pres, step_velo );

    nl_counter += 1;

    dot_d -> PlusAX( step_velo, -1.0 * val_1 );
    dot_p -> PlusAX( step_pres, -1.0 );
    dot_v -> PlusAX( step_velo, -1.0 );

    sol_d -> PlusAX( step_velo, -1.0 * val_1 * gamma * dt );
    sol_p -> PlusAX( step_pres, -1.0 * gamma * dt );
    sol_v -> PlusAX( step_velo, -1.0 * gamma * dt );

    dot_d_alpha -> PlusAX( step_velo, -1.0 * val_1 * alpha_m );
    dot_p_alpha -> PlusAX( step_pres, -1.0 * alpha_m );
    dot_v_alpha -> PlusAX( step_velo, -1.0 * alpha_m );

    sol_d_alpha -> PlusAX( step_velo, -1.0 * val_1 * gamma * dt * alpha_f );
    sol_p_alpha -> PlusAX( step_pres, -1.0 * gamma * dt * alpha_f );
    sol_v_alpha -> PlusAX( step_velo, -1.0 * gamma * dt * alpha_f );

    if( nl_counter >= nrenew_freq )
    {
      gassem_ptr -> Clear_KG();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(mat_assem_1_event, 0,0,0,0);
#endif
      gassem_ptr -> Assem_tangent_residual( dot_d_alpha, dot_p_alpha,
          dot_v_alpha, sol_d_alpha, sol_p_alpha, sol_v_alpha,
          curr_time, dt, lassem_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, feanode_ptr, nbc_part, ebc_part );
#ifdef PETSC_USE_LOG
      PetscLogEventEnd(mat_assem_1_event,0,0,0,0);
#endif 
     
      // optional, Symmetric Jacobi Scaling
      lsolver_ptr -> SymmJacobi_MatVec_Scale( gassem_ptr );

      // Destroy previously generated Schur complement object
      //lsolver_ptr -> DestroyOperators();

      // Set Lhs and Rhs
      lsolver_ptr -> SetLHS_RHS( gassem_ptr );
      
      PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
    }
    else
    {
      gassem_ptr -> Clear_G();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(vec_assem_1_event, 0,0,0,0);
#endif
      gassem_ptr -> Assem_residual( dot_d_alpha, dot_p_alpha, dot_v_alpha,
          sol_d_alpha, sol_p_alpha, sol_v_alpha, curr_time, dt,
          lassem_ptr, elementv, elements, quad_v, quad_s, lien_ptr,
          feanode_ptr, nbc_part, ebc_part );
#ifdef PETSC_USE_LOG
      PetscLogEventEnd(vec_assem_1_event,0,0,0,0);
#endif 

      // Optional, Symmetric Jacobi Scaling
      lsolver_ptr -> SymmJacobi_Vec_Scale( gassem_ptr );

      // Set Rhs
      lsolver_ptr -> SetRHS( gassem_ptr );
    }

    VecNorm( gassem_ptr->G_0, NORM_2, &res_0 );
    VecNorm( gassem_ptr->G_1, NORM_2, &res_1 );
    
    residual_norm = std::sqrt(res_0 * res_0 + res_1 * res_1);
    
    PetscPrintf(PETSC_COMM_WORLD, "  --- res norm: res0 %e res1 %e total %e \n", 
        res_0, res_1, residual_norm);

    relative_error = residual_norm / initial_norm;

    SYS_T::print_fatal_if( relative_error >= nd_tol, "Error: Nonlinear solver is diverging.\n");
  }while(nl_counter < nmaxits && relative_error > nr_tol && residual_norm > na_tol);

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol) conv_flag = true;
  else conv_flag = false;
}


// EOF
