#include "PNonlinear_Tet4_NS_3D_Solver.hpp"

PNonlinear_Tet4_NS_3D_Solver::PNonlinear_Tet4_NS_3D_Solver( 
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const double &input_nrtol,
    const double &input_natol, const double &input_ndtol,
    const int &input_max_iteration, const int &input_renew_freq )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq)
{
  // Generate a zero solution vector for storing the Newton update
  velo_step = new PDNSolution_Tet4_NS_3D( anode_ptr, feanode_ptr, 0, false );

  dot_sol_alpha = new PDNSolution_Tet4_NS_3D( anode_ptr, feanode_ptr, 0, false );
  
  sol_alpha = new PDNSolution_Tet4_NS_3D( anode_ptr, feanode_ptr, 0, false );
}


PNonlinear_Tet4_NS_3D_Solver::~PNonlinear_Tet4_NS_3D_Solver()
{
  delete velo_step; velo_step = NULL;
  delete dot_sol_alpha; dot_sol_alpha = NULL;
  delete sol_alpha; sol_alpha = NULL;
}


void PNonlinear_Tet4_NS_3D_Solver::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "relative tolerance: %e \n", nr_tol);
  PetscPrintf(PETSC_COMM_WORLD, "absolute tolerance: %e \n", na_tol);
  PetscPrintf(PETSC_COMM_WORLD, "divergence tolerance: %e \n", nd_tol);
  PetscPrintf(PETSC_COMM_WORLD, "maximum iteration: %d \n", nmaxits);
  PetscPrintf(PETSC_COMM_WORLD, "tangent matrix renew frequency: %d \n", nrenew_freq);
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PNonlinear_Tet4_NS_3D_Solver::GenAlpha_VMS_solve(
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
  nl_counter = 0;
  double residual_norm = 0.0, initial_norm = 0.0, relative_error = 0.0;

  // Gen-alpha parameters
  const double gamma   = tmga_ptr->get_gamma();
  const double alpha_m = tmga_ptr->get_alpha_m();
  const double alpha_f = tmga_ptr->get_alpha_f();

  // predictor
  sol -> ScaleValue(0.0);
  dot_sol -> ScaleValue(0.0);

  sol     -> PlusAX(*pre_sol, 1.0);
  dot_sol -> PlusAX(*pre_dot_sol, (gamma-1.0)/gamma);

  // define alpha_f dot_sol & alpha_m sol
  sol_alpha->Copy(pre_sol);
  
  dot_sol_alpha->Copy(pre_dot_sol);
  dot_sol_alpha->ScaleValue(1.0 - alpha_m);
  dot_sol_alpha->PlusAX(dot_sol, alpha_m);
  
  if( new_tangent_flag )
  {
    gassem_ptr->Clear_KG();
    gassem_ptr->Assem_tangent_residual( dot_sol_alpha, sol_alpha,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part );
    PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
    lsolver_ptr->SetOperator(gassem_ptr->K);
  }
  else
  {
    gassem_ptr->Clear_G();
    gassem_ptr->Assem_residual( dot_sol_alpha, sol_alpha,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part );
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  PetscPrintf(PETSC_COMM_WORLD, "  Init res 2-norm: %e \n", initial_norm);

  do
  {
    lsolver_ptr->Solve( gassem_ptr->G, velo_step );

    nl_counter += 1;

    // corrector
    sol->PlusAiX(*velo_step, -1.0, -1.0*gamma*dt, 1, 3);
    dot_sol->PlusAX( velo_step, -1.0 );

    sol_alpha->PlusAiX(*velo_step, -1.0, -1.0 * alpha_f * gamma * dt, 1, 3);
    dot_sol_alpha->PlusAX(velo_step, -1.0 * alpha_m);

    if( nl_counter >= nrenew_freq )
    {
      gassem_ptr->Clear_KG();
      gassem_ptr->Assem_tangent_residual( dot_sol_alpha, sol_alpha,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part );
      PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
      lsolver_ptr->SetOperator(gassem_ptr->K);
    }
    else
    {
      gassem_ptr->Clear_G();
      gassem_ptr->Assem_residual( dot_sol_alpha, sol_alpha,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part );
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
  }while(nl_counter < nmaxits && relative_error > nr_tol &&
      residual_norm > na_tol );

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol) conv_flag = true;
  else conv_flag = false;
}



void PNonlinear_Tet4_NS_3D_Solver::GenAlpha_VMS_solve(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &sol_base,
    const PDNSolution * const &pre_dot_sol,
    const PDNSolution * const &pre_sol,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Inflow_NodalBC * const &infnbc_part,
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
  PetscLogEvent lin_solve_event;
  PetscClassId classid_assembly;
  PetscClassIdRegister("mat_vec_assembly", &classid_assembly);
  PetscLogEventRegister("assembly mat 0", classid_assembly, &mat_assem_0_event);
  PetscLogEventRegister("assembly mat 1", classid_assembly, &mat_assem_1_event);
  PetscLogEventRegister("assembly vec 0", classid_assembly, &vec_assem_0_event);
  PetscLogEventRegister("assembly vec 1", classid_assembly, &vec_assem_1_event);
  PetscLogEventRegister("lin_solve", classid_assembly, &lin_solve_event);
#endif

  nl_counter = 0;
  double residual_norm = 0.0, initial_norm = 0.0, relative_error = 0.0;

  // Gen-alpha parameters
  const double gamma   = tmga_ptr->get_gamma();
  const double alpha_m = tmga_ptr->get_alpha_m();
  const double alpha_f = tmga_ptr->get_alpha_f();

  // predictor
  sol -> ScaleValue(0.0);
  dot_sol -> ScaleValue(0.0);

  sol     -> PlusAX(*pre_sol, 1.0);
  dot_sol -> PlusAX(*pre_dot_sol, (gamma-1.0)/gamma);


  // define alpha_f dot_sol & alpha_m sol
  sol_alpha->Copy(pre_sol);

  // -------------------------------------
  // Update the infow boundary value here
  rescale_inflow_value(curr_time + dt,           infnbc_part, sol_base, sol); 
  rescale_inflow_value(curr_time + alpha_f * dt, infnbc_part, sol_base, sol_alpha); 
  // -------------------------------------

  dot_sol_alpha->Copy(pre_dot_sol);
  dot_sol_alpha->ScaleValue(1.0 - alpha_m);
  dot_sol_alpha->PlusAX(dot_sol, alpha_m);

  if( new_tangent_flag )
  {
    gassem_ptr->Clear_KG();

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(mat_assem_0_event, 0,0,0,0);
#endif

    gassem_ptr->Assem_tangent_residual( dot_sol_alpha, sol_alpha,
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
    
    gassem_ptr->Assem_residual( dot_sol_alpha, sol_alpha,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part );
  
#ifdef PETSC_USE_LOG
    PetscLogEventEnd(vec_assem_0_event,0,0,0,0);
#endif
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  PetscPrintf(PETSC_COMM_WORLD, "  Init res 2-norm: %e \n", initial_norm);

  do
  {
#ifdef PETSC_USE_LOG
    PetscLogEventBegin(lin_solve_event, 0,0,0,0);
#endif
  
    lsolver_ptr->Solve( gassem_ptr->G, velo_step );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(lin_solve_event,0,0,0,0);
#endif

    nl_counter += 1;

    // corrector
    sol->PlusAiX(*velo_step, -1.0, -1.0*gamma*dt, 1, 3);
    dot_sol->PlusAX( velo_step, -1.0 );

    sol_alpha->PlusAiX(*velo_step, -1.0, -1.0 * alpha_f * gamma * dt, 1, 3);
    dot_sol_alpha->PlusAX(velo_step, -1.0 * alpha_m);

    if( nl_counter >= nrenew_freq )
    {
      gassem_ptr->Clear_KG();
      
#ifdef PETSC_USE_LOG
      PetscLogEventBegin(mat_assem_1_event, 0,0,0,0);
#endif
      
      gassem_ptr->Assem_tangent_residual( dot_sol_alpha, sol_alpha,
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
      
      gassem_ptr->Assem_residual( dot_sol_alpha, sol_alpha,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part );
    
#ifdef PETSC_USE_LOG
      PetscLogEventEnd(vec_assem_1_event,0,0,0,0);
#endif
    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);

    PetscPrintf(PETSC_COMM_WORLD, "  --- res norm: %e \n", residual_norm);

    // Avoid undefined division by zero
    if(initial_norm != 0.0)
      relative_error = residual_norm / initial_norm;
    else
      relative_error = residual_norm;

    // Break if diverged
    if( relative_error >= nd_tol )
    {
      PetscPrintf(PETSC_COMM_WORLD,
          "Warning: nonlinear solver is diverging with error %e \n",
          relative_error);
      break;
    }
  }while(nl_counter < nmaxits && relative_error > nr_tol &&
      residual_norm > na_tol );

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol) conv_flag = true;
  else conv_flag = false;
}



void PNonlinear_Tet4_NS_3D_Solver::rescale_inflow_value( const double &stime,
    const ALocal_Inflow_NodalBC * const &infbc,
    const PDNSolution * const &sol_base,
    PDNSolution * const &sol ) const
{
  // Get the number of nodes of the inflow bc surface belonging to the CPU
  const int numnode = infbc -> get_Num_LD();

  // ----------------------------------------------------------------
  // Define the fourier series as a function of time to scale the flow rate Q
  const double a0 = 90.72 ;
  const double a1 = 63.24 ;
  const double b1 = 136.2 ;
  const double a2 = -53.14;
  const double b2 = 70.93;
  const double w = 6.246;
  const double period = 1.0;

  const double tt = stime - floor( stime / period );
  
  //const double val = a0 + a1 * cos(tt* w) + b1 * sin(tt * w)
  //  + a2 * cos(tt * 2.0 * w) + b2 * sin(tt * 2.0 * w);
  const double val = a0 + a1 * sin(tt*w);
  // ----------------------------------------------------------------
 
  double base_vals[3];
  int base_idx[3];

  // Initialize random generator
  //srand ( (unsigned int)time(NULL) );
  //std::vector<double> pert;
  //MATH_T::gen_Gaussian(3*numnode, 0.0, 0.01, pert);

  for(int ii=0; ii<numnode; ++ii)
  {
    int node_index = infbc -> get_LDN(ii);

    base_idx[0] = node_index * 4 + 1;
    base_idx[1] = node_index * 4 + 2;
    base_idx[2] = node_index * 4 + 3;

    VecGetValues(sol_base->solution, 3, base_idx, base_vals);

    const double del_x = 0.0; //pert[3*ii];
    const double del_y = 0.0; //pert[3*ii+1];
    const double del_z = 1.0; //1.0 + pert[3*ii+2];

    VecSetValue(sol->solution, node_index*4+1, base_vals[2] * val * del_x, INSERT_VALUES);
    VecSetValue(sol->solution, node_index*4+2, base_vals[2] * val * del_y, INSERT_VALUES);
    VecSetValue(sol->solution, node_index*4+3, base_vals[2] * val * del_z, INSERT_VALUES);
  }

  VecAssemblyBegin(sol->solution); VecAssemblyEnd(sol->solution);
  sol->GhostUpdate();
}

// EOF
