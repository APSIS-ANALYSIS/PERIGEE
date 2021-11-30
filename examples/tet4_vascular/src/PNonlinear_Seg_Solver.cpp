#include "PNonlinear_Seg_Solver.hpp"

PNonlinear_Seg_Solver::PNonlinear_Seg_Solver(
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol,
    const int &input_max_iteration, const int &input_renew_freq )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq)
{
  // Generate the solution vector for P-V system as the step updater
  dot_P_V_step = new PDNSolution_P_V_Mixed_3D( anode_ptr, feanode_ptr, 0, false );

  // Generate the solution vector for the mesh displacement
  mesh_disp = new PDNSolution_U_Mixed_3D( anode_ptr, feanode_ptr, 0, false );
}


PNonlinear_Seg_Solver::~PNonlinear_Seg_Solver()
{
  delete dot_P_V_step; dot_P_V_step = nullptr;
  delete mesh_disp;    mesh_disp    = nullptr;
}


void PNonlinear_Seg_Solver::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("relative tolerance: %e \n", nr_tol);
  SYS_T::commPrint("absolute tolerance: %e \n", na_tol);
  SYS_T::commPrint("divergence tolerance: %e \n", nd_tol);
  SYS_T::commPrint("maximum iteration: %d \n", nmaxits);
  SYS_T::commPrint("tangent matrix renew frequency: %d \n", nrenew_freq);
  SYS_T::commPrint("----------------------------------------------------------- \n");
}


void PNonlinear_Seg_Solver::GenAlpha_Seg_solve_ALE_NS(
    const bool &new_tangent_flag,
    const bool &is_ale_flag,
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
    const ALocal_NodalBC * const &nbc_mesh_part,
    const ALocal_EBC * const &ebc_part,
    const ALocal_EBC * const &ebc_mesh_part,
    const IGenBC * const &gbc,
    const Matrix_PETSc * const &bc_mat,
    const Matrix_PETSc * const &bc_mesh_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_solid_ptr,
    IPLocAssem * const &lassem_mesh_ptr,
    IPGAssem * const &gassem_ptr,
    IPGAssem * const &gassem_mesh_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PLinear_Solver_PETSc * const &lsolver_mesh_ptr,
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

  // ------------------------------------------------- 
  // Update the inflow boundary values
  rescale_inflow_value(curr_time+dt,         infnbc_part, flr_ptr, sol_base, sol);
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
        curr_time, dt, alelem_ptr, lassem_solid_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part, gbc );
   
#ifdef PETSC_USE_LOG
    PetscLogEventEnd(mat_assem_0_event,0,0,0,0);
#endif

    SYS_T::commPrint("  --- M updated");
    lsolver_ptr->SetOperator(gassem_ptr->K);
  }
  else
  {
    gassem_ptr->Clear_G();

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(vec_assem_0_event, 0,0,0,0);
#endif

    gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
        curr_time, dt, alelem_ptr, lassem_solid_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part, gbc );

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

    lsolver_ptr->Solve( gassem_ptr->G, dot_P_V_step );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(lin_solve_event,0,0,0,0);
#endif

    bc_mat->MatMultSol( dot_P_V_step );

    nl_counter += 1;

    dot_step.ScaleValue(0.0);

    // In fluid we do not need to update the mesh disp by velocity
    //SEG_SOL_T::PlusAiPV( val_1 * (-1.0), -1.0, -1.0, dot_P_V_step, &dot_step );
    SEG_SOL_T::PlusAiPV( 0.0, -1.0, -1.0, dot_P_V_step, &dot_step );

    dot_sol->PlusAX( dot_step, 1.0 );
    sol->PlusAX( dot_step, gamma * dt );

    dot_sol_alpha.PlusAX( dot_step, alpha_m );
    sol_alpha.PlusAX( dot_step, alpha_f * gamma * dt );

    // --------------------------------------------------------------
    // --- if ale is on, solve the mesh equation and update the fluid mesh
    if(is_ale_flag)
    {
      gassem_mesh_ptr->Clear_G();
      gassem_mesh_ptr->Assem_residual( sol, sol,
          curr_time, dt, alelem_ptr, lassem_mesh_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_mesh_part, ebc_mesh_part );
      lsolver_mesh_ptr->Solve( gassem_mesh_ptr->G, mesh_disp ); 
      bc_mesh_mat -> MatMultSol( mesh_disp );
      SEG_SOL_T::UpdateU( -1.0, mesh_disp, sol );

      SEG_SOL_T::UpdateV( dt, gamma, pre_dot_sol, pre_sol, sol, dot_sol);
    }
    // --- Finish mesh update
    // --------------------------------------------------------------

    // Assembly ALE-NS residual (& tangent) 
    if( nl_counter % nrenew_freq == 0 || nl_counter >= 4 )
    {
      gassem_ptr->Clear_KG();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(mat_assem_1_event, 0,0,0,0);
#endif

      gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
          curr_time, dt, alelem_ptr, lassem_solid_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part, gbc );

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
          curr_time, dt, alelem_ptr, lassem_solid_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part, gbc );

#ifdef PETSC_USE_LOG
      PetscLogEventEnd(vec_assem_1_event,0,0,0,0);
#endif
    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    SYS_T::commPrint("  --- nl_res: %e \n", residual_norm);

    relative_error = residual_norm / initial_norm;

    if( relative_error >= nd_tol )
    {
      SYS_T::commPrint("Warning: nonlinear solver is diverging with error %e. \n", relative_error);
      break;
    }
  }while(nl_counter<nmaxits && relative_error > nr_tol && residual_norm > na_tol);

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol) conv_flag = true;
  else conv_flag = false;
}


void PNonlinear_Seg_Solver::GenAlpha_Seg_solve_FSI(
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
    const ALocal_NodalBC * const &nbc_mesh_part,
    const ALocal_EBC * const &ebc_part,
    const ALocal_EBC * const &ebc_mesh_part,
    const IGenBC * const &gbc,
    const Matrix_PETSc * const &bc_mat,
    const Matrix_PETSc * const &bc_mesh_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_fluid_ptr,
    IPLocAssem * const &lassem_solid_ptr,
    IPLocAssem * const &lassem_mesh_ptr,
    IPGAssem * const &gassem_ptr,
    IPGAssem * const &gassem_mesh_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PLinear_Solver_PETSc * const &lsolver_mesh_ptr,
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

  // Generate the U-P-V step update: dot_step
  PDNSolution dot_step(*dot_sol);

  // Define the dol_sol at alpha_m: dot_sol_alpha
  PDNSolution dot_sol_alpha(*pre_dot_sol);
  dot_sol_alpha.ScaleValue( 1.0 - alpha_m );
  dot_sol_alpha.PlusAX(*dot_sol, alpha_m);

  // Define the sol at alpha_f: sol_alpha
  PDNSolution sol_alpha(*pre_sol);
  sol_alpha.ScaleValue( 1.0 - alpha_f );
  sol_alpha.PlusAX( *sol, alpha_f );

  // Get the dot_step from the dot_P_V_step only for the solid sub-domain
  dot_step.ScaleValue(0.0);
  SEG_SOL_T::PlusAiUPV( val_2, 0.0, 0.0, anode_ptr, pre_dot_sol, &dot_step );
  SEG_SOL_T::PlusAiVPV( val_3, 0.0, 0.0, anode_ptr, pre_sol, &dot_step );

  // Now use dot_step to udpate solutions
  dot_sol->PlusAX( dot_step, 1.0 );
  sol->PlusAX( dot_step, gamma * dt );

  dot_sol_alpha.PlusAX( dot_step, alpha_m );
  sol_alpha.PlusAX( dot_step, alpha_f * gamma * dt );

  // ------------------------------------------------- 
  // Update the inflow boundary values
  rescale_inflow_value(curr_time+dt,         infnbc_part, flr_ptr, sol_base, sol);
  rescale_inflow_value(curr_time+alpha_f*dt, infnbc_part, flr_ptr, sol_base, &sol_alpha);
  // ------------------------------------------------- 

  // If new_tangent_flag == TRUE, update the tangent matrix;
  // otherwise, use the matrix from the previous time step
  if( new_tangent_flag )
  {
    gassem_ptr->Clear_KG();
    gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
        curr_time, dt, alelem_ptr, lassem_fluid_ptr, lassem_solid_ptr,
        elementv, elements, quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part, gbc );
    SYS_T::commPrint("  --- M updated");
    lsolver_ptr->SetOperator(gassem_ptr->K);
  }
  else
  {
    gassem_ptr->Clear_G();
    gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
        curr_time, dt, alelem_ptr, lassem_fluid_ptr, lassem_solid_ptr, 
        elementv, elements, quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part, gbc );
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  SYS_T::commPrint("  Init res 2-norm: %e \n", initial_norm);

  // Now do consistent Newton-Raphson iteration
  do
  {
    lsolver_ptr->Solve( gassem_ptr->G, dot_P_V_step );

    bc_mat->MatMultSol( dot_P_V_step );

    nl_counter += 1;

    dot_step.ScaleValue(0.0);

    // In fluid we do not need to update the mesh disp by velocity
    // the val_1 * -1.0 is only applied for the solid sub-domain.
    SEG_SOL_T::PlusAiPV( 0.0, -1.0, -1.0, dot_P_V_step, &dot_step );
    SEG_SOL_T::PlusAiPV( val_1 * (-1.0), 0.0, 0.0, anode_ptr, 
        dot_P_V_step, &dot_step );

    dot_sol->PlusAX( dot_step, 1.0 );
    sol->PlusAX( dot_step, gamma * dt );

    dot_sol_alpha.PlusAX( dot_step, alpha_m );
    sol_alpha.PlusAX( dot_step, alpha_f * gamma * dt );

    // --------------------------------------------------------------
    // --- (elastic) mesh equation solve
    gassem_mesh_ptr->Clear_KG();
    gassem_mesh_ptr->Assem_tangent_residual( pre_sol, sol,
        curr_time, dt, alelem_ptr, lassem_mesh_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_mesh_part, ebc_mesh_part );
    lsolver_mesh_ptr->Solve( gassem_mesh_ptr->K, gassem_mesh_ptr->G, mesh_disp ); 
    bc_mesh_mat -> MatMultSol( mesh_disp );

    // Use mesh disp to update the disp slots in the fluid domain 
    SEG_SOL_T::UpdateU( -1.0, mesh_disp, sol );
    SEG_SOL_T::UpdateU( -1.0 * alpha_f, mesh_disp, &sol_alpha );
    // Update the first three slots of the dot_sol vector in the fluid domain
    SEG_SOL_T::UpdateV( dt, gamma, pre_dot_sol, pre_sol, sol, dot_sol);
    SEG_SOL_T::UpdateV( dt, gamma/alpha_m, pre_dot_sol, pre_sol, sol, &dot_sol_alpha);
    // --- Finish mesh update
    // --------------------------------------------------------------

    // Assembly residual (& tangent) 
    if( nl_counter % nrenew_freq == 0 || nl_counter >= 4 )
    {
      gassem_ptr->Clear_KG();
      gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
          curr_time, dt, alelem_ptr, lassem_fluid_ptr, lassem_solid_ptr, 
          elementv, elements, quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part, gbc );
      SYS_T::commPrint("  - M");
      lsolver_ptr->SetOperator(gassem_ptr->K);
    }
    else
    {
      gassem_ptr->Clear_G();
      gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
          curr_time, dt, alelem_ptr, lassem_fluid_ptr, lassem_solid_ptr, 
          elementv, elements, quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part, gbc );
    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    SYS_T::commPrint("  --- nl_res: %e \n", residual_norm);

    relative_error = residual_norm / initial_norm;

    if( relative_error >= nd_tol )
    {
      SYS_T::commPrint("Warning: nonlinear solver is diverging with error %e. \n", relative_error);
      break;
    }
  }while(nl_counter<nmaxits && relative_error > nr_tol && residual_norm > na_tol);

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol) conv_flag = true;
  else conv_flag = false;
}


void PNonlinear_Seg_Solver::rescale_inflow_value( const double &stime,
    const ALocal_Inflow_NodalBC * const &infbc,
    const ICVFlowRate * const &flrate,
    const PDNSolution * const &sol_base,
    PDNSolution * const &sol ) const
{
  const int num_nbc = infbc -> get_num_nbc();

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    const int numnode = infbc -> get_Num_LD( nbc_id );

    const double val = flrate -> get_flow_rate( nbc_id, stime );

    for(int ii=0; ii<numnode; ++ii)
    {
      const int node_index = infbc -> get_LDN( nbc_id, ii );

      const int base_idx[3] = { node_index * 7 + 4, node_index * 7 + 5, node_index * 7 + 6 };

      double base_vals[3];

      VecGetValues(sol_base->solution, 3, base_idx, base_vals);

      VecSetValue(sol->solution, node_index*7+4, base_vals[0] * val, INSERT_VALUES);
      VecSetValue(sol->solution, node_index*7+5, base_vals[1] * val, INSERT_VALUES);
      VecSetValue(sol->solution, node_index*7+6, base_vals[2] * val, INSERT_VALUES);
    }
  }
  
  sol -> Assembly_GhostUpdate();
}

// EOF
