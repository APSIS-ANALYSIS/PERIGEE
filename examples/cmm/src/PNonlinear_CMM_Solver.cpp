#include "PNonlinear_CMM_Solver.hpp"

PNonlinear_CMM_Solver::PNonlinear_CMM_Solver(
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol,
    const int &input_max_iteration, 
    const int &input_renew_freq,
    const int &input_renew_threshold )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq),
  nrenew_threshold(input_renew_threshold)
{
}

PNonlinear_CMM_Solver::~PNonlinear_CMM_Solver()
{
}

void PNonlinear_CMM_Solver::print_info() const
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

void PNonlinear_CMM_Solver::GenAlpha_Solve_CMM(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &sol_base,
    const PDNSolution * const &pre_dot_sol,
    const PDNSolution * const &pre_sol,
    const PDNSolution * const &pre_dot_sol_wall_disp,
    const PDNSolution * const &pre_sol_wall_disp,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ICVFlowRate * const flr_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_InflowBC * const &infnbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_part,
    const ALocal_EBC * const &ebc_wall_part,
    const IGenBC * const &gbc,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    FEAElement * const &elementw,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PDNSolution * const &dot_sol,
    PDNSolution * const &sol,
    PDNSolution * const &dot_sol_wall_disp,
    PDNSolution * const &sol_wall_disp,
    int &nl_counter ) const
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
  // fluid solution
  sol->Copy(*pre_sol);
  dot_sol->Copy(*pre_dot_sol);
  dot_sol->ScaleValue( (gamma-1.0)/gamma );

  // wall solution
  sol_wall_disp->Copy(*pre_sol_wall_disp);
  dot_sol_wall_disp->Copy(*pre_dot_sol_wall_disp);
  dot_sol_wall_disp->ScaleValue( (gamma-1.0)/gamma );

  // Define the dol_sol at alpha_m: dot_sol_alpha, dot_wall_disp_alpha
  PDNSolution dot_sol_alpha(*pre_dot_sol);
  dot_sol_alpha.ScaleValue( 1.0 - alpha_m );
  dot_sol_alpha.PlusAX(*dot_sol, alpha_m);

  PDNSolution dot_wall_disp_alpha(*pre_dot_sol_wall_disp);
  dot_wall_disp_alpha.ScaleValue( 1.0 - alpha_m );
  dot_wall_disp_alpha.PlusAX(*dot_sol_wall_disp, alpha_m);

  // Define the sol at alpha_f: sol_alpha, wall_disp_alpha
  PDNSolution sol_alpha(*pre_sol);
  sol_alpha.ScaleValue( 1.0 - alpha_f );
  sol_alpha.PlusAX( *sol, alpha_f );

  PDNSolution wall_disp_alpha(*pre_sol_wall_disp);
  wall_disp_alpha.ScaleValue( 1.0 - alpha_f );
  wall_disp_alpha.PlusAX( *sol_wall_disp, alpha_f );

  // Compute kinematic residual = dot_wall_disp_alpha - velo_alpha
  PDNSolution G_kinematic(dot_wall_disp_alpha);
    
  update_wall(-1.0, &sol_alpha, ebc_wall_part, &G_kinematic);

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

    gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, &wall_disp_alpha,
        dot_sol, sol, curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        elementw, quad_v, quad_s, lien_ptr, feanode_ptr, 
        nbc_part, ringnbc_part, ebc_part, ebc_wall_part, gbc );
   
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

    gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha, &wall_disp_alpha,
        dot_sol, sol, curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        elementw, quad_v, quad_s, lien_ptr, feanode_ptr, 
        nbc_part, ringnbc_part, ebc_part, ebc_wall_part, gbc );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(vec_assem_0_event,0,0,0,0);
#endif
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  SYS_T::commPrint("  Init res 2-norm: %e \n", initial_norm);

  // Allocate the incremental solution vector
  PDNSolution * dot_step = new PDNSolution( sol );

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

    // Skew boundary conditions for in-plane motion of ring nodes:
    // Rotate ring node velo dofs back into the global Cartesian frame
    rotate_ringbc( ringnbc_part, dot_step );

    // Update dot_sol, dot_sol_wall_disp
    dot_sol->PlusAX( dot_step, -1.0 );
    update_wall( (-1.0) * alpha_f * gamma * dt / alpha_m,
        dot_step, ebc_wall_part, dot_sol_wall_disp ); 

    dot_sol_wall_disp->PlusAX(G_kinematic, (-1.0) / alpha_m );

    // Update sol, sol_wall_disp
    sol->PlusAX( dot_step, (-1.0) * gamma * dt );
    update_wall( (-1.0) * alpha_f * gamma * gamma * dt * dt / alpha_m,
        dot_step, ebc_wall_part, sol_wall_disp ); 
    sol_wall_disp->PlusAX(G_kinematic, (-1.0) * gamma * dt / alpha_m );

    // Update dol_sol at alpha_m: dot_sol_alpha, dot_wall_disp_alpha
    dot_sol_alpha.PlusAX( dot_step, (-1.0) * alpha_m );
    update_wall( (-1.0) * alpha_f * gamma * dt,
        dot_step, ebc_wall_part, &dot_wall_disp_alpha ); 
    dot_wall_disp_alpha.PlusAX(G_kinematic, -1.0 );

    // Update sol at alpha_f: sol_alpha, wall_disp_alpha
    sol_alpha.PlusAX( dot_step, (-1.0) * alpha_f * gamma * dt );
    update_wall( (-1.0) * alpha_f * alpha_f * gamma * gamma * dt * dt / alpha_m,
        dot_step, ebc_wall_part, &wall_disp_alpha ); 
    wall_disp_alpha.PlusAX(G_kinematic, (-1.0) * alpha_f * gamma * dt / alpha_m );

    // Update kinematic residual = dot_wall_disp_alpha - velo_alpha
    G_kinematic.Copy(dot_wall_disp_alpha);
    update_wall(-1.0, &sol_alpha, ebc_wall_part, &G_kinematic);

    nl_counter += 1;
    
    // Assembly residual (& tangent if condition satisfied) 
    if( nl_counter % nrenew_freq == 0 || nl_counter >= nrenew_threshold )
    {
      gassem_ptr->Clear_KG();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(mat_assem_1_event, 0,0,0,0);
#endif

      gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, &wall_disp_alpha,
          dot_sol, sol, curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          elementw, quad_v, quad_s, lien_ptr, feanode_ptr, 
          nbc_part, ringnbc_part, ebc_part, ebc_wall_part, gbc );

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

      gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha, &wall_disp_alpha,
          dot_sol, sol, curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          elementw, quad_v, quad_s, lien_ptr, feanode_ptr, 
          nbc_part, ringnbc_part, ebc_part, ebc_wall_part, gbc );

#ifdef PETSC_USE_LOG
      PetscLogEventEnd(vec_assem_1_event,0,0,0,0);
#endif
    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    
    SYS_T::commPrint("  --- nl_res: %e \n", residual_norm);

    // Print the residual norm of the kienmatic equation
    SYS_T::commPrint("  --- kinematic_res: %e \n", G_kinematic.Norm_2());

    relative_error = residual_norm / initial_norm;

    SYS_T::print_fatal_if( relative_error >= nd_tol, "Error: nonlinear solver is diverging with error %e. Job killed.\n", relative_error);

  }while(nl_counter<nmaxits && relative_error > nr_tol && residual_norm > na_tol);

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  // free the incremental solution vector
  delete dot_step; dot_step = nullptr;

  // Debugging: check ring BC constraints
  // compute_ringbc_constraints(sol, sol_wall_disp, ringnbc_part);
}


void PNonlinear_CMM_Solver::GenAlpha_Solve_Prestress(
    const bool &new_tangent_flag,
    const double &prestress_tol,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &sol_base,
    const PDNSolution * const &pre_dot_sol,
    const PDNSolution * const &pre_sol,
    const PDNSolution * const &pre_dot_sol_wall_disp,
    const PDNSolution * const &pre_sol_wall_disp,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ICVFlowRate * const flr_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_InflowBC * const &infnbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_part,
    ALocal_EBC * const &ebc_wall_part,
    const IGenBC * const &gbc,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    FEAElement * const &elementw,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PDNSolution * const &dot_sol,
    PDNSolution * const &sol,
    PDNSolution * const &dot_sol_wall_disp,
    PDNSolution * const &sol_wall_disp,
    bool &prestress_conv_flag, int &nl_counter ) const
{
  // Initialize the counter and error
  nl_counter = 0;
  double residual_norm = 0.0, initial_norm = 0.0, relative_error = 0.0;

  // Gen-alpha parameters
  const double gamma   = tmga_ptr->get_gamma();
  const double alpha_m = tmga_ptr->get_alpha_m();
  const double alpha_f = tmga_ptr->get_alpha_f();

  // Same-Y predictor
  // fluid solution
  sol->Copy(*pre_sol);
  dot_sol->Copy(*pre_dot_sol);
  dot_sol->ScaleValue( (gamma-1.0)/gamma );

  // wall solution
  sol_wall_disp->Copy(*pre_sol_wall_disp);
  dot_sol_wall_disp->Copy(*pre_dot_sol_wall_disp);
  dot_sol_wall_disp->ScaleValue( (gamma-1.0)/gamma );

  // Define the dol_sol at alpha_m: dot_sol_alpha, dot_wall_disp_alpha
  PDNSolution dot_sol_alpha(*pre_dot_sol);
  dot_sol_alpha.ScaleValue( 1.0 - alpha_m );
  dot_sol_alpha.PlusAX(*dot_sol, alpha_m);

  PDNSolution dot_wall_disp_alpha(*pre_dot_sol_wall_disp);
  dot_wall_disp_alpha.ScaleValue( 1.0 - alpha_m );
  dot_wall_disp_alpha.PlusAX(*dot_sol_wall_disp, alpha_m);

  // Define the sol at alpha_f: sol_alpha, wall_disp_alpha
  PDNSolution sol_alpha(*pre_sol);
  sol_alpha.ScaleValue( 1.0 - alpha_f );
  sol_alpha.PlusAX( *sol, alpha_f );

  PDNSolution wall_disp_alpha(*pre_sol_wall_disp);
  wall_disp_alpha.ScaleValue( 1.0 - alpha_f );
  wall_disp_alpha.PlusAX( *sol_wall_disp, alpha_f );

  // Compute kinematic residual = dot_wall_disp_alpha - velo_alpha
  PDNSolution G_kinematic(dot_wall_disp_alpha);
    
  update_wall(-1.0, &sol_alpha, ebc_wall_part, &G_kinematic);

  gassem_ptr->Clear_KG();

  gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, &wall_disp_alpha,
      dot_sol, sol, curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
      elementw, quad_v, quad_s, lien_ptr, feanode_ptr, 
      nbc_part, ringnbc_part, ebc_part, ebc_wall_part, gbc );

  SYS_T::commPrint("  --- M updated");

  // SetOperator will pass the tangent matrix to the linear solver and the
  // linear solver will generate the preconditioner based on the new matrix.
  lsolver_ptr->SetOperator( gassem_ptr->K );

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  SYS_T::commPrint("  Init res 2-norm: %e \n", initial_norm);

  // Allocate the incremental solution vector
  PDNSolution * dot_step = new PDNSolution( sol );

  // Now do consistent Newton-Raphson iteration
  do
  {
    // solve the equation K dot_step = G
    lsolver_ptr->Solve( gassem_ptr->G, dot_step );

    bc_mat->MatMultSol( dot_step );

    // Skew boundary conditions for in-plane motion of ring nodes:
    // Rotate ring node velo dofs back into the global Cartesian frame
    rotate_ringbc( ringnbc_part, dot_step );

    // Update dot_sol, dot_sol_wall_disp
    dot_sol->PlusAX( dot_step, -1.0 );
    update_wall( (-1.0) * alpha_f * gamma * dt / alpha_m,
        dot_step, ebc_wall_part, dot_sol_wall_disp ); 

    dot_sol_wall_disp->PlusAX(G_kinematic, (-1.0) / alpha_m );

    // Update sol, sol_wall_disp
    sol->PlusAX( dot_step, (-1.0) * gamma * dt );
    update_wall( (-1.0) * alpha_f * gamma * gamma * dt * dt / alpha_m,
        dot_step, ebc_wall_part, sol_wall_disp ); 
    sol_wall_disp->PlusAX(G_kinematic, (-1.0) * gamma * dt / alpha_m );

    // Update dol_sol at alpha_m: dot_sol_alpha, dot_wall_disp_alpha
    dot_sol_alpha.PlusAX( dot_step, (-1.0) * alpha_m );
    update_wall( (-1.0) * alpha_f * gamma * dt,
        dot_step, ebc_wall_part, &dot_wall_disp_alpha ); 
    dot_wall_disp_alpha.PlusAX(G_kinematic, -1.0 );

    // Update sol at alpha_f: sol_alpha, wall_disp_alpha
    sol_alpha.PlusAX( dot_step, (-1.0) * alpha_f * gamma * dt );
    update_wall( (-1.0) * alpha_f * alpha_f * gamma * gamma * dt * dt / alpha_m,
        dot_step, ebc_wall_part, &wall_disp_alpha ); 
    wall_disp_alpha.PlusAX(G_kinematic, (-1.0) * alpha_f * gamma * dt / alpha_m );

    // Update kinematic residual = dot_wall_disp_alpha - velo_alpha
    G_kinematic.Copy(dot_wall_disp_alpha);
    update_wall(-1.0, &sol_alpha, ebc_wall_part, &G_kinematic);

    nl_counter += 1;

    // Assembly residual & tangent
    gassem_ptr->Clear_KG();

    gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, &wall_disp_alpha,
        dot_sol, sol, curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        elementw, quad_v, quad_s, lien_ptr, feanode_ptr, 
        nbc_part, ringnbc_part, ebc_part, ebc_wall_part, gbc );

    SYS_T::commPrint("  --- M updated");
    lsolver_ptr->SetOperator(gassem_ptr->K);

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);

    SYS_T::commPrint("  --- nl_res: %e \n", residual_norm);

    // Print the residual norm of the kienmatic equation
    SYS_T::commPrint("  --- kinematic_res: %e \n", G_kinematic.Norm_2());

    relative_error = residual_norm / initial_norm;

    SYS_T::print_fatal_if( relative_error >= nd_tol, "Error: nonlinear solver is diverging with error %e. Job killed.\n", relative_error);

  }while(nl_counter<nmaxits && relative_error > nr_tol && residual_norm > na_tol);

  // free the incremental solution vector
  delete dot_step; dot_step = nullptr;

  // --------------------------------------------------------------------------
  // Update the prestress values
  gassem_ptr->Update_Wall_Prestress(sol_wall_disp, lassem_ptr, elementw, quad_s, ebc_wall_part);

  SYS_T::commPrint("  --- wall_disp_norm: %e \n", sol_wall_disp->Norm_2());

  if( sol_wall_disp->Norm_2() <= prestress_tol ) prestress_conv_flag = true;
  // --------------------------------------------------------------------------

  Print_convergence_info(nl_counter, relative_error, residual_norm);
}

void PNonlinear_CMM_Solver::rescale_inflow_value( const double &stime,
    const ALocal_InflowBC * const &infbc,
    const ICVFlowRate * const &flrate,
    const PDNSolution * const &sol_base,
    PDNSolution * const &sol ) const
{
  const int num_nbc = infbc -> get_num_nbc();

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    const int numnode = infbc -> get_Num_LD( nbc_id );
    const double factor = flrate -> get_flow_rate( nbc_id, stime );

    for(int ii=0; ii<numnode; ++ii)
    {
      const int node_index = infbc -> get_LDN( nbc_id, ii );

      const int base_idx[3] = { node_index*4+1, node_index*4+2, node_index*4+3 };

      double base_vals[3];
      
      VecGetValues(sol_base->solution, 3, base_idx, base_vals);

      const double vals[3] = { base_vals[0] * factor, base_vals[1] * factor,
          base_vals[2] * factor };

      VecSetValues(sol->solution, 3, base_idx, vals, INSERT_VALUES);
    }
  }

  sol->Assembly_GhostUpdate();
}


void PNonlinear_CMM_Solver::update_wall( const double &val,
    const PDNSolution * const &dot_step,
    const ALocal_EBC * const &ebc_wall_part,
    PDNSolution * const &wall_data ) const
{
  // Verify that the dof of dot_step is 4
  SYS_T::print_fatal_if(dot_step->get_dof_num() != 4,
      "Error in PNonlinear_CMM_Solver::update_dot_wall_disp: incorrect dimension of dot_step. \n");

  // Verify that the dof of wall_data is 3
  SYS_T::print_fatal_if(wall_data->get_dof_num() != 3,
      "Error in PNonlinear_CMM_Solver::update_dot_wall_disp: incorrect dimension of wall_data. \n");

  // Verify consistency in the number of local nodes
  SYS_T::print_fatal_if( !is_layout_equal(*dot_step, *wall_data), "Error in PNonlinear_CMM_Solver::update_dot_wall_disp: solution vector layout mismatch between dot_step and wall_data. \n");

  Vec ldotstep, lwalldata;
  double * array_dotstep, * array_walldata;

  VecGhostGetLocalForm(dot_step->solution, &ldotstep);
  VecGhostGetLocalForm(wall_data->solution, &lwalldata);

  VecGetArray(ldotstep, &array_dotstep);
  VecGetArray(lwalldata, &array_walldata);

  const int num_snode = ebc_wall_part -> get_num_local_node_on_sur();

  for(int ii=0; ii<num_snode; ++ii)
  {
    const int pos = ebc_wall_part -> get_local_node_on_sur_pos(ii);

    array_walldata[pos*3]   += val * array_dotstep[pos*4+1];
    array_walldata[pos*3+1] += val * array_dotstep[pos*4+2];
    array_walldata[pos*3+2] += val * array_dotstep[pos*4+3];
  }

  // Deallocation of the local copy
  VecRestoreArray(ldotstep, &array_dotstep);
  VecRestoreArray(lwalldata, &array_walldata);
  VecGhostRestoreLocalForm(dot_step->solution, &ldotstep);
  VecGhostRestoreLocalForm(wall_data->solution, &lwalldata);

  // Update ghost values
  wall_data->GhostUpdate();
}


void PNonlinear_CMM_Solver::rotate_ringbc(
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    PDNSolution * const &dot_step) const
{
  const int ringbc_type = ringnbc_part -> get_ringbc_type();

  // Clamped rings
  if( ringbc_type == 0 ) {}

  // Skew boundary conditions for in-plane motion of ring nodes
  else if( ringbc_type == 1)
  {
    double vals[3], rot_vals[3];

    const int num_ringnode = ringnbc_part -> get_Num_LD();

    for(int ii = 0; ii < num_ringnode; ++ii)
    {
      const int dnode = ringnbc_part -> get_LDN( ii );

      const int idx[3] = { dnode*4 + 1, dnode*4 + 2, dnode*4 + 3 };

      VecGetValues(dot_step->solution, 3, idx, vals);

      Matrix_3x3 Q = ringnbc_part->get_rotation_matrix( ii );
      Q.transpose(); // Skew-to-global transformation matrix 

      // rot_vals = Q * vals
      Q.VecMult( vals[0], vals[1], vals[2], rot_vals[0], rot_vals[1], rot_vals[2] );
     
      VecSetValues(dot_step->solution, 3,  idx, rot_vals, INSERT_VALUES);
    }

    dot_step->Assembly_GhostUpdate();
  }
  else
    SYS_T::print_fatal("Error: this ringbc_type is not supported in PNonlinear_CMM_Solver.\n");
}


void PNonlinear_CMM_Solver::compute_ringbc_constraints(
    const PDNSolution * const &sol,
    const PDNSolution * const &sol_wall_disp,
    const ALocal_Ring_NodalBC * const &ringnbc_part ) const
{
  const int num_ringnode = ringnbc_part -> get_Num_LD();

  // Check whether the ring BC constraints are properly satisfied
  // by printing their evaluations 
  if( SYS_T::get_MPI_rank() == 0 && num_ringnode > 0 )
  {
    for(int ii = 0; ii < num_ringnode; ++ii)
    {
      double velo_val[3], disp_val[3];

      const int dnode  = ringnbc_part -> get_LDN( ii );

      const double outvec[3] = {ringnbc_part -> get_outvec(ii, 0), ringnbc_part -> get_outvec(ii, 1), ringnbc_part -> get_outvec(ii, 2)};

      const int velo_idx[3] = {dnode*4 + 1, dnode*4 + 2, dnode*4 + 3};
      const int disp_idx[3] = {dnode*3 + 0, dnode*3 + 1, dnode*3 + 2};

      VecGetValues(sol->solution, 3, velo_idx, velo_val);
      VecGetValues(sol_wall_disp->solution, 3, disp_idx, disp_val);

      const double v_dot_n = velo_val[0] * outvec[0] + velo_val[1] * outvec[1] + velo_val[2] * outvec[2];
      const double u_dot_n = disp_val[0] * outvec[0] + disp_val[1] * outvec[1] + disp_val[2] * outvec[2];

      std::cout << std::scientific << std::setprecision(3) << "Ring node " << std::setw(8) << dnode << ": ";
      std::cout << "Velo=[" << std::setw(10) << velo_val[0] << ", " << std::setw(10) << velo_val[1] << ", " << std::setw(10) << velo_val[2] << "], "; 
      std::cout << "Disp=[" << std::setw(10) << disp_val[0] << ", " << std::setw(10) << disp_val[1] << ", " << std::setw(10) << disp_val[2] << "], ";  
      std::cout << "v_dot_n = " << std::setw(10) << v_dot_n << ", u_dot_n = " << std::setw(10) << u_dot_n << std::endl;
    }
  }
}

// EOF
