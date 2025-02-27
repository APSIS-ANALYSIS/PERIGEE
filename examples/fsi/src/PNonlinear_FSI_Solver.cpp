#include "PNonlinear_FSI_Solver.hpp"

PNonlinear_FSI_Solver::PNonlinear_FSI_Solver(
    std::unique_ptr<IPGAssem> in_gassem_mesh,  
    std::unique_ptr<PLinear_Solver_PETSc> in_lsolver,
    std::unique_ptr<PLinear_Solver_PETSc> in_lsolver_mesh,
    std::unique_ptr<Matrix_PETSc> in_bc_mat,
    std::unique_ptr<Matrix_PETSc> in_bc_mesh_mat,
    std::unique_ptr<TimeMethod_GenAlpha> in_tmga,
    std::unique_ptr<IFlowRate> in_flrate,
    std::unique_ptr<PDNSolution> in_sol_base,
    std::unique_ptr<APart_Node> in_pnode_v,
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol,
    const int &input_max_iteration, 
    const int &input_renew_freq,
    const int &input_renew_threshold )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq),
  nrenew_threshold(input_renew_threshold),
  gassem_mesh(std::move(in_gassem_mesh)),
  lsolver(std::move(in_lsolver)),
  lsolver_mesh(std::move(in_lsolver_mesh)),
  bc_mat(std::move(in_bc_mat)),
  bc_mesh_mat(std::move(in_bc_mesh_mat)),
  tmga(std::move(in_tmga)),
  flrate(std::move(in_flrate)),
  sol_base(std::move(in_sol_base)),
  pnode_v(std::move(in_pnode_v))
{}

PNonlinear_FSI_Solver::~PNonlinear_FSI_Solver()
{}

void PNonlinear_FSI_Solver::print_info() const
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

void PNonlinear_FSI_Solver::update_solid_kinematics( 
    const double &val, const Vec &input, 
    PDNSolution * const &output ) const
{
  const int nlocal = pnode_v -> get_nlocalnode_solid();
  
  Vec local_output;

  VecGhostGetLocalForm(output->solution, &local_output);

  double * array_input, * array_output;

  VecGetArray(input,        &array_input);
  VecGetArray(local_output, &array_output);

  for(int ii=0; ii<nlocal; ++ii)
  {
    const int ii3 = pnode_v->get_node_loc_solid(ii) * 3;

    array_output[ii3  ] += val * array_input[ii3  ];
    array_output[ii3+1] += val * array_input[ii3+1];
    array_output[ii3+2] += val * array_input[ii3+2];
  }

  VecRestoreArray(input,        &array_input);
  VecRestoreArray(local_output, &array_output);
  VecGhostRestoreLocalForm(output->solution, &local_output);

  output->GhostUpdate(); // update the ghost slots
}

void PNonlinear_FSI_Solver::rescale_inflow_value( const double &stime,
    const ALocal_InflowBC * const &infbc,
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

      const int base_idx[3] = { node_index * 3, node_index * 3 + 1, node_index * 3 + 2 };

      double base_vals[3];

      VecGetValues(sol_base->solution, 3, base_idx, base_vals);

      VecSetValue(sol->solution, node_index*3  , base_vals[0] * val, INSERT_VALUES);
      VecSetValue(sol->solution, node_index*3+1, base_vals[1] * val, INSERT_VALUES);
      VecSetValue(sol->solution, node_index*3+2, base_vals[2] * val, INSERT_VALUES);
    }
  }

  sol -> Assembly_GhostUpdate();
}

void PNonlinear_FSI_Solver::GenAlpha_Seg_solve_FSI(
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
    const ALocal_InflowBC * const &infnbc_part,
    const IGenBC * const &gbc,
    const Tissue_prestress * const &ps_ptr,
    IPGAssem * const &gassem_ptr,
    PDNSolution * const &dot_disp,
    PDNSolution * const &dot_velo,
    PDNSolution * const &dot_pres,
    PDNSolution * const &disp,
    PDNSolution * const &velo,
    PDNSolution * const &pres,
    bool &conv_flag, int &nl_counter ) const
{
#ifdef PETSC_USE_LOG
  PetscLogEvent assem_event_0, assem_event_1, assem_event_2;
  PetscLogEvent solve_mech_event, solve_mesh_event;
  PetscClassId classid;
  PetscClassIdRegister("user-log-info", &classid);
  PetscLogEventRegister("assembly_0", classid, &assem_event_0);
  PetscLogEventRegister("assembly_1", classid, &assem_event_1);
  PetscLogEventRegister("assembly_2", classid, &assem_event_2);
  PetscLogEventRegister("lin_solve_mech", classid, &solve_mech_event);
  PetscLogEventRegister("lin_solve_mesh", classid, &solve_mesh_event);
#endif

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
  std::unique_ptr<PDNSolution> Delta_dot_disp = SYS_T::make_unique<PDNSolution_V>(pnode_v.get(), 0, false, "delta_dot_disp");

  update_solid_kinematics( -1.0 / alpha_m, dot_disp_alpha->solution, Delta_dot_disp.get() );
  update_solid_kinematics(  1.0 / alpha_m, velo_alpha->solution, Delta_dot_disp.get() );

  // Now update displacement solutions for solid sub-domain
  dot_disp       -> PlusAX( Delta_dot_disp.get(), 1.0 );
  disp           -> PlusAX( Delta_dot_disp.get(), gamma * dt );
  dot_disp_alpha -> PlusAX( Delta_dot_disp.get(), alpha_m );
  disp_alpha     -> PlusAX( Delta_dot_disp.get(), alpha_f * gamma * dt );
  
  // Update inflow boundary values
  rescale_inflow_value( curr_time + dt,           infnbc_part, velo );
  rescale_inflow_value( curr_time + alpha_f * dt, infnbc_part, velo_alpha.get() );

#ifdef PETSC_USE_LOG
  PetscLogEventBegin(assem_event_0, 0,0,0,0);
#endif

  // If new_tangent_flag == TRUE, update the tangent matrix;
  // otherwise, use the matrix from the previous time step
  if( new_tangent_flag )
  {
    gassem_ptr -> Clear_KG();

    gassem_ptr->Assem_Tangent_Residual( curr_time, dt, 
        dot_disp_alpha.get(), dot_velo_alpha.get(), dot_pres_alpha.get(),
        disp_alpha.get(), velo_alpha.get(), pres_alpha.get(), 
        dot_velo, velo, disp, gbc, ps_ptr );

    SYS_T::commPrint("  --- M updated");
    lsolver->SetOperator(gassem_ptr->K);
  }
  else
  {
    gassem_ptr -> Clear_G();
    
    gassem_ptr->Assem_Residual( curr_time, dt, 
        dot_disp_alpha.get(), dot_velo_alpha.get(), dot_pres_alpha.get(),
        disp_alpha.get(), velo_alpha.get(), pres_alpha.get(), 
        dot_velo, velo, disp, gbc, ps_ptr );
  }

#ifdef PETSC_USE_LOG
  PetscLogEventEnd(assem_event_0, 0,0,0,0);
#endif

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  SYS_T::commPrint("  Init res 2-norm: %e \n", initial_norm);

  Vec sol_vp, sol_v, sol_p;
  VecDuplicate( gassem_ptr->G, &sol_vp );

  Vec sol_mesh;
  VecDuplicate( gassem_mesh->G, &sol_mesh );

  // Now we do consistent Newton-Raphson iteration
  do
  {
    // Check the residual dot_u_alpha - v_alpha
    Delta_dot_disp -> ScaleValue( 0.0 );
    update_solid_kinematics(  1.0, dot_disp_alpha->solution, Delta_dot_disp.get() );
    update_solid_kinematics( -1.0,     velo_alpha->solution, Delta_dot_disp.get() );
    const double solid_kinematics_residual = Delta_dot_disp -> Norm_2();
    // Finish calculating dot_u_alpha - v_alpha

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(solve_mech_event, 0,0,0,0);
#endif

    lsolver->Solve( gassem_ptr->G, sol_vp );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(solve_mech_event, 0,0,0,0);
#endif

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
    update_solid_kinematics( -1.0 * val_1 * alpha_f * gamma * dt, sol_v, disp_alpha.get());

    VecRestoreSubVector(sol_vp, is_v, &sol_v);
    VecRestoreSubVector(sol_vp, is_p, &sol_p);

    // Solve for mesh motion
    gassem_mesh -> Clear_G();

#ifdef PETSC_USE_LOG
  PetscLogEventBegin(assem_event_2, 0,0,0,0);
#endif

    gassem_mesh -> Assem_residual( pre_disp, disp, curr_time, dt );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(assem_event_2, 0,0,0,0);
#endif

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(solve_mesh_event, 0,0,0,0);
#endif
    
    lsolver_mesh -> Solve( gassem_mesh -> G, sol_mesh );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(solve_mesh_event,0,0,0,0);
#endif

    bc_mesh_mat -> MatMultSol( sol_mesh );

    // update the mesh displacement
    disp       -> PlusAX( sol_mesh, -1.0 );
    disp_alpha -> PlusAX( sol_mesh, -1.0 * alpha_f );

    // update the mesh velocity
    dot_disp -> Copy( pre_dot_disp );
    dot_disp -> ScaleValue( (gamma-1.0) / gamma );
    dot_disp -> PlusAX( pre_disp, -1.0 / (gamma * dt) );
    dot_disp -> PlusAX( disp,      1.0 / (gamma * dt) );

    dot_disp_alpha -> Copy( pre_dot_disp );
    dot_disp_alpha -> ScaleValue( 1.0 - alpha_m );
    dot_disp_alpha -> PlusAX( dot_disp, alpha_m );

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(assem_event_1, 0,0,0,0);
#endif

    // Assemble residual & tangent
    if( nl_counter % nrenew_freq == 0 || nl_counter >= nrenew_threshold )
    {
      gassem_ptr -> Clear_KG();

      gassem_ptr->Assem_Tangent_Residual( curr_time, dt,
          dot_disp_alpha.get(), dot_velo_alpha.get(), dot_pres_alpha.get(),
          disp_alpha.get(), velo_alpha.get(), pres_alpha.get(), 
          dot_velo, velo, disp, gbc, ps_ptr );

      SYS_T::commPrint("  --- M updated");
      lsolver->SetOperator(gassem_ptr->K);
    }
    else
    {
      gassem_ptr -> Clear_G();

      gassem_ptr->Assem_Residual( curr_time, dt,
          dot_disp_alpha.get(), dot_velo_alpha.get(), dot_pres_alpha.get(),
          disp_alpha.get(), velo_alpha.get(), pres_alpha.get(),
          dot_velo, velo, disp, gbc, ps_ptr );
    }

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(assem_event_1,0,0,0,0);
#endif

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    SYS_T::commPrint("  --- nl_res: %e \n", residual_norm);

    SYS_T::commPrint("  --- solid kinematics residual: %e \n", solid_kinematics_residual);

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

  VecDestroy(&sol_vp);
  VecDestroy(&sol_mesh);
}

void PNonlinear_FSI_Solver::GenAlpha_Seg_solve_Prestress(
    const bool &new_tangent_flag,
    const double &prestress_tol,
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
    Tissue_prestress * const &ps_ptr,
    IPGAssem * const &gassem_ptr,
    PDNSolution * const &dot_disp,
    PDNSolution * const &dot_velo,
    PDNSolution * const &dot_pres,
    PDNSolution * const &disp,
    PDNSolution * const &velo,
    PDNSolution * const &pres,
    bool &prestress_conv_flag, int &nl_counter ) const
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
  std::unique_ptr<PDNSolution> Delta_dot_disp = SYS_T::make_unique<PDNSolution_V>(pnode_v.get(), 0, false, "delta_dot_disp");

  update_solid_kinematics( -1.0 / alpha_m, dot_disp_alpha->solution, Delta_dot_disp.get() );
  update_solid_kinematics(  1.0 / alpha_m,     velo_alpha->solution, Delta_dot_disp.get() );

  // Now update displacement solutions for solid sub-domain
  dot_disp       -> PlusAX( Delta_dot_disp.get(), 1.0 );
  disp           -> PlusAX( Delta_dot_disp.get(), gamma * dt );
  dot_disp_alpha -> PlusAX( Delta_dot_disp.get(), alpha_m );
  disp_alpha     -> PlusAX( Delta_dot_disp.get(), alpha_f * gamma * dt );

  // If new_tangent_flag == TRUE, update the tangent matrix;
  // otherwise, use the matrix from the previous time step
  if( new_tangent_flag )
  {
    gassem_ptr -> Clear_KG();

    gassem_ptr->Assem_Tangent_Residual( curr_time, dt,
        dot_disp_alpha.get(), dot_velo_alpha.get(), dot_pres_alpha.get(),
        disp_alpha.get(), velo_alpha.get(), pres_alpha.get(), ps_ptr );

    SYS_T::commPrint("  --- M updated");
    lsolver->SetOperator(gassem_ptr->K);
  }
  else
  {
    gassem_ptr -> Clear_G();

    gassem_ptr->Assem_Residual( curr_time, dt,
        dot_disp_alpha.get(), dot_velo_alpha.get(), dot_pres_alpha.get(),
        disp_alpha.get(), velo_alpha.get(), pres_alpha.get(), ps_ptr );
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  SYS_T::commPrint("  Init res 2-norm: %e \n", initial_norm);

  Vec sol_vp, sol_v, sol_p;
  VecDuplicate( gassem_ptr->G, &sol_vp );

  // Now we do consistent Newton-Raphson iteration
  do
  {
    lsolver->Solve( gassem_ptr->G, sol_vp );

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

    // Assemble residual & tangent
    if( nl_counter % nrenew_freq == 0 || nl_counter >= nrenew_threshold )
    {
      gassem_ptr -> Clear_KG();

      gassem_ptr->Assem_Tangent_Residual( curr_time, dt,
          dot_disp_alpha.get(), dot_velo_alpha.get(), dot_pres_alpha.get(),
          disp_alpha.get(), velo_alpha.get(), pres_alpha.get(), ps_ptr );

      SYS_T::commPrint("  --- M updated");
      lsolver->SetOperator(gassem_ptr->K);
    }
    else
    {
      gassem_ptr -> Clear_G();

      gassem_ptr->Assem_Residual( curr_time, dt,
          dot_disp_alpha.get(), dot_velo_alpha.get(), dot_pres_alpha.get(),
          disp_alpha.get(), velo_alpha.get(), pres_alpha.get(), ps_ptr );
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

  // --------------------------------------------------------------------------
  // Calculate teh Cauchy stress in solid element and update the prestress
  gassem_ptr -> Update_Wall_Prestress( disp, pres, ps_ptr );

  const double solid_disp_norm = disp -> Norm_inf();

  SYS_T::commPrint("  --- solid disp l_inf norm: %e.\n", solid_disp_norm );

  if( solid_disp_norm < prestress_tol ) prestress_conv_flag = true;
  // --------------------------------------------------------------------------

  Print_convergence_info(nl_counter, relative_error, residual_norm);
  VecDestroy(&sol_vp);
}

// EOF
