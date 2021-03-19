#include "PNonlinear_CMM_Solver.hpp"

PNonlinear_CMM_Solver::PNonlinear_CMM_Solver(
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


PNonlinear_CMM_Solver::~PNonlinear_CMM_Solver()
{
  delete dot_step; dot_step = nullptr;
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
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Inflow_NodalBC * const &infnbc_part,
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
  update_wall(-1.0, &sol_alpha, &G_kinematic, ebc_wall_part);

  // ==== WOMERSLEY CHANGES BEGIN ====
  update_nodal_bc(curr_time + dt, nbc_part, anode_ptr, feanode_ptr, sol);
  update_nodal_bc(curr_time + alpha_f * dt, nbc_part, anode_ptr, feanode_ptr, &sol_alpha);

  // ------------------------------------------------- 
  // // Update the inflow boundary values
  // rescale_inflow_value(curr_time+dt, infnbc_part, flr_ptr, sol_base, sol);
  // rescale_inflow_value(curr_time+alpha_f*dt, infnbc_part, flr_ptr, sol_base, &sol_alpha);
  // // ------------------------------------------------- 
  // ==== WOMERSLEY CHANGES END ====

  // **** PRESTRESS TODO: if (prestress_flag), set wall disp and velo to zero

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
        elementw, quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part, ebc_wall_part, gbc );
   
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
        elementw, quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part, ebc_wall_part, gbc );

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

    // Update dot_sol, dot_sol_wall_disp
    dot_sol->PlusAX( dot_step, -1.0 );
    update_wall( (-1.0) * alpha_f * gamma * dt / alpha_m,
        dot_step, dot_sol_wall_disp, ebc_wall_part ); 
    dot_sol_wall_disp->PlusAX(G_kinematic, (-1.0) / alpha_m );

    // Update sol, sol_wall_disp
    sol->PlusAX( dot_step, (-1.0) * gamma * dt );
    update_wall( (-1.0) * alpha_f * gamma * gamma * dt * dt / alpha_m,
        dot_step, sol_wall_disp, ebc_wall_part ); 
    sol_wall_disp->PlusAX(G_kinematic, (-1.0) * gamma * dt / alpha_m );

    // Update dol_sol at alpha_m: dot_sol_alpha, dot_wall_disp_alpha
    dot_sol_alpha.PlusAX( dot_step, (-1.0) * alpha_m );
    update_wall( (-1.0) * alpha_f * gamma * dt,
        dot_step, &dot_wall_disp_alpha, ebc_wall_part ); 
    dot_wall_disp_alpha.PlusAX(G_kinematic, -1.0 );

    // Update sol at alpha_f: sol_alpha, wall_disp_alpha
    sol_alpha.PlusAX( dot_step, (-1.0) * alpha_f * gamma * dt );
    update_wall( (-1.0) * alpha_f * alpha_f * gamma * gamma * dt * dt / alpha_m,
        dot_step, &wall_disp_alpha, ebc_wall_part ); 
    wall_disp_alpha.PlusAX(G_kinematic, (-1.0) * alpha_f * gamma * dt / alpha_m );

    // Update kinematic residual = dot_wall_disp_alpha - velo_alpha
    G_kinematic.Copy(dot_wall_disp_alpha);
    update_wall(-1.0, &sol_alpha, &G_kinematic, ebc_wall_part);

    // Assembly residual (& tangent if condition satisfied) 
    if( nl_counter % nrenew_freq == 0 || nl_counter >= nrenew_threshold )
    {
      gassem_ptr->Clear_KG();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(mat_assem_1_event, 0,0,0,0);
#endif

      gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, &wall_disp_alpha,
          dot_sol, sol, curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          elementw, quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part, ebc_wall_part, gbc );

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
          elementw, quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part, ebc_wall_part, gbc );

#ifdef PETSC_USE_LOG
      PetscLogEventEnd(vec_assem_1_event,0,0,0,0);
#endif
    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    
    SYS_T::commPrint("  --- nl_res: %e \n", residual_norm);

    // DEBUGGING ISL 3/6/2021
    double kinematic_residual_norm = 0.0;
    VecNorm(G_kinematic.solution, NORM_2, &kinematic_residual_norm);
    SYS_T::commPrint("  --- kinematic_res: %e \n", kinematic_residual_norm);

    relative_error = residual_norm / initial_norm;

    if( relative_error >= nd_tol )
      SYS_T::print_fatal( "Error: nonlinear solver is diverging with error %e. Job killed.\n", relative_error);

  }while(nl_counter<nmaxits && relative_error > nr_tol && residual_norm > na_tol);

  // **** PRESTRESS TODO: if (prestress_flag), call PGAssem_Tet_CMM_GenAlpha::Update_Wall_Prestress() 
  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol) conv_flag = true;
  else conv_flag = false;
}


void PNonlinear_CMM_Solver::rescale_inflow_value( const double &stime,
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


// ==== WOMERSLEY CHANGES BEGIN ====
void PNonlinear_CMM_Solver::update_nodal_bc( const double &stime,
    const ALocal_NodalBC * const &nbc_part,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    PDNSolution * const &sol ) const
{
  // Verify that the dof of sol is 4
  SYS_T::print_fatal_if(sol->get_dof_num() != 4,
      "Error in PNonlinear_CMM_Solver::update_nodal_bc: incorrect dimension of sol. \n");

  const double rho    = 1.0;                                                 // fluid density
  const double vis_mu = 4.0e-2;                                              // fluid viscosity

  const double R     = 0.3;                                                  // pipe radius
  const double omega = MATH_T::PI * 2.0 / 1.1;                               // freqency
  const std::complex<double> i1(0.0, 1.0);
  const std::complex<double> i1_1d5(-0.707106781186547, 0.707106781186547);
  const auto Omega   = std::sqrt(rho * omega / vis_mu) * R;                  // womersley number 
  const auto Lambda  = i1_1d5 * Omega;

  const double k0 = -21.0469;                                                // mean pressure gradient
  const std::complex<double> B1(-4.926286624202966e3, -4.092542965905093e3); // pressure Fourier coeff
  const std::complex<double> c1(8.863128942479001e2,   2.978553160539686e1); // wave speed
  const std::complex<double> G1(0.829733473284180,      -0.374935589823809); // elasticity factor

  const int nlocalnode = anode_ptr->get_nlocalnode();

  std::vector<int> LDN; 

  // Update velocities
  for(int ii=1; ii<4; ++ii)
  {
    const int num_LD = nbc_part -> get_Num_LD(ii);

    // Store a copy of LDN
    LDN.resize(num_LD); 

    for(int jj=0; jj<num_LD; ++jj)
    {
      LDN[jj] = nbc_part -> get_LDN(ii, jj);
    }

    for(int jj=0; jj<nlocalnode; ++jj)
    {
      const int node_index = anode_ptr->get_node_loc(jj);

      if( VEC_T::is_invec(LDN, node_index) )
      {
        const double x = feanode_ptr->get_ctrlPts_x(jj);
        const double y = feanode_ptr->get_ctrlPts_y(jj);
        const double z = feanode_ptr->get_ctrlPts_z(jj);

        const double r  = std::sqrt(x*x + y*y);
        const auto   xi = Lambda * r / R;

        const auto bes0_xi     = sp_bessel::besselJ(0, xi);
        const auto bes1_xi     = sp_bessel::besselJ(1, xi);
        const auto bes0_Lambda = sp_bessel::besselJ(0, Lambda);

        // axial velocity
        const double w = k0 * (x*x + y*y - R*R) / (4.0*vis_mu)
            + std::real( B1 / (rho * c1) * (1.0 - G1 * bes0_xi / bes0_Lambda) * exp(i1*omega*(stime-z/c1)) );

        // radial velocity
        const double vr = std::real( i1 * omega * R * B1 / ( 2.0 * rho * c1 * c1 )
            * ( r / R - 2.0 * G1 * bes1_xi / (Lambda * bes0_Lambda) ) * exp(i1*omega*(stime-z/c1)) );

        // polar to cartesian transformation
        const double theta = std::atan2(y, x);

        VecSetValue(sol->solution, node_index*4+1, vr * std::cos(theta), INSERT_VALUES);
        VecSetValue(sol->solution, node_index*4+2, vr * std::sin(theta), INSERT_VALUES);
        VecSetValue(sol->solution, node_index*4+3, w, INSERT_VALUES);
      }
    } // end nlocalnode loop
  }

  VecAssemblyBegin(sol->solution); VecAssemblyEnd(sol->solution);
  sol->GhostUpdate();
}
// ==== WOMERSLEY CHANGES END ====


void PNonlinear_CMM_Solver::update_wall( const double &val,
    const PDNSolution * const &dot_step,
    PDNSolution * const &wall_data,
    const ALocal_EBC * const &ebc_wall_part ) const
{
  // Verify that the dof of dot_step is 4
  SYS_T::print_fatal_if(dot_step->get_dof_num() != 4,
      "Error in PNonlinear_CMM_Solver::update_dot_wall_disp: incorrect dimension of dot_step. \n");

  // Verify that the dof of wall_data is 3
  SYS_T::print_fatal_if(wall_data->get_dof_num() != 3,
      "Error in PNonlinear_CMM_Solver::update_dot_wall_disp: incorrect dimension of wall_data. \n");

  // Verify consistency in the number of local nodes
  SYS_T::print_fatal_if( !is_layout_equal(*dot_step, *wall_data), "Error in PNonlinear_CMM_Solver::update_dot_wall_disp: solution vector layout mismatch between dot_step and wall_data. \n");

  // wall has only one surface per the assumption in wall ebc
  const int ebc_id = 0;
  const int num_snode = ebc_wall_part->get_num_local_node(ebc_id);

  Vec ldotstep, lwalldata;
  double * array_dotstep, * array_walldata;

  VecGhostGetLocalForm(dot_step->solution, &ldotstep);
  VecGhostGetLocalForm(wall_data->solution, &lwalldata);

  VecGetArray(ldotstep, &array_dotstep);
  VecGetArray(lwalldata, &array_walldata);

  for(int ii=0; ii<num_snode; ++ii)
  {
    const int pos = ebc_wall_part->get_local_node_pos(ebc_id, ii);

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

// EOF
