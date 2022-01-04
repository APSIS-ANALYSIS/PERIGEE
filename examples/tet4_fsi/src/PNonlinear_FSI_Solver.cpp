#include "PNonlinear_FSI_Solver.hpp"

PNonlinear_FSI_Solver::PNonlinear_FSI_Solver(
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol, const int &input_max_iteration,
    const int &input_renew_freq )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq)
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
  SYS_T::print_sep_line();
}

void PNonlinear_FSI_Solver::rescale_inflow_value( const double &stime,
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
    const PDNSolution * const &sol_base,
    const PDNSolution * const &pre_dot_disp,
    const PDNSolution * const &pre_dot_velo,
    const PDNSolution * const &pre_dot_pres,
    const PDNSolution * const &pre_disp,
    const PDNSolution * const &pre_velo,
    const PDNSolution * const &pre_pres,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ICVFlowRate * const flr_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_v,
    const ALocal_IEN * const &lien_p,
    const FEANode * const &feanode_ptr,
    const APart_Node * const &pnode_v,
    const APart_Node * const &pnode_p,
    const ALocal_NodalBC * const &nbc_v,
    const ALocal_NodalBC * const &nbc_p,
    const ALocal_Inflow_NodalBC * const &infnbc_part,
    const ALocal_NodalBC * const &nbc_mesh,
    const ALocal_EBC * const &ebc_part,
    const ALocal_EBC * const &ebc_mesh_part,
    const IGenBC * const &gbc,
    const Matrix_PETSc * const &bc_mat,
    const Matrix_PETSc * const &bc_mesh_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const Prestress_solid * const &ps_ptr,
    IPLocAssem_2x2Block * const &lassem_fluid_ptr,
    IPLocAssem_2x2Block * const &lassem_solid_ptr,
    IPLocAssem * const &lassem_mesh_ptr,
    IPGAssem * const &gassem_ptr,
    IPGAssem * const &gassem_mesh_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PLinear_Solver_PETSc * const &lsolver_mesh_ptr,
    PDNSolution * const &dot_disp,
    PDNSolution * const &dot_velo,
    PDNSolution * const &dot_pres,
    PDNSolution * const &disp,
    PDNSolution * const &velo,
    PDNSolution * const &pres,
    bool &conv_flag, int &nl_counter ) const
{
  // Initialization
  nl_counter = 0;
  double residual_norm = 0.0, initial_norm = 0.0, relative_error = 0.0;

  const double gamma   = tmga_ptr->get_gamma();
  const double alpha_m = tmga_ptr->get_alpha_m();
  const double alpha_f = tmga_ptr->get_alpha_f();

  const double val_1 = alpha_f * gamma * dt / alpha_m;
  const double val_2 = (1.0/gamma) - (1.0/alpha_m);
  const double val_3 = 1.0 / alpha_m;

  // Same-Y predictor
  dot_disp -> Copy( pre_dot_disp ); dot_disp -> ScaleValue( (gamma-1.0)/gamma );
  dot_velo -> Copy( pre_dot_velo ); dot_velo -> ScaleValue( (gamma-1.0)/gamma );
  dot_pres -> Copy( pre_dot_pres ); dot_pres -> ScaleValue( (gamma-1.0)/gamma );

  disp -> Copy( pre_disp );
  velo -> Copy( pre_velo );
  pres -> Copy( pre_pres );

  // Define intermediate solutions
  PDNSolution dot_disp_alpha( pre_dot_disp ); 
  dot_disp_alpha.ScaleValue( 1.0 - alpha_m );
  dot_disp_alpha.PlusAX( dot_disp, alpha_m );

  PDNSolution dot_velo_alpha( pre_dot_velo );
  dot_velo_alpha.ScaleValue( 1.0 - alpha_m );
  dot_velo_alpha.PlusAX( dot_velo, alpha_m );

  PDNSolution dot_pres_alpha( pre_dot_pres );
  dot_pres_alpha.ScaleValue( 1.0 - alpha_m );
  dot_pres_alpha.PlusAX( dot_pres, alpha_m );

  PDNSolution disp_alpha( pre_disp );
  disp_alpha.ScaleValue( 1.0 - alpha_f );
  disp_alpha.PlusAX( disp, alpha_f );

  PDNSolution velo_alpha( pre_velo );
  velo_alpha.ScaleValue( 1.0 - alpha_f );
  velo_alpha.PlusAX( velo, alpha_f );

  PDNSolution pres_alpha( pre_pres );
  pres_alpha.ScaleValue( 1.0 - alpha_f );
  pres_alpha.PlusAX( pres, alpha_f );

  // Get Delta_dot_disp by assuming Delta_v is zero
  PDNSolution * Delta_dot_disp = new PDNSolution_V( pnode_v, 0, false, "delta_dot_disp" );

  update_solid_kinematics( -1.0 * val_3, pnode_v, dot_disp_alpha.solution, Delta_dot_disp );
  update_solid_kinematics(        val_3, pnode_v,     velo_alpha.solution, Delta_dot_disp );

  // Now update displacement solutions
  dot_disp->PlusAX( Delta_dot_disp, 1.0 );
  disp    ->PlusAX( Delta_dot_disp, gamma * dt );
  dot_disp_alpha.PlusAX( Delta_dot_disp, alpha_m );
  disp_alpha.PlusAX(     Delta_dot_disp, alpha_f * gamma * dt );
  
  // Update inflow boundary values
  rescale_inflow_value( curr_time + dt,           infnbc_part, flr_ptr, sol_base, velo );
  rescale_inflow_value( curr_time + alpha_f * dt, infnbc_part, flr_ptr, sol_base, &velo_alpha );






















  delete Delta_dot_disp;
}

void PNonlinear_FSI_Solver::update_solid_kinematics( 
    const double &val, const APart_Node * const &pnode,
    const Vec &input, PDNSolution * const &output ) const
{
  const int nlocal = pnode -> get_nlocalnode_solid();
  
  Vec local_input, local_output;

  VecGhostGetLocalForm(input,  &local_input);
  VecGhostGetLocalForm(output->solution, &local_output);

  double * array_input, * array_output;

  VecGetArray(local_input, &array_input);
  VecGetArray(local_output, &array_output);

  for(int ii=0; ii<nlocal; ++ii)
  {
    const int ii3 = pnode->get_node_loc_solid(ii) * 3;

    array_output[ii3  ] += val * array_input[ii3  ];
    array_output[ii3+1] += val * array_input[ii3+1];
    array_output[ii3+2] += val * array_input[ii3+2];
  }

  VecRestoreArray(local_input, &array_input);
  VecRestoreArray(local_output, &array_output);

  VecGhostRestoreLocalForm(input, &local_input);
  VecGhostRestoreLocalForm(output->solution, &local_output);

  output->GhostUpdate(); // update the ghost slots
}






// EOF
