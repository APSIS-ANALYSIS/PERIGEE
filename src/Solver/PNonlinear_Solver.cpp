#include "PNonlinear_Solver.hpp"

PNonlinear_Solver::PNonlinear_Solver(
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol,
    const int &input_max_iteration, const int &input_renew_freq )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq)
{
}


PNonlinear_Solver::~PNonlinear_Solver()
{}


void PNonlinear_Solver::Info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "relative tolerance: %e \n", nr_tol);
  PetscPrintf(PETSC_COMM_WORLD, "absolute tolerance: %e \n", na_tol);
  PetscPrintf(PETSC_COMM_WORLD, "divergence tolerance: %e \n", nd_tol);
  PetscPrintf(PETSC_COMM_WORLD, "maximum iteration: %d \n", nmaxits);
  PetscPrintf(PETSC_COMM_WORLD, "tangent matrix renew frequency: %d \n", nrenew_freq);
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PNonlinear_Solver::Gen_alpha_solve(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &pre_velo,
    const PDNSolution * const &pre_disp,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const IALocal_BC * const &bc_part,
    const AInt_Weight * const &wei_ptr,
    const std::vector<FEAElement *> &ele_ptr,
    IPLocAssem * const &lassem_ptr,
    PGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PDNSolution * const &velo,
    PDNSolution * const &disp,
    bool &conv_flag,
    int &nl_counter ) const
{
  nl_counter = 0;
  double ksp_its_num;
  double ksp_max_its_num = (double) lsolver_ptr->get_ksp_maxits();
  double ksp_its_check = 0.3 * ksp_max_its_num;
  double residual_norm = 0.0;
  double initial_norm = 0.0;
  double relative_error = 0.0;

  const double gamma   = tmga_ptr->get_gamma();
  const double alpha_m = tmga_ptr->get_alpha_m();
  const double alpha_f = tmga_ptr->get_alpha_f();

  // predictor
  disp->ScaleValue(0.0);
  velo->ScaleValue(0.0);
  
  PDNSolution velo_step(*velo);
  
  disp->PlusAX(*pre_disp, 1.0);
  velo->PlusAX(*pre_velo, (gamma-1.0)/gamma);

  // define displacement at alpha_f and velocity at alpha_m
  PDNSolution disp_alpha(*pre_disp);
  PDNSolution velo_alpha(*pre_velo);

  disp_alpha.ScaleValue((1.0 - alpha_f));
  velo_alpha.ScaleValue((1.0 - alpha_m));
  
  disp_alpha.PlusAX(*disp, alpha_f);
  velo_alpha.PlusAX(*velo, alpha_m);
  
  // if new tanget flag is true, update the tangent matrix,
  // otherwise, keep using the tangent matrix from the previous
  // time step
  if( new_tangent_flag )
  {
    gassem_ptr->Clear_KG();
    gassem_ptr->Assem_tangent_residual( &velo_alpha, &disp_alpha,
     curr_time, dt, alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
     feanode_ptr, wei_ptr, ele_ptr, bc_part );
    lsolver_ptr->SetOperator(gassem_ptr->K);
  }
  else
  {
    gassem_ptr->Clear_G();
    gassem_ptr->Assem_residual( &velo_alpha, &disp_alpha,
     curr_time, dt, alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
     feanode_ptr, wei_ptr, ele_ptr, bc_part );
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  PetscPrintf(PETSC_COMM_WORLD, "  Initial residual 2-norm: %e \n", initial_norm);

  do
  {
    lsolver_ptr->Solve( gassem_ptr->G, &velo_step );
    
    nl_counter += 1;

    ksp_its_num = (double) lsolver_ptr->get_ksp_it_num();
   
    if(ksp_its_num > ksp_its_check)
    {
      PetscPrintf(PETSC_COMM_WORLD, "Warning: linear solver converges slowly! \n");
      break;
    }

    // corrector
    disp->PlusAX( velo_step, -1.0 * gamma * dt );
    velo->PlusAX( velo_step, -1.0 );

    disp_alpha.PlusAX(velo_step, -1.0 * alpha_f * gamma * dt);
    velo_alpha.PlusAX(velo_step, -1.0 * alpha_m);

    if(nl_counter % nrenew_freq == 0)
    {
      gassem_ptr->Clear_KG();
      gassem_ptr->Assem_tangent_residual( &velo_alpha, &disp_alpha,
          curr_time, dt, alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
          feanode_ptr, wei_ptr, ele_ptr, bc_part );
      lsolver_ptr->SetOperator(gassem_ptr->K);
    }
    else
    {
      gassem_ptr->Clear_G();
      gassem_ptr->Assem_residual( &velo_alpha, &disp_alpha,
          curr_time, dt, alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
          feanode_ptr, wei_ptr, ele_ptr, bc_part );
    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    PetscPrintf(PETSC_COMM_WORLD, "  --- residual norm: %e \n",
                residual_norm);

    relative_error = residual_norm / initial_norm;

    if( relative_error >= nd_tol )
    {
      PetscPrintf(PETSC_COMM_WORLD, "Warning: nonlinear solver is diverging with error %e \n",
          relative_error);
      break;
    }

  }while(nl_counter < nmaxits && relative_error > nr_tol && 
      residual_norm > na_tol );

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol)
    conv_flag = true;
  else
    conv_flag = false;
}


void PNonlinear_Solver::NewtonRaphson_solve(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &curr,
    PDNSolution * const &step,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const IALocal_BC * const &bc_part,
    const AInt_Weight * const &wei_ptr,
    const std::vector<FEAElement *> &ele_ptr,
    IPLocAssem * const &lassem_ptr,
    PGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PDNSolution * const &next,
    bool &conv_flag, 
    int &nl_counter ) const
{
  nl_counter = 0;
  //double ksp_its_num;
  //double ksp_max_its_num = (double) lsolver_ptr->get_ksp_maxits();
  //double ksp_its_check = 0.3 * ksp_max_its_num;
  double residual_norm = 0.0;
  double initial_norm = 0.0;
  double relative_error = 0.0;

  // set next as curr as an initial guess
  //next->ScaleValue(0.0);
  //next->PlusAX(*curr, 1.0); 
  next->Copy(*curr);

  // if new tanget flag is true, update the tangent matrix,
  // otherwise, keep using the tangent matrix from the previous
  // time step
  if( new_tangent_flag )
  {
    gassem_ptr->Clear_KG();
    gassem_ptr->Assem_tangent_residual( curr, next,
     curr_time, dt, alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
     feanode_ptr, wei_ptr, ele_ptr, bc_part );
    PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
    lsolver_ptr->SetOperator(gassem_ptr->K);
  }
  else
  {
    gassem_ptr->Clear_G();
    gassem_ptr->Assem_residual( curr, next,
     curr_time, dt, alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
     feanode_ptr, wei_ptr, ele_ptr, bc_part );
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  PetscPrintf(PETSC_COMM_WORLD, "  Init res 2-norm: %e \n", initial_norm);

  //lsolver_ptr->SetOperator(gassem_ptr->K);

  do
  {
    lsolver_ptr->Solve( gassem_ptr->G, step );
    
    nl_counter += 1;

    //ksp_its_num = (double) lsolver_ptr->get_ksp_it_num();
   
    //if(ksp_its_num > ksp_its_check)
    //{
    //  PetscPrintf(PETSC_COMM_WORLD, "Warning: linear solver converges slowly! \n");
    //  break;
    //}

    // update next
    next->PlusAX( *step, -1.0 );

    if(nl_counter > nrenew_freq)
    {
      gassem_ptr->Clear_KG();
      gassem_ptr->Assem_tangent_residual( curr, next,
          curr_time, dt, alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
          feanode_ptr, wei_ptr, ele_ptr, bc_part );
      lsolver_ptr->SetOperator(gassem_ptr->K);
      PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
    }
    else
    {
      gassem_ptr->Clear_G();
      gassem_ptr->Assem_residual( curr, next,
          curr_time, dt, alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
          feanode_ptr, wei_ptr, ele_ptr, bc_part );
    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    PetscPrintf(PETSC_COMM_WORLD, "  --- res norm: %e \n",
                residual_norm);

    relative_error = residual_norm / initial_norm;

    if( relative_error >= nd_tol )
    {
      PetscPrintf(PETSC_COMM_WORLD, "Warning: nonlinear solver is diverging with error %e \n",
          relative_error);
      break;
    }

  }while(nl_counter < nmaxits && relative_error > nr_tol && 
      residual_norm > na_tol );

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol)
    conv_flag = true;
  else
    conv_flag = false;
}


void PNonlinear_Solver::Gen_alpha_VMS_solve(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &pre_velo,
    const PDNSolution * const &pre_disp,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const IALocal_BC * const &bc_part,
    const AInt_Weight * const &wei_ptr,
    const std::vector<FEAElement *> &ele_ptr,
    IPLocAssem * const &lassem_ptr,
    PGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PDNSolution * const &velo,
    PDNSolution * const &disp,
    bool &conv_flag,
    int &nl_counter ) const
{
  nl_counter = 0;
  double ksp_its_num;
  double ksp_max_its_num = (double) lsolver_ptr->get_ksp_maxits();
  double ksp_its_check = 0.7 * ksp_max_its_num;
  double residual_norm = 0.0;
  double initial_norm = 0.0;
  double relative_error = 0.0;

  const int dof_per_node = anode_ptr->get_dof();
  const int dof_minus_one = dof_per_node - 1;

  const double gamma   = tmga_ptr->get_gamma();
  const double alpha_m = tmga_ptr->get_alpha_m();
  const double alpha_f = tmga_ptr->get_alpha_f();

  // predictor
  disp->ScaleValue(0.0);
  velo->ScaleValue(0.0);
  
  PDNSolution velo_step(*velo);
  
  disp->PlusAX(*pre_disp, 1.0);
  velo->PlusAX(*pre_velo, (gamma-1.0)/gamma);

  // define displacement at alpha_f and velocity at alpha_m
  PDNSolution velo_alpha(*pre_velo);
  velo_alpha.ScaleValue((1.0 - alpha_m));
  velo_alpha.PlusAX(*velo, alpha_m);

  PDNSolution disp_alpha(*pre_disp);
  
  // if new tanget flag is true, update the tangent matrix,
  // otherwise, keep using the tangent matrix from the previous
  // time step
  if( new_tangent_flag )
  {
    gassem_ptr->Clear_KG();
    gassem_ptr->Assem_tangent_residual( &velo_alpha, &disp_alpha,
     curr_time, dt, alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
     feanode_ptr, wei_ptr, ele_ptr, bc_part );
    lsolver_ptr->SetOperator(gassem_ptr->K);
  }
  else
  {
    gassem_ptr->Clear_G();
    gassem_ptr->Assem_residual( &velo_alpha, &disp_alpha,
     curr_time, dt, alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
     feanode_ptr, wei_ptr, ele_ptr, bc_part );
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  PetscPrintf(PETSC_COMM_WORLD, "  Initial residual 2-norm: %e \n", initial_norm);

  do
  {
    lsolver_ptr->Solve( gassem_ptr->G, &velo_step );
    
    nl_counter += 1;

    ksp_its_num = (double) lsolver_ptr->get_ksp_it_num();
   
    if(ksp_its_num > ksp_its_check)
    {
      PetscPrintf(PETSC_COMM_WORLD, "Warning: linear solver converges slowly! \n");
      break;
    }

    // corrector
    disp->PlusAiX( velo_step, -1.0 * gamma * dt, -1.0, dof_minus_one, 1);
    velo->PlusAX( velo_step, -1.0 );

    disp_alpha.PlusAiX(velo_step, -1.0 * alpha_f * gamma * dt, -1.0, dof_minus_one, 1);
    velo_alpha.PlusAX(velo_step, -1.0 * alpha_m);

    if(nl_counter % nrenew_freq == 0)
    {
      gassem_ptr->Clear_KG();
      gassem_ptr->Assem_tangent_residual( &velo_alpha, &disp_alpha,
          curr_time, dt, alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
          feanode_ptr, wei_ptr, ele_ptr, bc_part );
      lsolver_ptr->SetOperator(gassem_ptr->K);
    }
    else
    {
      gassem_ptr->Clear_G();
      gassem_ptr->Assem_residual( &velo_alpha, &disp_alpha,
          curr_time, dt, alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
          feanode_ptr, wei_ptr, ele_ptr, bc_part );
    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    
    PetscPrintf(PETSC_COMM_WORLD, " --- nl_res_norm: %e \n", residual_norm);

    relative_error = residual_norm / initial_norm;

    if( relative_error >= nd_tol )
    {
      PetscPrintf(PETSC_COMM_WORLD, "Warning: nonlinear solver is diverging with error %e \n",
          relative_error);
      break;
    }

  }while(nl_counter < nmaxits && relative_error > nr_tol && 
      residual_norm > na_tol );

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol)
    conv_flag = true;
  else
    conv_flag = false;
}


void PNonlinear_Solver::GenAlpha_solve(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &pre_velo,
    const PDNSolution * const &pre_disp,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_EBC * const &ebc_part,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PDNSolution * const &velo,
    PDNSolution * const &disp,
    bool &conv_flag, int &nl_counter ) const
{
  // Set counter and error-related quantities
  nl_counter = 0;
  double residual_norm = 0.0, initial_norm = 0.0, relative_error = 0.0;

  // Gen-alpha parameters
  const double gamma   = tmga_ptr->get_gamma();
  const double alpha_m = tmga_ptr->get_alpha_m();
  const double alpha_f = tmga_ptr->get_alpha_f();

  // Predictor stage
  disp->Copy(*pre_disp);
  velo->Copy(*pre_velo);
  velo->ScaleValue( (gamma-1.0)/gamma );

  PDNSolution velo_step(*velo);
  
  // Define the acce at alpha_m and disp at alpha_f
  PDNSolution velo_alpha(*pre_velo);
  velo_alpha.ScaleValue(1.0 - alpha_m);
  velo_alpha.PlusAX(*velo, alpha_m);

  PDNSolution disp_alpha(*pre_disp);
  disp_alpha.ScaleValue(1.0 - alpha_f);
  disp_alpha.PlusAX(*disp, alpha_f);

  // If new_tangent_flag == true, update the tangent matrix,
  // otherwise, use the tangent matrix from the previous time step.
  if( new_tangent_flag )
  {
    gassem_ptr->Clear_KG();
    gassem_ptr->Assem_tangent_residual( &velo_alpha, &disp_alpha,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements, 
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part );
    PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
    lsolver_ptr->SetOperator(gassem_ptr->K);
  }
  else
  {
    gassem_ptr->Clear_G();
    gassem_ptr->Assem_residual( &velo_alpha, &disp_alpha,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements, 
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part );
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  PetscPrintf(PETSC_COMM_WORLD, "  Init res 2-norm: %e \n", initial_norm);

  do
  {
    lsolver_ptr->Solve( gassem_ptr->G, &velo_step );

    // bc_mat is a square matrix that maps the dir node to 0 and slave node
    // to its master node's value
    bc_mat->MatMultSol( &velo_step );

    nl_counter += 1;

    // Corrector
    velo->PlusAX( velo_step, -1.0 );
    disp->PlusAX( velo_step, (-1.0) * gamma * dt );

    velo_alpha.PlusAX( velo_step, (-1.0) * alpha_m );
    disp_alpha.PlusAX( velo_step, (-1.0) * alpha_f * gamma * dt );

    if(nl_counter >= nrenew_freq)
    {
      gassem_ptr->Clear_KG();
      gassem_ptr->Assem_tangent_residual( &velo_alpha, &disp_alpha,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements, 
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part );
      lsolver_ptr->SetOperator(gassem_ptr->K);
      PetscPrintf(PETSC_COMM_WORLD, "  --- M updated");
    }
    else
    {
      gassem_ptr->Clear_G();
      gassem_ptr->Assem_residual( &velo_alpha, &disp_alpha,
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

  }while(nl_counter < nmaxits && relative_error > nr_tol && residual_norm > na_tol );

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol) conv_flag = true;
  else conv_flag = false;
}


// EOF
