#include "PTime_NS_HERK_Solver.hpp"

PTime_NS_HERK_Solver::PTime_NS_HERK_Solver(
    std::unique_ptr<PGAssem_Block_NS_FEM_HERK> in_gassem,
    std::unique_ptr<PLinear_Solver_PETSc> in_lsolver,
    std::unique_ptr<Matrix_PETSc> in_bc_mat,
    std::unique_ptr<ITimeMethod_RungeKutta> in_tmRK,
    std::unique_ptr<IFlowRate> in_flrate,
    std::unique_ptr<IFlowRate> in_dot_flrate,
    std::unique_ptr<PDNSolution> in_sol_base,
    std::unique_ptr<ALocal_InflowBC> in_infnbc, 
    const std::string &input_name, const int &in_nlocalnode,
    const int &input_record_freq, const double &input_final_time )
: final_time(input_final_time), sol_record_freq(input_record_freq),
  pb_name(input_name), nlocalnode(in_nlocalnode), gassem(std::move(in_gassem)), 
  lsolver(std::move(in_lsolver)), bc_mat(std::move(in_bc_mat)), 
  tmRK(std::move(in_tmRK)), flrate(std::move(in_flrate)), 
  dot_flrate(std::move(in_dot_flrate)), sol_base(std::move(in_sol_base)),
  infnbc(std::move(in_infnbc))
{}

std::string PTime_NS_HERK_Solver::Name_Generator(const int &counter) const
{
  std::ostringstream temp;
  temp.str("");
  temp<<900000000 + counter;

  return pb_name + temp.str();
}

void PTime_NS_HERK_Solver::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("Time stepping solver setted up:\n");
  SYS_T::commPrint("  final time: %e \n", final_time);
  SYS_T::commPrint("  solution record frequency : %d \n", sol_record_freq);
  SYS_T::commPrint("  solution base name: %s \n", pb_name.c_str());
  SYS_T::commPrint("----------------------------------------------------------- \n");
}

void PTime_NS_HERK_Solver::Write_restart_file(const PDNTimeStep * const &timeinfo,
    const std::string &solname ) const
{
  std::ofstream restart_file("restart_file.txt", std::ofstream::out | std::ofstream::trunc);
  if( restart_file.is_open() )
  {
    restart_file<<timeinfo->get_index()<<std::endl;
    restart_file<<timeinfo->get_time()<<std::endl;
    restart_file<<timeinfo->get_step()<<std::endl;
    restart_file<<solname.c_str()<<std::endl;
    restart_file.close();
  }
  else
    SYS_T::print_fatal("Error: PTimeSolver cannot open restart_file.txt");
}

void PTime_NS_HERK_Solver::TM_NS_HERK(
    const bool &restart_init_assembly_flag, 
    std::unique_ptr<PDNSolution> init_sol,
    std::unique_ptr<PDNSolution> init_velo,
    std::unique_ptr<PDNSolution> init_dot_velo,
    std::unique_ptr<PDNSolution> init_pres,
    std::unique_ptr<PDNTimeStep> time_info ) const
{
  const int ss = tmRK->get_RK_step();

  // The velo solutions in each sub-step at the (n+1)-th time step
  PDNSolution** cur_velo_sols = new PDNSolution*[ss];
  for(int ii = 0; ii < ss; ++ii)
    cur_velo_sols[ii] = new PDNSolution(*init_velo);

  // The velo solution in the final step at the (n+1)-th time step
  PDNSolution * cur_velo = new PDNSolution(*init_velo);

  // The dot_velo solution in the final step at the (n+1)-th time step
  PDNSolution * cur_dot_velo = new PDNSolution(*init_dot_velo);

  // The pres solutions in each sub-step at the (n+1)-th time step
  PDNSolution** cur_pres_sols = new PDNSolution*[ss];
  for(int ii = 0; ii < ss; ++ii)
    cur_pres_sols[ii] = new PDNSolution(*init_pres);

  // The pres solution in the final step at the (n+1)-th time step
  PDNSolution * cur_pres = new PDNSolution(*init_pres);

  // The solution in the final step at the (n+1)-th time step
  PDNSolution * cur_sol = new PDNSolution(*init_sol);

  // The velo solutions in each sub-step at the n-th time step
  PDNSolution** pre_velo_sols = new PDNSolution*[ss];
  for(int ii = 0; ii < ss; ++ii)
    pre_velo_sols[ii] = new PDNSolution(*init_velo);

  // The pres solutions in each sub-step at the n-th time step
  PDNSolution** pre_pres_sols = new PDNSolution*[ss];
  for(int ii = 0; ii < ss; ++ii)
    pre_pres_sols[ii] = new PDNSolution(*init_pres);    

  // The velo solution in the final step at the n-th time step
  PDNSolution * pre_velo = new PDNSolution(*init_velo);

  // The pres solution in the final step at the n-th time step
  PDNSolution * pre_pres = new PDNSolution(*init_pres);

  // The velo solution in the final step at the (n-1)-th time step
  PDNSolution * pre_velo_before = new PDNSolution(*init_velo);

  // If this is a restart run, do not re-write the solution binaries
  if(restart_init_assembly_flag == false)
  {
    const auto sol_name = Name_Generator(time_info->get_index());
    cur_sol->WriteBinary(sol_name);
  }

  SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
      time_info->get_time(), time_info->get_step(), time_info->get_index(),
      SYS_T::get_time().c_str());

  while (time_info->get_time() < final_time)
  {   
    HERK_Solve_NS(
        time_info->get_time(), time_info->get_step(),
        cur_velo_sols, cur_velo, cur_dot_velo, 
        cur_pres_sols, cur_pres, pre_velo_sols, pre_velo,
        pre_pres_sols, pre_pres, pre_velo_before, cur_sol );

    // Update the time step information
    time_info->TimeIncrement();

    SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index(),
        SYS_T::get_time().c_str());

    // Record solution if meets criteria
    if( time_info->get_index()%sol_record_freq == 0 )
    {
      const auto sol_name = Name_Generator( time_info->get_index() );
      cur_sol->WriteBinary(sol_name);
    }

    // Prepare for next time step
    pre_velo_before->Copy(pre_velo);
    pre_velo->Copy(cur_velo);
    pre_pres->Copy(cur_pres);

    for(int ii = 0; ii < ss; ++ii)
    {
      pre_velo_sols[ii]->Copy(cur_velo_sols[ii]);
      pre_pres_sols[ii]->Copy(cur_pres_sols[ii]);
    }  
  }

  for (int ii = 0; ii < ss; ++ii) 
  {
      delete cur_velo_sols[ii]; delete pre_velo_sols[ii]; 
      delete cur_pres_sols[ii]; delete pre_pres_sols[ii];
  }
  delete[] cur_velo_sols; delete[] pre_velo_sols; delete[] cur_pres_sols; delete[] pre_pres_sols;
  delete cur_velo; delete cur_dot_velo; delete cur_pres; delete cur_sol; 
  delete pre_velo; delete pre_pres; delete pre_velo_before;
}

void PTime_NS_HERK_Solver::HERK_Solve_NS(
    const double &curr_time, const double &dt,
    PDNSolution ** const &cur_velo_sols,
    PDNSolution * const &cur_velo,
    PDNSolution * const &cur_dot_velo,
    PDNSolution ** const &cur_pres_sols,
    PDNSolution * const &cur_pres,
    PDNSolution ** const &pre_velo_sols,
    PDNSolution * const &pre_velo,
    PDNSolution ** const &pre_pres_sols,
    PDNSolution * const &pre_pres,
    PDNSolution * const &pre_velo_before,
    PDNSolution * const &cur_sol ) const
{ 
  auto dot_step = SYS_T::make_unique<PDNSolution>( cur_sol );

  // HERK's number of steps
  const int ss = tmRK->get_RK_step();
  
  // The first sub-step, u_1 = u_n
  cur_velo_sols[0] -> Copy(*pre_velo);
    
  SYS_T::commPrint(" ==> Start solving the SubStep -- substep = 1 is solved. \n");
  
  // Sub-step (starting from the second sub-step)
  for(int ii = 1; ii < ss; ++ii)
  {
    // Make the velo in each sub step meet the Dirchlet boundary
    rescale_inflow_velo(curr_time + tmRK->get_RK_c(ii) * dt, cur_velo_sols[ii]);
  
    gassem->Clear_G();  // K use Matrix-free
      
    gassem->Assem_tangent_residual_substep( ii, cur_velo_sols, cur_pres_sols,
      pre_velo_sols, pre_velo, pre_pres_sols, pre_velo_before, tmRK.get(), 
      curr_time, dt );     
    
    // lsolver->SetOperator(gassem->K);
    lsolver->Solve(gassem->G, dot_step.get());
  
    bc_mat->MatMultSol( dot_step.get() );
  
    SYS_T::commPrint(" --- substep = %d is solved. \n", ii+1);
  
    Update_pressure_velocity(cur_velo_sols[ii], cur_pres_sols[ii-1], dot_step.get());
  }
  // Final step
  SYS_T::commPrint(" ==> Start solving the LastStep: \n");
  
   // Make the velo in the last step meet the Dirchlet boundary
   rescale_inflow_velo(curr_time + dt, cur_velo);
  
   gassem->Clear_G();
  
   gassem->Assem_tangent_residual_finalstep( cur_velo_sols, cur_velo, 
      cur_pres_sols, pre_velo_sols, pre_velo, pre_pres_sols, pre_velo_before,
      tmRK.get(), curr_time, dt );

    // lsolver->SetOperator(gassem->K);
    lsolver->Solve(gassem->G, dot_step.get());
  
    bc_mat->MatMultSol( dot_step.get() );
  
    SYS_T::commPrint(" --- laststep is solved. \n");
  
    Update_pressure_velocity(cur_velo, cur_pres_sols[ss-1], dot_step.get());
  
    // Pressure stage
    SYS_T::commPrint(" ==> Start solving the FinalStep: \n");
  
    //Make the dot_velo in the final step meet the Dirchlet boundary
    rescale_inflow_dot_velo(curr_time + dt, cur_dot_velo);
  
    gassem->Clear_G();
  
    gassem->Assem_tangent_residual_presstage( cur_dot_velo, cur_velo_sols, 
      cur_velo, cur_pres_sols, pre_velo, cur_pres, tmRK.get(), curr_time, dt );
  
    // lsolver->SetOperator(gassem->K);
    lsolver->Solve(gassem->G, dot_step.get());
  
    bc_mat->MatMultSol( dot_step.get() );
  
    SYS_T::commPrint(" --- finalstep is solved. \n");
  
    Update_pressure_velocity(cur_dot_velo, cur_pres, dot_step.get());
  
    // Assemble velo and pres at the (n+1)-th time step into a solution vector
    Update_solutions(cur_velo, cur_pres, cur_sol);
}

void PTime_NS_HERK_Solver::rescale_inflow_velo( const double &stime,
    PDNSolution * const &velo ) const
{
  const int num_nbc = infnbc -> get_num_nbc();

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    const int numnode = infnbc -> get_Num_LD( nbc_id );

    const double factor  = flrate -> get_flow_rate( nbc_id, stime );
    const double std_dev = flrate -> get_flow_TI_std_dev( nbc_id );

    for(int ii=0; ii<numnode; ++ii)
    {
      const int node_index = infnbc -> get_LDN( nbc_id, ii );
      
      const int base_idx[3] = { node_index*4+1, node_index*4+2, node_index*4+3 };

      double base_vals[3];

      VecGetValues(sol_base->solution, 3, base_idx, base_vals);

      const double perturb_x = MATH_T::gen_double_rand_normal(0, std_dev);
      const double perturb_y = MATH_T::gen_double_rand_normal(0, std_dev);
      const double perturb_z = MATH_T::gen_double_rand_normal(0, std_dev);

      const double vals[3] = { base_vals[0] * factor * (1.0 + perturb_x), 
        base_vals[1] * factor * (1.0 + perturb_y),
        base_vals[2] * factor * (1.0 + perturb_z) };

      const int velo_idx[3] = { node_index*3, node_index*3+1, node_index*3+2 };

      VecSetValues(velo->solution, 3, velo_idx, vals, INSERT_VALUES);
    }
  }

  velo->Assembly_GhostUpdate();
}

void PTime_NS_HERK_Solver::rescale_inflow_dot_velo( const double &stime,
    PDNSolution * const &dot_velo ) const
{
  const int num_nbc = infnbc -> get_num_nbc();

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    const int numnode = infnbc -> get_Num_LD( nbc_id );

    const double factor  = dot_flrate -> get_flow_rate( nbc_id, stime );
    const double std_dev = dot_flrate -> get_flow_TI_std_dev( nbc_id );

    for(int ii=0; ii<numnode; ++ii)
    {
      const int node_index = infnbc -> get_LDN( nbc_id, ii );
      
      const int base_idx[3] = { node_index*4+1, node_index*4+2, node_index*4+3 };

      double base_vals[3];

      VecGetValues(sol_base->solution, 3, base_idx, base_vals);

      const double perturb_x = MATH_T::gen_double_rand_normal(0, std_dev);
      const double perturb_y = MATH_T::gen_double_rand_normal(0, std_dev);
      const double perturb_z = MATH_T::gen_double_rand_normal(0, std_dev);

      const double vals[3] = { base_vals[0] * factor * (1.0 + perturb_x), 
        base_vals[1] * factor * (1.0 + perturb_y),
        base_vals[2] * factor * (1.0 + perturb_z) };

      const int dot_velo_idx[3] = { node_index*3, node_index*3+1, node_index*3+2 };

      VecSetValues(dot_velo->solution, 3, dot_velo_idx, vals, INSERT_VALUES);
    }
  }

  dot_velo->Assembly_GhostUpdate();
}

void PTime_NS_HERK_Solver::Update_pressure_velocity(     
    PDNSolution * const &velo,
    PDNSolution * const &pres,
    const PDNSolution * const &step) const
  {
    Vec lvelo, lpres, lstep;
    double * array_velo, * array_pres, * array_step;

    VecGhostGetLocalForm(velo->solution, &lvelo);    
    VecGhostGetLocalForm(pres->solution, &lpres);    
    VecGhostGetLocalForm(step->solution, &lstep);

    VecGetArray(lvelo, &array_velo);
    VecGetArray(lpres, &array_pres); 
    VecGetArray(lstep, &array_step);

    // for(int ii=0; ii<nlocalnode; ++ii)
    // {
    //   array_pres[ii       ] -= array_step[ii*4 + 0];
    //   array_velo[ii*3 + 0 ] -= array_step[ii*4 + 1];
    //   array_velo[ii*3 + 1 ] -= array_step[ii*4 + 2];
    //   array_velo[ii*3 + 2 ] -= array_step[ii*4 + 3];
    // }

    for(int ii=0; ii<nlocalnode; ++ii)
    {
      array_velo[ii*3 + 0 ] -= array_step[ii*3 + 0];
      array_velo[ii*3 + 1 ] -= array_step[ii*3 + 1];
      array_velo[ii*3 + 2 ] -= array_step[ii*3 + 2];
      array_pres[ii       ] -= array_step[nlocalnode*3 + ii];
    }

    VecRestoreArray(lvelo, &array_velo);    
    VecRestoreArray(lpres, &array_pres);

    VecGhostRestoreLocalForm(velo->solution, &lvelo);
    VecGhostRestoreLocalForm(pres->solution, &lpres);
    
    velo->GhostUpdate();
    pres->GhostUpdate();  
  }

void PTime_NS_HERK_Solver::Update_solutions(     
    const PDNSolution * const &velo,
    const PDNSolution * const &pres,
    PDNSolution * const &sol) const
  {
    Vec lvelo, lpres, lsol;
    double * array_velo, * array_pres, * array_sol;

    VecGhostGetLocalForm(velo->solution, &lvelo);    
    VecGhostGetLocalForm(pres->solution, &lpres);    
    VecGhostGetLocalForm(sol->solution, &lsol);

    VecGetArray(lvelo, &array_velo);
    VecGetArray(lpres, &array_pres); 
    VecGetArray(lsol, &array_sol);

    for(int ii=0; ii<nlocalnode; ++ii)
    {
      array_sol[ii*4 + 0] = array_pres[ii       ];
      array_sol[ii*4 + 1] = array_velo[ii*3 + 0 ];
      array_sol[ii*4 + 2] = array_velo[ii*3 + 1 ];
      array_sol[ii*4 + 3] = array_velo[ii*3 + 2 ];
    }

    VecRestoreArray(lsol, &array_sol);    

    VecGhostRestoreLocalForm(sol->solution, &lsol);
    
    sol->GhostUpdate();
  }

// EOF
