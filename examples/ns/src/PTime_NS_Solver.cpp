#include "PTime_NS_Solver.hpp"

PTime_NS_Solver::PTime_NS_Solver(
    const std::string &input_name, const int &input_record_freq,
    const int &input_renew_tang_freq, const double &input_final_time )
: final_time(input_final_time), sol_record_freq(input_record_freq),
  renew_tang_freq(input_renew_tang_freq), pb_name(input_name)
{}

std::string PTime_NS_Solver::Name_Generator(const int &counter) const
{
  std::ostringstream temp;
  temp.str("");
  temp<<900000000 + counter;

  std::string out_name(pb_name);
  out_name.append(temp.str());
  return out_name;
}

std::string PTime_NS_Solver::Name_dot_Generator(const int &counter) const
{
  std::ostringstream temp;
  temp.str("");
  temp<<900000000 + counter;

  std::string out_name("dot_");
  out_name.append(pb_name);
  out_name.append(temp.str());
  return out_name;
}

void PTime_NS_Solver::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("Time stepping solver setted up:\n");
  SYS_T::commPrint("  final time: %e \n", final_time);
  SYS_T::commPrint("  solution record frequency : %d \n", sol_record_freq);
  SYS_T::commPrint("  tangent update frequency over time steps: %d \n", renew_tang_freq);
  SYS_T::commPrint("  solution base name: %s \n", pb_name.c_str());
  SYS_T::commPrint("----------------------------------------------------------- \n");
}

void PTime_NS_Solver::Write_restart_file(const PDNTimeStep * const &timeinfo,
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

void PTime_NS_Solver::TM_NS_GenAlpha( 
    const bool &restart_init_assembly_flag,
    PDNSolution * const &sol_base,
    const PDNSolution * const &init_dot_sol,
    const PDNSolution * const &init_sol,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    PDNTimeStep * const &time_info,
    const ICVFlowRate * const flr_ptr,
    const APart_Node * const &pNode_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_InflowBC * const &infnbc_part,
    const ALocal_EBC * const &ebc_part,
    IGenBC * const &gbc,
    const ALocal_WeakBC * const &wbc_part,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    FEAElement * const &elementvs,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_fluid_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_NS_Solver * const &nsolver_ptr ) const
{
  PDNSolution * pre_sol = new PDNSolution(*init_sol);
  PDNSolution * cur_sol = new PDNSolution(*init_sol);
  PDNSolution * pre_dot_sol = new PDNSolution(*init_dot_sol);
  PDNSolution * cur_dot_sol = new PDNSolution(*init_dot_sol);

  // If this is a restart run, do not re-write the solution binaries
  if(restart_init_assembly_flag == false)
  {
    const auto sol_name = Name_Generator(time_info->get_index());
    cur_sol->WriteBinary(sol_name);
    
    const auto sol_dot_name = Name_dot_Generator(time_info->get_index());
    cur_dot_sol->WriteBinary(sol_dot_name);
  }

  bool conv_flag, renew_flag;
  int nl_counter = 0;

  bool rest_flag = restart_init_assembly_flag;

  SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
      time_info->get_time(), time_info->get_step(), time_info->get_index(),
      SYS_T::get_time().c_str());

  // Enter into time integration
  while( time_info->get_time() < final_time )
  {
    if(time_info->get_index() % renew_tang_freq == 0 || rest_flag )
    {
      renew_flag = true;
      rest_flag = false;
    }
    else renew_flag = false;

    // If the previous step is solved in ONE Newton iteration, we do not update
    // the tangent matrix
    if( nl_counter == 1 ) renew_flag = false;

    // Call the nonlinear equation solver
    nsolver_ptr->GenAlpha_Solve_NS( renew_flag, 
        time_info->get_time(), time_info->get_step(),
        sol_base, pre_dot_sol, pre_sol, tmga_ptr, flr_ptr,
        alelem_ptr, lien_ptr, feanode_ptr, nbc_part, infnbc_part,
        ebc_part, gbc, wbc_part, bc_mat, elementv, elements, elementvs, quad_v, quad_s, lassem_fluid_ptr,
        gassem_ptr, lsolver_ptr, cur_dot_sol, cur_sol, conv_flag, nl_counter );

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

      const auto sol_dot_name = Name_dot_Generator(time_info->get_index());
      cur_dot_sol->WriteBinary(sol_dot_name);
    }

    // Calculate the flow rate & averaged pressure on all outlets
    for(int face=0; face<ebc_part -> get_num_ebc(); ++face)
    {
      // Calculate the 3D dot flow rate on the outlet
      const double dot_face_flrate = gassem_ptr -> Assem_surface_flowrate( 
          cur_dot_sol, lassem_fluid_ptr, elements, quad_s, ebc_part, face); 

      // Calculate the 3D flow rate on the outlet
      const double face_flrate = gassem_ptr -> Assem_surface_flowrate( 
          cur_sol, lassem_fluid_ptr, elements, quad_s, ebc_part, face); 

      // Calculate the 3D averaged pressure on the outlet
      const double face_avepre = gassem_ptr -> Assem_surface_ave_pressure( 
          cur_sol, lassem_fluid_ptr, elements, quad_s, ebc_part, face);

      // Calculate the 0D pressure from LPN model
      const double dot_lpn_flowrate = dot_face_flrate;
      const double lpn_flowrate = face_flrate;
      const double lpn_pressure = gbc -> get_P( face, dot_lpn_flowrate, lpn_flowrate,
         time_info -> get_time() );

      // Update the initial values in genbc
      gbc -> reset_initial_sol( face, lpn_flowrate, lpn_pressure, time_info->get_time(), false );

      // On the CPU 0, write the time, flow rate, averaged pressure, and 0D
      // calculated pressure into the txt file, which is first generated in the
      // driver
      if( SYS_T::get_MPI_rank() == 0 )
      {
        std::ofstream ofile;
        ofile.open( ebc_part->gen_flowfile_name(face).c_str(), std::ofstream::out | std::ofstream::app );
        ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'<<dot_face_flrate<<'\t'<<face_flrate<<'\t'<<face_avepre<<'\t'<<lpn_pressure<<'\n';
        ofile.close();
      }
      MPI_Barrier(PETSC_COMM_WORLD);
    }
   
    // Calcualte the inlet data
    for(int face=0; face<infnbc_part -> get_num_nbc(); ++face)
    {
      const double inlet_face_flrate = gassem_ptr -> Assem_surface_flowrate(
          cur_sol, lassem_fluid_ptr, elements, quad_s, infnbc_part, face ); 

      const double inlet_face_avepre = gassem_ptr -> Assem_surface_ave_pressure(
          cur_sol, lassem_fluid_ptr, elements, quad_s, infnbc_part, face );

      if( SYS_T::get_MPI_rank() == 0 )
      {
        std::ofstream ofile;
        ofile.open( infnbc_part->gen_flowfile_name(face).c_str(), std::ofstream::out | std::ofstream::app );
        ofile<<time_info->get_index()<<'\t'<<time_info->get_time()<<'\t'<<inlet_face_flrate<<'\t'<<inlet_face_avepre<<'\n';
        ofile.close();
      } 
      MPI_Barrier(PETSC_COMM_WORLD);
    }

    // Prepare for next time step
    pre_sol->Copy(*cur_sol);
    pre_dot_sol->Copy(*cur_dot_sol);
  }

  delete pre_sol; delete cur_sol; delete pre_dot_sol; delete cur_dot_sol;
}

void PTime_NS_Solver::TM_NS_HERK( 
    const bool &restart_init_assembly_flag,
    const PDNSolution * const &sol_base,
    const PDNSolution * const &dot_sol_base,
    const PDNSolution * const &init_sol,
    const PDNSolution * const &init_velo,
    const PDNSolution * const &init_dot_velo,
    const PDNSolution * const &init_pres,
    const Runge_Kutta_Butcher * const &tm_RK_ptr,
    PDNTimeStep * const &time_info,
    const ICVFlowRate * const flr_ptr,
    const ICVFlowRate * const dot_flr_ptr,
    const APart_Node * const &pNode_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NBC * const &nbc_part,
    const ALocal_InflowBC * const &infnbc_part,
    const ALocal_EBC * const &ebc_part,
    IGenBC * const &gbc,
    const ALocal_WeakBC * const &wbc_part,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    FEAElement * const &elementvs,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_fluid_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_NS_Solver * const &nsolver_ptr ) const
{
  const int ss = tm_RK_ptr->get_RK_step();

  // n+1 步的子步速度解
  PDNSolution** cur_velo_sols = new PDNSolution*[ss];
  for(int ii = 0; ii < ss; ++ii)
    cur_velo_sols[ii] = new PDNSolution(*init_velo);

  // n+1 步的终步速度解
  PDNSolution * cur_velo = new PDNSolution(*init_velo);

  // n+1 步的最终步dot速度解
  PDNSolution * cur_dot_velo = new PDNSolution(*init_dot_velo);

  // n+1 步的子步压强解
  PDNSolution** cur_pres_sols = new PDNSolution*[ss];
  for(int ii = 0; ii < ss; ++ii)
    cur_pres_sols[ii] = new PDNSolution(*init_pres);

  // n+1 步的终步压强解
  PDNSolution * cur_pres = new PDNSolution(*init_pres);

  // n+1 步的解
  PDNSolution * cur_sol = new PDNSolution(*init_sol);

  // n 步的子步速度解
  PDNSolution** pre_velo_sols = new PDNSolution*[ss];
  for(int ii = 0; ii < ss; ++ii)
    pre_velo_sols[ii] = new PDNSolution(*init_velo);

  // n 步的子步压强解
  PDNSolution** pre_pres_sols = new PDNSolution*[ss];
  for(int ii = 0; ii < ss; ++ii)
    pre_pres_sols[ii] = new PDNSolution(*init_pres);    

  // n 步的终步速度解
  PDNSolution * pre_velo = new PDNSolution(*init_velo);

  // n 步的终步压强解
  PDNSolution * pre_pres = new PDNSolution(*init_pres);

  // n-1 步的终步速度解
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
    // Call the nonlinear equation solver
    nsolver_ptr->HERK_Solve_NS(
        time_info->get_time(), time_info->get_step(),
        sol_base, dot_sol_base, cur_velo_sols, cur_velo, cur_dot_velo, 
        cur_pres_sols, cur_pres, pre_velo_sols, pre_velo,
        pre_pres_sols, pre_pres, pre_velo_before, tm_RK_ptr, flr_ptr, dot_flr_ptr,
        alelem_ptr, lien_ptr, pNode_ptr, feanode_ptr, nbc_part, infnbc_part,
        ebc_part, gbc, wbc_part, bc_mat, elementv, elements, elementvs, quad_v, quad_s, lassem_fluid_ptr,
        gassem_ptr, lsolver_ptr, cur_sol);

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

    // for(int ii = 0; ii < ss; ++ii)
    //   VecView(cur_velo_sols[ii]->solution, PETSC_VIEWER_STDOUT_WORLD);

    // // n+1 步的子步速度解
    for(int ii = 0; ii < ss; ++ii)
      cur_velo_sols[ii]->Copy(*init_velo);

    // n+1 步的终步速度解
    cur_velo->Copy(*init_velo);

    // n+1 步的最终步dot速度解
    cur_dot_velo->Copy(*init_dot_velo);

    // n+1 步的子步压强解
    for(int ii = 0; ii < ss; ++ii)
      cur_pres_sols[ii]->Copy(*init_pres);

    // n+1 步的终步压强解
    cur_pres->Copy(*init_pres);

    // n+1 步的解
    cur_sol->Copy(*init_sol);    
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

// EOF
