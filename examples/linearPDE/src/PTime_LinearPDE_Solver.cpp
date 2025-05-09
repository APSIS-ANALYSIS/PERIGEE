#include "PTime_LinearPDE_Solver.hpp"

PTime_LinearPDE_Solver::PTime_LinearPDE_Solver( 
    std::unique_ptr<PNonlinear_LinearPDE_Solver> in_nsolver,
    const std::string &input_name,
    const int &input_record_freq, const int &input_renew_tang_freq,
    const double &input_final_time )
: final_time(input_final_time), sol_record_freq(input_record_freq),
  renew_tang_freq(input_renew_tang_freq), pb_name(input_name),
  nsolver(std::move(in_nsolver))
{}

void PTime_LinearPDE_Solver::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("Time stepping solver setted up:\n");
  SYS_T::commPrint("  final time: %e \n", final_time);
  SYS_T::commPrint("  solution record frequency : %d \n", sol_record_freq);
  SYS_T::commPrint("  tangent update frequency over time steps: %d \n", renew_tang_freq);
  SYS_T::commPrint("  solution base name: %s \n", pb_name.c_str());
  SYS_T::commPrint("----------------------------------------------------------- \n");
}

std::string PTime_LinearPDE_Solver::Name_Generator( const std::string &middle_name,
    const int &counter ) const
{
  std::ostringstream temp;
  temp.str("");
  temp<<900000000 + counter;

  std::string out_name(pb_name);
  out_name.append(middle_name);
  out_name.append(temp.str());
  return out_name;
}

std::string PTime_LinearPDE_Solver::Name_dot_Generator( const std::string &middle_name,
    const int &counter ) const
{
  std::ostringstream temp;
  temp.str("");
  temp<<900000000 + counter;

  std::string out_name("dot_");
  out_name.append(pb_name);
  out_name.append(middle_name);
  out_name.append(temp.str());
  return out_name;
}

void PTime_LinearPDE_Solver::TM_GenAlpha_Transport(
    const bool &restart_init_assembly_flag,
    std::unique_ptr<PDNSolution> init_dot_sol,
    std::unique_ptr<PDNSolution> init_sol,
    std::unique_ptr<PDNTimeStep> time_info ) const
{
  auto pre_sol     = SYS_T::make_unique<PDNSolution>(*init_sol);
  auto cur_sol     = SYS_T::make_unique<PDNSolution>(*init_sol);
  auto pre_dot_sol = SYS_T::make_unique<PDNSolution>(*init_dot_sol);
  auto cur_dot_sol = SYS_T::make_unique<PDNSolution>(*init_dot_sol);

  // If this is a restart run, do not re-write the solution binaries
  if(restart_init_assembly_flag == false)
  {
    const auto sol_name = Name_Generator("temp_", time_info->get_index());
    cur_sol->WriteBinary(sol_name.c_str());

    const auto sol_dot_name = Name_dot_Generator("temp_", time_info->get_index());
    cur_dot_sol->WriteBinary(sol_dot_name.c_str());
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
    nsolver->GenAlpha_Solve_Transport( renew_flag,
        time_info->get_time(), time_info->get_step(),
        pre_dot_sol.get(), pre_sol.get(),
        cur_dot_sol.get(), cur_sol.get(), 
        conv_flag, nl_counter );

    // Update the time step information
    time_info->TimeIncrement();

    SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index(),
        SYS_T::get_time().c_str());

    // Record solution if meets criteria
    if( time_info->get_index()%sol_record_freq == 0 )
    {
      const auto sol_name = Name_Generator("temp_", time_info->get_index() );
      cur_sol->WriteBinary(sol_name.c_str());

      const auto sol_dot_name = Name_dot_Generator("temp_", time_info->get_index());
      cur_dot_sol->WriteBinary(sol_dot_name.c_str());
    }

    // Prepare for next time step
    pre_sol->Copy(*cur_sol);
    pre_dot_sol->Copy(*cur_dot_sol);
  } 
}

void PTime_LinearPDE_Solver::TM_GenAlpha_Elastodynamics(
    const bool &restart_init_assembly_flag,
    std::unique_ptr<PDNSolution> init_dot_disp,
    std::unique_ptr<PDNSolution> init_dot_velo,
    std::unique_ptr<PDNSolution> init_disp,
    std::unique_ptr<PDNSolution> init_velo,
    std::unique_ptr<PDNTimeStep> time_info ) const
{
  auto pre_disp = SYS_T::make_unique<PDNSolution>( *init_disp );
  auto cur_disp = SYS_T::make_unique<PDNSolution>( *init_disp );
  auto pre_velo = SYS_T::make_unique<PDNSolution>( *init_velo );
  auto cur_velo = SYS_T::make_unique<PDNSolution>( *init_velo );
  auto pre_dot_disp = SYS_T::make_unique<PDNSolution>( *init_dot_disp );
  auto cur_dot_disp = SYS_T::make_unique<PDNSolution>( *init_dot_disp );
  auto pre_dot_velo = SYS_T::make_unique<PDNSolution>( *init_dot_velo );
  auto cur_dot_velo = SYS_T::make_unique<PDNSolution>( *init_dot_velo );

  // If this is a restart run, do not re-write the solution binaries
  if(restart_init_assembly_flag == false)
  {
    auto sol_name = Name_Generator("disp_", time_info->get_index());
    cur_disp->WriteBinary(sol_name.c_str());

    sol_name = Name_Generator("velo_", time_info->get_index());
    cur_velo->WriteBinary(sol_name.c_str());

    auto sol_dot_name = Name_dot_Generator("disp_", time_info->get_index());
    cur_dot_disp->WriteBinary(sol_dot_name.c_str());

    sol_dot_name = Name_dot_Generator("velo_", time_info->get_index());
    cur_dot_velo->WriteBinary(sol_dot_name.c_str());
  }

  bool conv_flag;
  bool renew_flag = true;
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
    nsolver->GenAlpha_Solve_Elastodynamics( renew_flag,
        time_info->get_time(), time_info->get_step(),
        pre_dot_disp.get(), pre_dot_velo.get(), pre_disp.get(), pre_velo.get(),
        cur_dot_disp.get(), cur_dot_velo.get(), cur_disp.get(), cur_velo.get(),
        conv_flag, nl_counter );

    // Update the time step information
    time_info->TimeIncrement();

    SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index(),
        SYS_T::get_time().c_str());

    // Record solution if meets criteria
    if( time_info->get_index()%sol_record_freq == 0 )
    {
      auto sol_name = Name_Generator("disp_", time_info->get_index());
      cur_disp->WriteBinary(sol_name.c_str());

      sol_name = Name_Generator("velo_", time_info->get_index());
      cur_velo->WriteBinary(sol_name.c_str());

      auto sol_dot_name = Name_dot_Generator("disp_", time_info->get_index());
      cur_dot_disp->WriteBinary(sol_dot_name.c_str());

      sol_dot_name = Name_dot_Generator("velo_", time_info->get_index());
      cur_dot_velo->WriteBinary(sol_dot_name.c_str());
    }

    // Prepare for next time step
    pre_disp -> Copy( cur_disp.get() );
    pre_velo -> Copy( cur_velo.get() );
    pre_dot_disp -> Copy( cur_dot_disp.get() );
    pre_dot_velo -> Copy( cur_dot_velo.get() );
  }
}

// EOF
