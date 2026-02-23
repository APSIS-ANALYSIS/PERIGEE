#include "PTime_Solid_Solver.hpp"

PTime_Solid_Solver::PTime_Solid_Solver(
    std::unique_ptr<PNonlinear_Solid_Solver> in_nsolver,
    const std::string &input_name,
    const int &input_record_freq,
    const int &input_renew_tang_freq,
    const double &input_final_time )
: final_time(input_final_time), sol_record_freq(input_record_freq),
  renew_tang_freq(input_renew_tang_freq), pb_name(input_name),
  nsolver(std::move(in_nsolver))
{}

void PTime_Solid_Solver::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("Time stepping solver setted up:\n");
  SYS_T::commPrint("  final time: %e \n", final_time);
  SYS_T::commPrint("  solution record frequency : %d \n", sol_record_freq);
  SYS_T::commPrint("  tangent update frequency over time steps: %d \n", renew_tang_freq);
  SYS_T::commPrint("  solution base name: %s \n", pb_name.c_str());
  SYS_T::commPrint("----------------------------------------------------------- \n");
}

std::string PTime_Solid_Solver::Name_Generator( const std::string &middle_name,
    const int &counter ) const
{
  std::ostringstream temp;
  temp.str("");
  temp << 900000000 + counter;
  return pb_name + middle_name + temp.str();
}

std::string PTime_Solid_Solver::Name_dot_Generator( const std::string &middle_name,
    const int &counter ) const
{
  std::ostringstream temp;
  temp.str("");
  temp << 900000000 + counter;
  return std::string("dot_") + pb_name + middle_name + temp.str();
}

void PTime_Solid_Solver::TM_Solid_GenAlpha(
    const bool &restart_init_assembly_flag,
    const IS &is_v,
    const IS &is_p,
    std::unique_ptr<PDNSolution> init_dot_disp,
    std::unique_ptr<PDNSolution> init_dot_velo,
    std::unique_ptr<PDNSolution> init_dot_pres,
    std::unique_ptr<PDNSolution> init_disp,
    std::unique_ptr<PDNSolution> init_velo,
    std::unique_ptr<PDNSolution> init_pres,
    std::unique_ptr<PDNTimeStep> time_info ) const
{
  auto pre_dot_disp = SYS_T::make_unique<PDNSolution>(*init_dot_disp);
  auto pre_dot_velo = SYS_T::make_unique<PDNSolution>(*init_dot_velo);
  auto pre_dot_pres = SYS_T::make_unique<PDNSolution>(*init_dot_pres);

  auto pre_disp = SYS_T::make_unique<PDNSolution>(*init_disp);
  auto pre_velo = SYS_T::make_unique<PDNSolution>(*init_velo);
  auto pre_pres = SYS_T::make_unique<PDNSolution>(*init_pres);

  auto cur_dot_disp = SYS_T::make_unique<PDNSolution>(*init_dot_disp);
  auto cur_dot_velo = SYS_T::make_unique<PDNSolution>(*init_dot_velo);
  auto cur_dot_pres = SYS_T::make_unique<PDNSolution>(*init_dot_pres);

  auto cur_disp = SYS_T::make_unique<PDNSolution>(*init_disp);
  auto cur_velo = SYS_T::make_unique<PDNSolution>(*init_velo);
  auto cur_pres = SYS_T::make_unique<PDNSolution>(*init_pres);

  if( restart_init_assembly_flag == false )
  {
    std::string sol_name = Name_Generator("disp_", time_info->get_index());
    cur_disp->WriteBinary(sol_name.c_str());

    sol_name = Name_Generator("velo_", time_info->get_index());
    cur_velo->WriteBinary(sol_name.c_str());

    sol_name = Name_Generator("pres_", time_info->get_index());
    cur_pres->WriteBinary(sol_name.c_str());

    std::string sol_dot_name = Name_dot_Generator("disp_", time_info->get_index());
    cur_dot_disp->WriteBinary(sol_dot_name.c_str());

    sol_dot_name = Name_dot_Generator("velo_", time_info->get_index());
    cur_dot_velo->WriteBinary(sol_dot_name.c_str());

    sol_dot_name = Name_dot_Generator("pres_", time_info->get_index());
    cur_dot_pres->WriteBinary(sol_dot_name.c_str());
  }

  bool conv_flag, renew_flag;
  int nl_counter = 0;

  bool rest_flag = restart_init_assembly_flag;

  SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
      time_info->get_time(), time_info->get_step(), time_info->get_index(),
      SYS_T::get_time().c_str());

  while( time_info->get_time() < final_time )
  {
    if( time_info->get_index() % renew_tang_freq == 0 || rest_flag )
    {
      renew_flag = true;
      rest_flag = false;
    }
    else renew_flag = false;

    if( nl_counter == 1 ) renew_flag = false;

    nsolver->GenAlpha_Seg_solve_Solid( renew_flag,
        time_info->get_time(), time_info->get_step(),
        is_v, is_p,
        pre_dot_disp.get(), pre_dot_velo.get(), pre_dot_pres.get(),
        pre_disp.get(), pre_velo.get(), pre_pres.get(),
        cur_dot_disp.get(), cur_dot_velo.get(), cur_dot_pres.get(),
        cur_disp.get(), cur_velo.get(), cur_pres.get(),
        conv_flag, nl_counter );

    time_info->TimeIncrement();

    SYS_T::commPrint("Time = %e, dt = %e, index = %d, %s \n",
        time_info->get_time(), time_info->get_step(), time_info->get_index(),
        SYS_T::get_time().c_str());

    if( time_info->get_index() % sol_record_freq == 0 )
    {
      std::string sol_name = Name_Generator("disp_", time_info->get_index());
      cur_disp->WriteBinary(sol_name.c_str());

      sol_name = Name_Generator("velo_", time_info->get_index());
      cur_velo->WriteBinary(sol_name.c_str());

      sol_name = Name_Generator("pres_", time_info->get_index());
      cur_pres->WriteBinary(sol_name.c_str());

      std::string sol_dot_name = Name_dot_Generator("disp_", time_info->get_index());
      cur_dot_disp->WriteBinary(sol_dot_name.c_str());

      sol_dot_name = Name_dot_Generator("velo_", time_info->get_index());
      cur_dot_velo->WriteBinary(sol_dot_name.c_str());

      sol_dot_name = Name_dot_Generator("pres_", time_info->get_index());
      cur_dot_pres->WriteBinary(sol_dot_name.c_str());
    }

    pre_dot_disp->Copy( cur_dot_disp.get() );
    pre_dot_velo->Copy( cur_dot_velo.get() );
    pre_dot_pres->Copy( cur_dot_pres.get() );

    pre_disp->Copy( cur_disp.get() );
    pre_velo->Copy( cur_velo.get() );
    pre_pres->Copy( cur_pres.get() );
  }
}

// EOF
