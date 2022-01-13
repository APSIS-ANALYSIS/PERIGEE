#include "PTime_FSI_Solver.hpp"

PTime_FSI_Solver::PTime_FSI_Solver( const std::string &input_name, 
    const int &input_record_freq, const int &input_renew_tang_freq, 
    const double &input_final_time )
: final_time(input_final_time), sol_record_freq(input_record_freq),
  renew_tang_freq(input_renew_tang_freq), pb_name(input_name)
{}


PTime_FSI_Solver::~PTime_FSI_Solver()
{}


void PTime_FSI_Solver::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint( "final time: %e \n", final_time);
  SYS_T::commPrint( "solution record frequency : %d \n", sol_record_freq);
  SYS_T::commPrint( "tangent update frequency over time steps: %d \n", renew_tang_freq);
  SYS_T::commPrint( "solution base name: %s \n", pb_name.c_str());
  SYS_T::print_sep_line();
}


std::string PTime_FSI_Solver::Name_Generator( const int &counter ) const
{
  const int aux = 900000000 + counter;
  std::ostringstream temp;
  temp<<aux;

  std::string out_name(pb_name);
  out_name.append(temp.str());
  return out_name;
}


std::string PTime_FSI_Solver::Name_dot_Generator( const int &counter ) const
{
  const int aux = 900000000 + counter;
  std::ostringstream temp;
  temp<<aux;

  std::string out_name("dot_");
  out_name.append(pb_name);
  out_name.append(temp.str());
  return out_name;
}


void PTime_FSI_Solver::Write_restart_file(const PDNTimeStep * const &timeinfo,
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

// EOF
