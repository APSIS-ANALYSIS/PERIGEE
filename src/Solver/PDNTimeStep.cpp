#include "PDNTimeStep.hpp"

PDNTimeStep::PDNTimeStep( const int &input_index, const double &input_time,
    const double &input_step )
{
  time = input_time;
  time_index = input_index;
  time_step = input_step;

  time_history.clear();
  time_step_history.clear();
  index_history.clear();

  time_history.push_back(time);
  time_step_history.push_back(time_step);
  index_history.push_back(time_index);
}


PDNTimeStep::~PDNTimeStep()
{}

void PDNTimeStep::TimeIncrement()
{
  time = time + time_step;
  time_index += 1;
  
  time_history.push_back(time);
  time_step_history.push_back(time_step);
  index_history.push_back(time_index);
}

void PDNTimeStep::TimeIncrement(const double &input_time_step)
{
  time_step = input_time_step;
  time = time + time_step;
  time_index += 1;
  
  time_history.push_back(time);
  time_step_history.push_back(time_step);
  index_history.push_back(time_index);
}


void PDNTimeStep::WriteTimeInfo() const
{
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(rank == 0)
  {
    std::ofstream timefile;
    timefile.open("Time_log.txt", std::ios::app);
    if(timefile.is_open())
    {
      for(unsigned int ii=0; ii<index_history.size(); ++ii)
      {
        timefile<<index_history[ii];
        timefile<<"  "<<time_step_history[ii];
        timefile<<"  "<<time_history[ii]<<std::endl;
      }
      timefile.close();
    }
    else
      std::cerr<<"Unable to open/create Time_log.txt file. \n";
  }
}


void PDNTimeStep::WriteTimeInfo_step() const
{
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(rank == 0)
  {
    std::ofstream timefile;
    timefile.open("Time_log.txt", std::ios::app);
    if(timefile.is_open())
    {
      timefile<<time_index;
      timefile<<"  "<<time_step;
      timefile<<"  "<<time<<std::endl;
      timefile.close();
    }
    else
      std::cerr<<"Unable to open/create Time_log.txt file. \n";
  }
}


// EOF
