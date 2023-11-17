#ifndef PDNTIMESTEP_HPP
#define PDNTIMESTEP_HPP
// ==================================================================
// PDNTimeStep.hpp
//
// This class defines (1) time, (2) time step, (3) time index. 
// We can perform time increment, time step update, time index record
// within this class.
//
// Date: Dec 9 2013
// ==================================================================
#include <vector>
#include "Sys_Tools.hpp"

class PDNTimeStep
{
  public:
    PDNTimeStep( const int &input_index, const double &input_time,
        const double &input_step );

    ~PDNTimeStep() = default;

    // Functions that give access to the class member data
    double get_time() const {return time;}

    double get_step() const {return time_step;}
    
    int get_index() const {return time_index;}

    // Perform time increment with the default time step
    void TimeIncrement();

    // Perform time increment with the given time step
    void TimeIncrement(const double &input_time_step);

    // Change the time step size
    void UpdateTimeStep(const double &new_time_step)
    {time_step = new_time_step;}

    // ! Write the three _history vectors to file on proc 0
    void WriteTimeInfo() const;

    // ! Append the time-time_step-time_index on Time_log.txt
    void WriteTimeInfo_step() const;

  private:
    int time_index;
    double time, time_step;

    std::vector<int> index_history;
    std::vector<double> time_history, time_step_history;
};

#endif
