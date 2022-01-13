#ifndef PTIME_FSI_SOLVER_HPP
#define PTIME_FSI_SOLVER_HPP
// ============================================================================
// PTime_FSI_Solver.hpp
//
// Parallel time solver for time integration.
//
// Date: Jan. 13 2022
// Author: Ju Liu
// ============================================================================
#include "PDNTimeStep.hpp"
#include "PNonlinear_FSI_Solver.hpp"

class PTime_FSI_Solver
{
  public:
    PTime_FSI_Solver( const std::string &input_name, const int &input_record_freq, 
        const int &input_renew_tang_freq, const double &input_final_time );

    ~PTime_FSI_Solver();

    void print_info() const;

  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const int renew_tang_freq; // the frequency for renewing tangents
    const std::string pb_name; // the problem base name for the solution

    std::string Name_Generator( const int &counter ) const;

    std::string Name_dot_Generator( const int &counter ) const;

    void Write_restart_file(const PDNTimeStep * const &timeinfo,
        const std::string &solname ) const;
};

#endif
