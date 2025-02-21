#ifndef PTIME_NS_SOLVER_HPP
#define PTIME_NS_SOLVER_HPP
// ==================================================================
// PTime_NS_Solver.hpp
//
// Parallel time solver for NS.
//
// Author: Ju Liu
// Date: May 23 2017
// ==================================================================
#include "PDNTimeStep.hpp"
#include "PNonlinear_NS_Solver.hpp"

class PTime_NS_Solver
{
  public:
    PTime_NS_Solver( 
        std::unique_ptr<PNonlinear_NS_Solver> in_nsolver,
        const std::string &input_name,      
        const int &input_record_freq, const int &input_renew_tang_freq, 
        const double &input_final_time );

    ~PTime_NS_Solver() = default;

    void print_info() const;

    void print_lsolver_info() const {nsolver -> print_lsolver_info();}

    // ------------------------------------------------------------------------
    // Generate a file name for inlet/outlet face as prefix_xxx_data.txt
    // ------------------------------------------------------------------------
    std::string gen_flowfile_name(const std::string &prefix, const int &id) const
    {
      std::ostringstream ss;
      ss << prefix;

      if(id < 10) ss << "00";
      else if(id < 100) ss << "0";

      ss << id << "_data.txt";
      return ss.str();
    }

    void TM_NS_GenAlpha(
        const bool &restart_init_assembly_flag,
        std::unique_ptr<PDNSolution> init_dot_sol,
        std::unique_ptr<PDNSolution> init_sol,
        std::unique_ptr<PDNTimeStep> time_info,
        const ALocal_InflowBC * const &infnbc_part,
        IGenBC * const &gbc, 
        IPGAssem * const &gassem_ptr ) const;

  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const int renew_tang_freq; // the frequency for renewing tangents
    const std::string pb_name; // the problem base name for the solution

    const std::unique_ptr<PNonlinear_NS_Solver> nsolver;

    std::string Name_Generator( const int &counter ) const;

    std::string Name_dot_Generator( const int &counter ) const;

    void Write_restart_file(const PDNTimeStep * const &timeinfo,
        const std::string &solname ) const;
};

#endif
