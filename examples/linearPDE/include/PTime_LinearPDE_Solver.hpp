#ifndef PTIME_LINEARPDE_SOLVER_HPP
#define PTIME_LINEARPDE_SOLVER_HPP
// ==================================================================
// PTime_LinearPDE_Solver.hpp
//
// Parallel time solver for linear PDE.
//
// Date: Oct. 26 2023
// ==================================================================
#include "PDNTimeStep.hpp"
#include "PNonlinear_LinearPDE_Solver.hpp"

class PTime_LinearPDE_Solver
{
  public:
    PTime_LinearPDE_Solver( 
        std::unique_ptr<PNonlinear_LinearPDE_Solver> in_nsolver,
        const std::string &input_name,
        const int &input_record_freq, const int &input_renew_tang_freq,
        const double &input_final_time );

    ~PTime_LinearPDE_Solver() = default;

    void print_info() const;

    void print_lsolver_info() const {nsolver -> print_lsolver_info();}

    void TM_GenAlpha_Transport(
        const bool &restart_init_assembly_flag,
        std::unique_ptr<PDNSolution> init_dot_sol,
        std::unique_ptr<PDNSolution> init_sol,
        std::unique_ptr<PDNTimeStep> time_info ) const;
   
    void TM_GenAlpha_Elastodynamics(
        const bool &restart_init_assembly_flag,
        std::unique_ptr<PDNSolution> init_dot_disp,
        std::unique_ptr<PDNSolution> init_dot_velo,
        std::unique_ptr<PDNSolution> init_disp,
        std::unique_ptr<PDNSolution> init_velo,
        std::unique_ptr<PDNTimeStep> time_info ) const;

  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const int renew_tang_freq; // the frequency for renewing tangents
    const std::string pb_name; // the problem base name for the solution

    const std::unique_ptr<PNonlinear_LinearPDE_Solver> nsolver;

    std::string Name_Generator( const std::string &middle_name,
        const int &counter ) const;

    std::string Name_dot_Generator( const std::string &middle_name,
        const int &counter ) const;
};

#endif
