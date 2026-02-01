#ifndef PTIME_SOLID_SOLVER_HPP
#define PTIME_SOLID_SOLVER_HPP
// ============================================================================
// PTime_Solid_Solver.hpp
//
// Parallel time solver for hyperelastic solid problems.
//
// Date: Jan 29 2026
// ==========================================================================
#include "PDNTimeStep.hpp"
#include "PNonlinear_Solid_Solver.hpp"

class PTime_Solid_Solver
{
  public:
    PTime_Solid_Solver(
        std::unique_ptr<PNonlinear_Solid_Solver> in_nsolver,
        const std::string &input_name,
        const int &input_record_freq,
        const int &input_renew_tang_freq,
        const double &input_final_time );

    ~PTime_Solid_Solver() = default;

    void print_info() const;

    void print_lsolver_info() const {nsolver->print_lsolver_info();}

    void TM_Solid_GenAlpha(
        const bool &restart_init_assembly_flag,
        const IS &is_v,
        const IS &is_p,
        std::unique_ptr<PDNSolution> init_dot_disp,
        std::unique_ptr<PDNSolution> init_dot_velo,
        std::unique_ptr<PDNSolution> init_dot_pres,
        std::unique_ptr<PDNSolution> init_disp,
        std::unique_ptr<PDNSolution> init_velo,
        std::unique_ptr<PDNSolution> init_pres,
        std::unique_ptr<PDNTimeStep> time_info ) const;

  private:
    const double final_time;
    const int sol_record_freq;
    const int renew_tang_freq;
    const std::string pb_name;

    const std::unique_ptr<PNonlinear_Solid_Solver> nsolver;

    std::string Name_Generator( const std::string &middle_name,
        const int &counter ) const;

    std::string Name_dot_Generator( const std::string &middle_name,
        const int &counter ) const;

    PTime_Solid_Solver() = delete;
};

#endif
