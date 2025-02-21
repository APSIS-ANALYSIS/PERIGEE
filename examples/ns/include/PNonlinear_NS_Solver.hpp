#ifndef PNONLINEAR_NS_SOLVER_HPP
#define PNONLINEAR_NS_SOLVER_HPP
// ==================================================================
// PNonlinear_NS_Solver.hpp
// 
// Parallel nonlinear solver for Navier-Stokes equations. 
//
// Author: Ju Liu
// Date: Feb 11 2020
// ==================================================================
#include "IFlowRate.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "IPGAssem.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"
#include "PDNSolution_NS.hpp"

class PNonlinear_NS_Solver
{
  public:
    PNonlinear_NS_Solver( 
        std::unique_ptr<PLinear_Solver_PETSc> in_lsolver,
        std::unique_ptr<Matrix_PETSc> in_bc_mat,
        std::unique_ptr<TimeMethod_GenAlpha> in_tmga,
        std::unique_ptr<IFlowRate> in_flrate,
        std::unique_ptr<PDNSolution> in_sol_base,
        const double &input_nrtol, const double &input_natol, 
        const double &input_ndtol, const int &input_max_iteration, 
        const int &input_renew_freq, 
        const int &input_renew_threshold = 4 );

    ~PNonlinear_NS_Solver() = default;

    int get_non_max_its() const {return nmaxits;}

    void print_info() const;

    void print_lsolver_info() const {lsolver->print_info();}

    // --------------------------------------------------------------
    // GenAlpha_Solve_NS:
    // This is a solver for fluid dynamics.
    //
    // This solver solves the Navier-Stokes using 2nd-order Generalized
    // alpha method.
    // --------------------------------------------------------------
    void GenAlpha_Solve_NS(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &pre_dot_sol,
        const PDNSolution * const &pre_sol,
        PDNSolution * const &dot_sol,
        PDNSolution * const &sol,
        const ALocal_InflowBC * const &infnbc_part,
        const IGenBC * const &gbc,
        IPGAssem * const &gassem_ptr,
        bool &conv_flag, int &nl_counter ) const;

  private:
    const double nr_tol, na_tol, nd_tol;
    const int nmaxits, nrenew_freq, nrenew_threshold;

    const std::unique_ptr<PLinear_Solver_PETSc> lsolver;
    const std::unique_ptr<Matrix_PETSc> bc_mat;
    const std::unique_ptr<TimeMethod_GenAlpha> tmga;
    const std::unique_ptr<IFlowRate> flrate;
    const std::unique_ptr<PDNSolution> sol_base;

    void Print_convergence_info( const int &count, const double rel_err,
        const double abs_err ) const
    {
      SYS_T::commPrint("  === NR ite: %d, r_error: %e, a_error: %e \n",
          count, rel_err, abs_err);
    }

    void rescale_inflow_value( const double &stime,
        const ALocal_InflowBC * const &infbc,
        const PDNSolution * const &sol_base,
        PDNSolution * const &sol ) const;
};

#endif
