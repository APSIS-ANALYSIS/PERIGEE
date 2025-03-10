#ifndef PNONLINEAR_LINEARPDE_SOLVER_HPP
#define PNONLINEAR_LINEARPDE_SOLVER_HPP
// ==================================================================
// PNonlinear_LinearPDE_Solver.hpp
// 
// Parallel nonlinear solver for linear PDE. 
//
// Date: Oct. 23 2023
// ==================================================================
#include "TimeMethod_GenAlpha.hpp"
#include "IPGAssem.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"
#include "PDNSolution_Transport.hpp"
#include "PDNSolution_Elastodynamics.hpp"

class PNonlinear_LinearPDE_Solver
{
  public:
    PNonlinear_LinearPDE_Solver( 
        std::unique_ptr<IPGAssem> in_gassem,
        std::unique_ptr<PLinear_Solver_PETSc> in_lsolver,
        std::unique_ptr<Matrix_PETSc> in_bc_mat,
        std::unique_ptr<TimeMethod_GenAlpha> in_tmga,
        const double &input_nrtol, const double &input_natol,
        const double &input_ndtol, const int &input_max_iteration,
        const int &input_renew_freq,
        const int &input_renew_threshold = 4 );

    ~PNonlinear_LinearPDE_Solver() = default;

    int get_non_max_its() const {return nmaxits;}

    void print_info() const;

    void print_lsolver_info() const {lsolver->print_info();}

    // --------------------------------------------------------------
    // GenAlpha_Solve_Transport:
    // This is a solver for transport equation.
    // --------------------------------------------------------------
    void GenAlpha_Solve_Transport(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &pre_dot_sol,
        const PDNSolution * const &pre_sol,
        PDNSolution * const &dot_sol,
        PDNSolution * const &sol,
        bool &conv_flag, int &nl_counter ) const;

    // --------------------------------------------------------------
    // GenAlpha_Solve_Elastodynamics:
    // This is a solver for elastodynamics equation.
    // --------------------------------------------------------------
    void GenAlpha_Solve_Elastodynamics(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &pre_dot_disp,
        const PDNSolution * const &pre_dot_velo,
        const PDNSolution * const &pre_disp,
        const PDNSolution * const &pre_velo,
        PDNSolution * const &dot_disp,
        PDNSolution * const &dot_velo,
        PDNSolution * const &disp,
        PDNSolution * const &velo,
        bool &conv_flag, int &nl_counter ) const;

  private:
    const double nr_tol, na_tol, nd_tol;
    const int nmaxits, nrenew_freq, nrenew_threshold;

    const std::unique_ptr<IPGAssem> gassem;
    const std::unique_ptr<PLinear_Solver_PETSc> lsolver;
    const std::unique_ptr<Matrix_PETSc> bc_mat;
    const std::unique_ptr<TimeMethod_GenAlpha> tmga;

    void Print_convergence_info( const int &count, const double rel_err,
        const double abs_err ) const
    {
      SYS_T::commPrint("  === NR ite: %d, r_error: %e, a_error: %e \n",
          count, rel_err, abs_err);
    }

};

#endif
