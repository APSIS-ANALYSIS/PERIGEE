#ifndef PNONLINEAR_FSI_SOLVER_HPP
#define PNONLINEAR_FSI_SOLVER_HPP
// ============================================================================
// PNonlinear_FSI_Solver.hpp
//
// This is the nonlienar solver for FSI problems.
//
// Author: Ju Liu
// Date: Jan 4 2022
// ============================================================================
#include "ICVFlowRate.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"
#include "PDNSolution_V.hpp"
#include "PDNSolution_P.hpp"

class PNonlinear_FSI_Solver
{
  public:
    PNonlinear_FSI_Solver( const double &input_nrtol, const double &input_natol,
        const double &input_ndtol, const int &input_max_iteration,
        const int &input_renew_freq );

    ~PNonlinear_FSI_Solver();

    int get_non_max_its() const {return nmaxits;}

    void print_info() const;

  private:
    const double nr_tol, na_tol, nd_tol;
    const int nmaxits, nrenew_freq;

    void Print_convergence_info( const int &count, const double &rel_err,
        const double &abs_err ) const
    {
      SYS_T::commPrint( "  === NR ite: %d, r_error: %e, a_error: %e \n", count, rel_err, abs_err);
    }

    void rescale_inflow_value( const double &stime,
        const ALocal_Inflow_NodalBC * const &infbc,
        const ICVFlowRate * const &flrate,
        const PDNSolution * const &sol_base,
        PDNSolution * const &sol ) const;
};

#endif
