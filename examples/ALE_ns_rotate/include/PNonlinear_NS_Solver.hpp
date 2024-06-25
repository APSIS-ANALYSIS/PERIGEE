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
#include "ICVFlowRate.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "IPGAssem.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"
#include "PDNSolution_NS.hpp"
#include "PDNSolution_V.hpp"

class PNonlinear_NS_Solver
{
  public:
    PNonlinear_NS_Solver( const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const double &input_nrtol, const double &input_natol, 
        const double &input_ndtol, const int &input_max_iteration, 
        const int &input_renew_freq, 
        const int &input_renew_threshold = 4 );

    ~PNonlinear_NS_Solver();

    int get_non_max_its() const {return nmaxits;}

    void print_info() const;

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
        const PDNSolution * const &sol_base,
        const PDNSolution * const &pre_dot_sol,
        const PDNSolution * const &pre_sol,
        const PDNSolution * const &pre_disp_mesh,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ICVFlowRate * const flr_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_InflowBC * const &infnbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc,
        const ALocal_WeakBC * const &wbc_part,
        const Matrix_PETSc * const &bc_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        FEAElement * const &elementvs,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PDNSolution * const &dot_sol,
        PDNSolution * const &sol,
        PDNSolution * const &disp_mesh,
        bool &conv_flag, int &nl_counter ) const;

  private:
    const double nr_tol, na_tol, nd_tol;
    const int nmaxits, nrenew_freq, nrenew_threshold;

    // vector container for the step update in the smaller matrix problem
    PDNSolution * dot_step;

    void Print_convergence_info( const int &count, const double rel_err,
        const double abs_err ) const
    {PetscPrintf(PETSC_COMM_WORLD,
        "  === NR ite: %d, r_error: %e, a_error: %e \n",
        count, rel_err, abs_err);}

    void rescale_inflow_value( const double &stime,
        const ALocal_InflowBC * const &infbc,
        const ICVFlowRate * const &flrate,
        const PDNSolution * const &sol_base,
        PDNSolution * const &sol ) const;
};

#endif
