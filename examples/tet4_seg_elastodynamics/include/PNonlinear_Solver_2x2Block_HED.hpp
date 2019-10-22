#ifndef PNONLINEAR_SOLVER_2X2BLOCK_HED_HPP
#define PNONLINEAR_SOLVER_2X2BLOCK_HED_HPP
// ==================================================================
// PNonlinear_Solver_2x2Block_HED.hpp
//
// This class implements the segregated algorithm for elastodynamics
// with the linear solver as a block-decomposed solver.
//
// Author: Ju Liu
// Date: Feb. 21 2018
// ==================================================================
#include "TimeMethod_GenAlpha.hpp"
#include "IPGAssem_2x2Block.hpp"
#include "IPLinear_Solver_2x2Block.hpp"
#include "PDNSolution_Disp_3D.hpp"
#include "PDNSolution_Pres_3D.hpp"

class PNonlinear_Solver_2x2Block_HED
{
  public:
    PNonlinear_Solver_2x2Block_HED( const APart_Node * const &pNode_ptr,
        const FEANode * const &feanode_ptr, const double &input_nrtol,
        const double &input_natol, const double &input_ndtol,
        const int &input_max_iteration, const int &input_renew_freq );

    ~PNonlinear_Solver_2x2Block_HED();
  
    int get_non_max_its() const {return nmaxits;}

    void print_info() const;

    void GenAlpha_solve(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &pre_dot_d,
        const PDNSolution * const &pre_dot_p,
        const PDNSolution * const &pre_dot_v,
        const PDNSolution * const &pre_sol_d,
        const PDNSolution * const &pre_sol_p,
        const PDNSolution * const &pre_sol_v,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IPLocAssem_2x2Block * const &lassem_ptr,
        IPGAssem_2x2Block * const &gassem_ptr,
        IPLinear_Solver_2x2Block * const &lsolver_ptr,
        PDNSolution * const &dot_d,
        PDNSolution * const &dot_p,
        PDNSolution * const &dot_v,
        PDNSolution * const &sol_d,
        PDNSolution * const &sol_p,
        PDNSolution * const &sol_v,
        bool &conv_flag, int &nl_counter ) const;


  private:
    const double nr_tol, na_tol, nd_tol;
    const int nmaxits, nrenew_freq;

    PDNSolution * step_velo, * sol_v_alpha, * dot_v_alpha;
    PDNSolution * step_pres, * sol_p_alpha, * dot_p_alpha;

    PDNSolution * sol_d_alpha, * dot_d_alpha;

    void Print_convergence_info( const int &count, const double rel_err,
        const double abs_err ) const
    { PetscPrintf(PETSC_COMM_WORLD, 
        "  === NR ite: %d, r_error: %e, a_error: %e \n",
        count, rel_err, abs_err); }
};

#endif
