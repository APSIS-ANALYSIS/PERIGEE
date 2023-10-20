#ifndef PNONLINEAR_LINEARPDE_SOLVER_HPP
#define PNONLINEAR_LINEARPDE_SOLVER_HPP

#include "TimeMethod_GenAlpha.hpp"
#include "IPGAssem.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"
#include "PDNSolution_LinearPDE.hpp"

class PNonlinear_LinearPDE_Solver
{
  public:
    PNonlinear_LinearPDE_Solver( 
        const double &input_nrtol, const double &input_natol,
        const double &input_ndtol, const int &input_max_iteration,
        const int &input_renew_freq,
        const int &input_renew_threshold = 4 );

    ~PNonlinear_LinearPDE_Solver();

    int get_non_max_its() const {return nmaxits;}

    void print_info() const;

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
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const Matrix_PETSc * const &bc_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PDNSolution * const &dot_sol,
        PDNSolution * const &sol,
        bool &conv_flag, int &nl_counter ) const;
    
    void GenAlpha_Solve_Elastodynamics(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &pre_dot_disp,
        const PDNSolution * const &pre_dot_velo,
        const PDNSolution * const &pre_disp,
        const PDNSolution * const &pre_velo,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const Matrix_PETSc * const &bc_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PDNSolution * const &dot_disp,
        PDNSolution * const &dot_velo,
        PDNSolution * const &disp,
        PDNSolution * const &velo,
        bool &conv_flag, int &nl_counter ) const;

  private:
    const double nr_tol, na_tol, nd_tol;
    const int nmaxits, nrenew_freq, nrenew_threshold;

    void Print_convergence_info( const int &count, const double rel_err,
        const double abs_err ) const
    {PetscPrintf(PETSC_COMM_WORLD,
        "  === NR ite: %d, r_error: %e, a_error: %e \n",
        count, rel_err, abs_err);}


};

#endif
