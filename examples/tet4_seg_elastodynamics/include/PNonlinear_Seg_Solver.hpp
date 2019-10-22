#ifndef PNONLINEAR_SEG_SOLVER_HPP
#define PNONLINEAR_SEG_SOLVER_HPP
// ==================================================================
// PNonlinear_Seg_Solver.hpp
// 
// This class implements the segregated algorithm for elastodynamics.
// This file resides in the tet4_seg_elastodynamic project folder.
//
// Author: Ju Liu
// Date: May 23 2017
// ==================================================================
#include "TimeMethod_GenAlpha.hpp"
#include "PGAssem_Seg_FEM.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "PLinear_Solver_DiagScale.hpp"
#include "Matrix_PETSc.hpp"
#include "PDNSolution_P_V_Mixed_Hyperelastic_3D.hpp"
#include "Seg_Sol_Tools.hpp"

class PNonlinear_Seg_Solver
{
  public:
    PNonlinear_Seg_Solver( const APart_Node * const &anode_ptr,
        const IPLocAssem * const &lassem_ptr,
        const FEANode * const &feanode_ptr,
        const double &input_nrtol, 
        const double &input_natol, const double &input_ndtol,
        const int &input_max_iteration, const int &input_renew_freq );

    ~PNonlinear_Seg_Solver();

    int get_non_max_its() const {return nmaxits;}

    void print_info() const;

    // --------------------------------------------------------------
    // Generalized-alpha Segregated algorithm for elastodynamics
    //
    // This solver solves the elastodynamics in a segregated algorithm.
    // The pressure-velocity pair is solved first and the displacement
    // is updated based on the velocity.
    // --------------------------------------------------------------
    void GenAlpha_Seg_solve(
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
        const ALocal_NodalBC * const &nbc_part,
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


    // --------------------------------------------------------------
    // Generalized-alpha Segregated algorithm 2 for elastodynamics
    //
    // This solver solves the elastodynamics in a segregated algorithm.
    // The pressure-velocity pair is solved first and the displacement
    // is updated based on the velocity.
    //
    // Compared with the previous solver, this one does not solve
    // dp-dv for the disp update. We simply set dp-dv = 0 in the first
    // stage. This one is much cheaper.
    // --------------------------------------------------------------
    void GenAlpha_Seg_solve_2(
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
        const ALocal_NodalBC * const &nbc_part,
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


    // --------------------------------------------------------------
    // This is the same nonlinear solver algorithm as GenAlpha_Seg_solve_2
    // with explicit diagonal scaling of the matrix problem.
    // For a new tangent-residual, we call SymmJacobi_MatVec_Scale
    // to scale the matrix and vector;
    // For a new residual, we call SymmJacobi_Vec_Scale to scale the
    // vector.
    // --------------------------------------------------------------
    void GenAlpha_Seg_solve_DiagScale(
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
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const Matrix_PETSc * const &bc_mat,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IPLocAssem * const &lassem_ptr,
        IPGAssem * const &gassem_ptr,
        PLinear_Solver_DiagScale * const &lsolver_ptr,
        PDNSolution * const &dot_sol,
        PDNSolution * const &sol,
        bool &conv_flag, int &nl_counter ) const;


  private:
    const double nr_tol;
    const double na_tol;
    const double nd_tol;
    const int nmaxits;
    const int nrenew_freq;

    // vector container for the step update in the smaller matrix problem
    PDNSolution * dot_P_V_step;

    void Print_convergence_info( const int &count, const double rel_err,
        const double abs_err ) const
    {PetscPrintf(PETSC_COMM_WORLD,
        "  === NR ite: %d, r_error: %e, a_error: %e \n",
        count, rel_err, abs_err);}
};

#endif
