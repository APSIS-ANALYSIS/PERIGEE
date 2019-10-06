#ifndef PLINEAR_SOLVER_BIPN_2X2BLOCK_HPP
#define PLINEAR_SOLVER_BIPN_2X2BLOCK_HPP
// ==================================================================
// PLinear_Solver_BIPN_2x2Block.hpp
//
// This is the solver for 2x2 Block matrix based on Bi-partitioned
// iterative algorithm. The matrix problem is 
// [ K_11, K_10   [ yu    = [ G_1
//   K_01, K_00 ]   yp ]      G_0 ]
//
// Author: Ju Liu
// Date: Jan 29 2018
// ==================================================================
#include "PLinear_Solver_PETSc.hpp"
#include "IPGAssem_2x2Block.hpp"
#include "Math_Tools.hpp"

class PLinear_Solver_BIPN_2x2Block
{
  public:
    PLinear_Solver_BIPN_2x2Block( const double &rtol0,
        const double &atol0, const double &dtol0, const int &maxit0,
        const double &rtol1, const double &atol1, const double &dtol1,
        const int &maxit1, const double &b_rtol, const int &b_maxits,
        const IPGAssem_2x2Block * const &gAssem_ptr );


    ~PLinear_Solver_BIPN_2x2Block();

    
    // Print the info of the solvers
    void print_info() const
    {
      solver_0 -> Info();
      solver_1 -> Info();
    }


    // Symmetric Jacobi scaling on the 4 mat objects and 2 vec objects
    // D0 = diag(K00)^{-0.5}, D1 = diag(K11)^{-0.5}.
    // K00 <-- D0  K00  D0  ,  K01 <-- D0  K01  D1
    // K10 <-- D1  K10  D0  ,  K11 <-- D1  K11  D1
    // G0  <-- D0  G0       ,  G1  <-- D1  G1
    // This scaling will change the mat objects in IPGAssem_2x2Block.
    void SymmJacobi_Scale( IPGAssem_2x2Block * const &gAssem_ptr );

    
    // Restore the Matrices and Vectors in gAssem_ptr;
    // D0 = 1.0 / D0, D1 = 1.0 / D1,
    // K00 <-- D0  K00  D0  ,  K01 <-- D0  K01  D1
    // K10 <-- D1  K10  D0  ,  K11 <-- D1  K11  D1
    // G0  <-- D0  G0       ,  G1  <-- D1  G1
    void SymmJacobi_Restore( IPGAssem_2x2Block * const &gAssem_ptr );


    // --------------------------------------------------------------
    // Solve using Bi-partitioned two-stage algorithm.
    // Within this algorithm, we use Symmetric Jacobi preconditioner
    // for the K_00, K_01, K_10, K_11 matrices and G0, G1 vectors.
    // So the globassem content will be modified within this function.
    // Ref. M. Esmaily-Moghadam, et al. CMAME 2015
    // Note: This is one-pass solve, meaning we assign matrices inside
    // this function.
    // --------------------------------------------------------------
    void Solve_BIPN( IPGAssem_2x2Block * const &gAssem_ptr,
        PDNSolution * const &sol_0, PDNSolution * const &sol_1 );


  private:
    const double bipn_tol; // tolerance for BIPN conv. criterion
    const int bipn_maxits; // maximum BIPN iteration number

    Vec D0, D1; // Diagonal entries to the power -0.5 for K00 and K11

    Vec G0, G1;

    // BIPN vectors for the residual in each component.
    Vec * Rmu, * Rmp, * Rcu, * Rcp;

    // Solution at each BIPN iteration steps
    Vec * y0, * y1;

    Mat S; // Schur Complement Mat context

    //Mat Id; // Identity matrix with same layout with K_11

    PLinear_Solver_PETSc * solver_0;

    PLinear_Solver_PETSc * solver_1;

    // Least square equation and RHS
    MATH_T::Matrix_SymPos_dense * A;
    MATH_T::Matrix_SymPos_dense * Asol;

    double * b; // RHS for the least square problem
    double * alpha; // solution for the least square problem
};

#endif
