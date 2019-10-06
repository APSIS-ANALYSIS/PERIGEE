#ifndef PLINEAR_SOLVER_2X2BLOCK_BIPN_HPP
#define PLINEAR_SOLVER_2X2BLOCK_BIPN_HPP
// ==================================================================
// PLinear_Solver_2x2Block_BIPN.hpp
//
// This is the solver for 2x2 Block matrix using the Bi-partitioned
// iterative algorithm.
//
// Author: Ju Liu
// Date: Feb. 22 2018
// ==================================================================
#include "IPLinear_Solver_2x2Block.hpp"
#include "Math_Tools.hpp"

class PLinear_Solver_2x2Block_BIPN : public IPLinear_Solver_2x2Block
{
  public:
    PLinear_Solver_2x2Block_BIPN( const double &rtol0,
        const double &atol0, const double &dtol0, const int &maxit0,
        const double &rtol1, const double &atol1, const double &dtol1,
        const int &maxit1, const double &b_rtol, const int &b_maxits,
        const IPGAssem_2x2Block * const &gAssem_ptr );

    virtual ~PLinear_Solver_2x2Block_BIPN();

    virtual void print_info() const
    {
      PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
      PetscPrintf(PETSC_COMM_WORLD, "BIPN solver: \n");
      PetscPrintf(PETSC_COMM_WORLD, "BIPN rtol: %e \n", bipn_tol);
      PetscPrintf(PETSC_COMM_WORLD, "BIPN max its: %d \n", bipn_maxits);
      solver_0 -> Info();
      solver_1 -> Info();
      PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
    }

    virtual void SetLHS_RHS( IPGAssem_2x2Block const * const &gAssem_ptr );

    virtual void DestroyOperators();

    virtual void Solve_RHS( IPGAssem_2x2Block const * const &gAssem_ptr,
        PDNSolution * const &sol_0, PDNSolution * const &sol_1 );

  private:
    const double bipn_tol; // tolerance for BIPN conv. criterion
    const int bipn_maxits; // maximum BIPN iteration number

    // BIPN vectors forthe residual in each component
    Vec * Rmu, * Rmp, * Rcu, * Rcp;

    // Solutions at each BIPN iteration step
    Vec * y0, * y1;

    // Schur complement Mat
    Mat S;

    // Least square equation
    MATH_T::Matrix_SymPos_dense * A;
    MATH_T::Matrix_SymPos_dense * Asol;

    double * b; // RHS for the least square problem
    double * alpha; // solution for the least square problem
};

#endif
