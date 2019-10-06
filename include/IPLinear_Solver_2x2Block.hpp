#ifndef IPLINEAR_SOLVER_2X2BLOCK_HPP
#define IPLINEAR_SOLVER_2X2BLOCK_HPP
// ==================================================================
// IPLinear_Solver_2x2Block.hpp
//
// Interface for Parallel linear solver context for 2x2 Block system.
//
// This solver is used together with the matrices in IPGAssem_2x2Block.
// 
// Usage:
// case 1: The tangent matrix and residual vectors are generated.
//         step 1 (optional) SymmJacobi_MatVec_Scale
//         step 2 DestroyOperators
//         step 3 SetLHS_RHS
//         step 4 SolveRHS
//
// case 2: The residual vectors are generated
//         step 1 (optional) SymmJacobi_Vec_Scale
//         step 2 SetRHS
//         step 3 SolveRHS
//
// Author: Ju Liu
// Date: Feb. 21 2018
// ==================================================================
#include "IPGAssem_2x2Block.hpp"
#include "PLinear_Solver_2x2Block.hpp"

class IPLinear_Solver_2x2Block
{
  public:
    // --------------------------------------------------------------
    // Constructor:
    // Setup the linear solver solver_0 & solver_1.
    // Allocate the vectors G0, G1, D0, & D1.
    // --------------------------------------------------------------
    IPLinear_Solver_2x2Block(
        const double &rtol0, const double &atol0,
        const double &dtol0, const int &maxit0,
        const double &rtol1, const double &atol1,
        const double &dtol1, const int &maxit1,
        const IPGAssem_2x2Block * const &gAssem_ptr );


    // --------------------------------------------------------------
    // Destructor
    // --------------------------------------------------------------
    virtual ~IPLinear_Solver_2x2Block();


    // --------------------------------------------------------------
    // Print information of the solver context
    // --------------------------------------------------------------
    virtual void print_info() const = 0;


    // --------------------------------------------------------------
    // Assign the matrix to the solvers
    // Typically, the solver_0 is associated with an approximated
    // Schur matrix. We leave this function to the derived class for
    // specific choice of the construction of the Schur matrix.
    // --------------------------------------------------------------
    virtual void SetLHS_RHS( IPGAssem_2x2Block const * const &gAssem_ptr )
    {
      SYS_T::commPrint("Warning: SetOperators() is not implemented.\n");
    }


    // --------------------------------------------------------------
    // DestroyOperators.
    // In the set operator function, if one generate new Mat object
    // to approximating certain matrices, one may need to destroy
    // these Mat objects after usage.
    // --------------------------------------------------------------
    virtual void DestroyOperators()
    {
      SYS_T::commPrint("Warning: DestroyOperators() is not implemented.\n");
    } 


    // --------------------------------------------------------------
    // Assign the vectors to the solver context.
    // Copy the content of G_0 & G_1 from the global assembly class
    // to the G0 & G1 vectors, which belong to the solver class.
    // --------------------------------------------------------------
    void SetRHS( IPGAssem_2x2Block const * const &gAssem_ptr )
    { 
      VecCopy(gAssem_ptr->G_0, G0);
      VecCopy(gAssem_ptr->G_1, G1);
    }


    // --------------------------------------------------------------
    // Obtain the diagonal and do symmetric Jacobi scaling on the
    // four Mat objects.
    // D0 = diag(K00)^{-0.5}, D1 = diag(K11)^{-0.5}.
    // K00 <-- D0  K00  D0  ,  K01 <-- D0  K01  D1
    // K10 <-- D1  K10  D0  ,  K11 <-- D1  K11  D1
    // G0  <-- D0  G0       ,  G1  <-- D1  G1
    // Note: This will modify the gAssem_ptr object.
    // --------------------------------------------------------------
    void SymmJacobi_MatVec_Scale( IPGAssem_2x2Block * const &gAssem_ptr );


    // --------------------------------------------------------------
    // Scale the RHS vectors by
    // G0  <-- D0  G0       ,  G1  <-- D1  G1
    // One should be sure that the D0, D1 has been obtained before
    // calling this function.
    // --------------------------------------------------------------
    void SymmJacobi_Vec_Scale( IPGAssem_2x2Block * const &gAssem_ptr )
    {
      PETSc_T::DiagonalScale( D0, gAssem_ptr->G_0 );
      PETSc_T::DiagonalScale( D1, gAssem_ptr->G_1 );
    }


    // --------------------------------------------------------------
    // Restore the Matrices
    // D0 = 1.0 / D0, D1 = 1.0 / D1,
    // K00 <-- D0  K00  D0  ,  K01 <-- D0  K01  D1
    // K10 <-- D1  K10  D0  ,  K11 <-- D1  K11  D1
    // Note: 1. This will modify the gAssem_ptr object.
    //       2. In practical computation, this function is rarely needed, 
    //          since the original matrix, in most time, will not be
    //          needed after the solution is obtained. This may be used
    //          if one wants to calculate the actual residual of the
    //          original linear system.
    // --------------------------------------------------------------
    void SymmJacobi_MatVec_Restore( IPGAssem_2x2Block * const &gAssem_ptr );


    // --------------------------------------------------------------
    // Solve the block matrix problem with the RHS given in G0 and G1.
    // One has to call SetOperators before this function to setup
    // matrices; One has to call setRHS before this function to setup
    // the right-hand side vectors G0, G1. 
    // --------------------------------------------------------------
    virtual void Solve_RHS( IPGAssem_2x2Block const * const &gAssem_ptr,
        PDNSolution * const &sol_0, PDNSolution * const &sol_1 )
    {
      SYS_T::commPrint("Warning: Solve_RHS() is not implemented.\n");
    }


  protected:
    // Diagonal entries to the power -0.5 for K00 and K11
    Vec D0, D1;

    // Right-hand side vectors ( copied from the global assembly class )
    Vec G0, G1;

    // Linear solver for the K_00 matrix
    PLinear_Solver_PETSc * solver_0;

    // Linear solver for the K_11 matrix
    PLinear_Solver_PETSc * solver_1;
};

#endif
