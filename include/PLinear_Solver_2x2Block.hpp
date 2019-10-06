#ifndef PLINEAR_SOLVER_2X2BLOCK_HPP
#define PLINEAR_SOLVER_2X2BLOCK_HPP
// ==================================================================
// PLinear_Solver_2X2Block.hpp
//
// This is the solver for 2x2 Block matrix based on exact LU block
// factorization.
// [ K_11, K_10   [ yu    = [ G_1
//   K_01, K_00 ]   yp ]      G_0 ]
//
// Author: Ju Liu
// Date: Jan 24 2018
// ==================================================================
#include "PLinear_Solver_PETSc.hpp"
#include "IPGAssem_2x2Block.hpp"

class PLinear_Solver_2x2Block
{
  public:
    // --------------------------------------------------------------
    // Input the relative tolerance for the _0 solver and the _1 solver
    // _0 is associated with the Schur complement.
    // _1 is associated with the K_11 matrix
    // PGAssem_2x2Block_FEM provides the Vec formats for G0 and G1
    // --------------------------------------------------------------
    PLinear_Solver_2x2Block( const double &rtol0, 
        const double &atol0, const double &dtol0, const int &maxit0,
        const double &rtol1, const double &atol1, const double &dtol1, 
        const int &maxit1,
        const IPGAssem_2x2Block * const &gAssem_ptr );


    ~PLinear_Solver_2x2Block();


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
    // Solve for two segregated problems:
    //   S yp    = G_0 - K_01 inv(K_11) G_1
    //   K_11 yu = G1 - K_10 yp
    // wherein S := K_00 - K_01 inv(K_11) K_10 is the pressure Schur
    // complement.
    // Note:
    // By default, we solve the S matrix problem with cg and no 
    // preconditioner. One may specify the solver options with prefix
    // ls0_. We do not change the pc here since the S is not explicitly
    // assembled.
    // The inv(K_11) inside S can be specified by options with prefix
    // ls0_inner.
    // One may control the solver and preconditioner for inv K_11
    // by giving options from the command line with a prefix ls1_
    // --------------------------------------------------------------
    void Solve_ExaFac( const IPGAssem_2x2Block * const &gAssem_ptr,
        PDNSolution * const &sol_0, PDNSolution * const &sol_1 );

  private:
    Vec D0, D1; // Diagonal entries to the power -0.5 for K00 and K11

    Vec G0, G1;

    Mat S; // Schur Complement Mat context

    PLinear_Solver_PETSc * solver_0;

    PLinear_Solver_PETSc * solver_1;
};

#endif
