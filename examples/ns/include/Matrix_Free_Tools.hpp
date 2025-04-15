#ifndef MATRIX_FREE_TOOLS_HPP
#define MATRIX_FREE_TOOLS_HPP

#include "PGAssem_Block_NS_FEM_HERK.hpp"
#include "PLinear_Solver_PETSc.hpp"

namespace MF_T
{
  struct SolverContext
  {
    PGAssem_Block_NS_FEM_HERK *const gloAssem;
    const std::unique_ptr<PLinear_Solver_PETSc> lsolver_A;
    const std::unique_ptr<PLinear_Solver_PETSc> lsolver_S;

    // Constructor
    SolverContext(PGAssem_Block_NS_FEM_HERK * in_gloAssem,
        std::unique_ptr<PLinear_Solver_PETSc> in_lsolver_A, 
        std::unique_ptr<PLinear_Solver_PETSc> in_lsolver_S)
    : gloAssem(in_gloAssem), 
      lsolver_A(std::move(in_lsolver_A)), 
      lsolver_S(std::move(in_lsolver_S))
    {}
  };

  PetscErrorCode MF_MatMult(Mat shell, Vec x, Vec y)
  {
    void *ptr;
    PGAssem_Block_NS_FEM_HERK *user;
    MatShellGetContext(shell, &ptr);
    user = (PGAssem_Block_NS_FEM_HERK*) ptr;

    Vec x1, x2, y1, y2, tmp1, tmp2;

    const PetscScalar coef = user->Get_tangent_alpha_RK();

    // Split the VectNest into subVec
    VecNestGetSubVec(x, 0, &x1);
    VecNestGetSubVec(x, 1, &x2);
    VecNestGetSubVec(y, 0, &y1);
    VecNestGetSubVec(y, 1, &y2);

    VecDuplicate(x1, &tmp1);
    VecDuplicate(x2, &tmp2);

    MatMult(user->subK[3], x1, y1);      // y1 = A * x1
    MatMult(user->subK[4], x1, tmp1);    // tmp1 = A_tilde * x1
    VecAXPY(y1, coef, tmp1);             // y1 = A * x1 + coef * A_tilde * x1
    MatMult(user->subK[2], x2, tmp1);    // tmp1 = B * x2
    VecAXPY(y1, coef, tmp1);             // y1 = (A + coef * A_tilde) * x1 + coef * B * x2

    MatMult(user->subK[1], x1, y2);      // y2 = C * x1
    MatMult(user->subK[0], x2, tmp2);    // tmp2 = D * x2
    VecAXPY(y2, coef, tmp2);             // y2 = C * x1 + coef * D * x2
        
    // Destruction of sub vectors
    VecDestroy(&tmp1);
    VecDestroy(&tmp2);

    return 0;
  }

  PetscErrorCode SetupApproxSchur(PGAssem_Block_NS_FEM_HERK *const user, Mat &S_approx)
  {
    // Schur complement approximation: S = D - C inv(DIAGFORM(A)) B
    Vec diag;
    PetscInt mA_local;

    MatGetLocalSize(user->subK[3], &mA_local, NULL);

    // inverse of diagonal of A
    VecCreate(PETSC_COMM_WORLD, &diag);   
    VecSetSizes(diag, mA_local, PETSC_DETERMINE);
    VecSetType(diag, VECMPI);
    // MatGetRowSum(user->subK[3], diag); // Replace MatGetDiagonal() with MatGetRowSum()
    MatGetDiagonal(user->subK[3], diag);
    VecReciprocal(diag);

    MatDiagonalScale(user->subK[2], diag, NULL); // overwrites B = subK[2]) 
    MatMatMult(user->subK[1], user->subK[2], MAT_INITIAL_MATRIX, PETSC_DETERMINE, &S_approx);

    MatScale(S_approx, -1.0);
    MatAXPY(S_approx, 1.0, user->subK[0], DIFFERENT_NONZERO_PATTERN);  // S_approx = D - S_approx

    // restore B = subK[2]
    // MatGetRowSum(user->subK[3], diag); // Restore with row sums again
    MatGetDiagonal(user->subK[3], diag);
    MatDiagonalScale(user->subK[2], diag, NULL);

    VecDestroy(&diag);

    return 0;
  }
  
  PetscErrorCode MF_PCSchurApply(PC pc, Vec x, Vec y)
  {
  #ifdef PETSC_USE_LOG
    PetscLogEvent A_solve, S_solve;
    PetscClassId classid_solve;
    PetscClassIdRegister("PCsolve", &classid_solve);
    PetscLogEventRegister("A_solve", classid_solve, &A_solve);
    PetscLogEventRegister("S_solve", classid_solve, &S_solve);
  #endif

    void *ptr;
    SolverContext *ctx;
    PCShellGetContext(pc, &ptr);
    ctx = (SolverContext*) ptr;    

    Mat B = ctx->gloAssem->subK[2], C = ctx->gloAssem->subK[1];

    Vec x1, x2, y1, y2, z1, z2, tmp1, tmp2;

    // Split x into x1, x2
    VecNestGetSubVec(x, 0, &x1);
    VecNestGetSubVec(x, 1, &x2);
    
    // Split y into y1, y2
    VecNestGetSubVec(y, 0, &y1);
    VecNestGetSubVec(y, 1, &y2);

    VecDuplicate(x1, &z1);
    VecDuplicate(x2, &z2);
    VecDuplicate(x1, &tmp1);
    VecDuplicate(x2, &tmp2);

  #ifdef PETSC_USE_LOG
    PetscLogEventBegin(A_solve, 0,0,0,0);
  #endif 
    // Step 1: Compute z1 = A^{-1} x1
    ctx->lsolver_A->Solve(x1, z1, false);
  #ifdef PETSC_USE_LOG
    PetscLogEventEnd(A_solve,0,0,0,0);
  #endif

    // Step 2: Compute z2 = x2 - C * z1
    MatMult(C, z1, tmp2);
    VecWAXPY(z2, -1.0, tmp2, x2);

  #ifdef PETSC_USE_LOG
    PetscLogEventBegin(S_solve, 0,0,0,0);
  #endif 
    // Step 3: Compute y2 = S^{-1} z2
    ctx->lsolver_S->Solve(z2, y2, false);
  #ifdef PETSC_USE_LOG
    PetscLogEventEnd(S_solve,0,0,0,0);
  #endif

    // Step 4: Compute y1 = z1 - A^{-1} B y2 = A^{-1} x1 - A^{-1} B y2
    MatMult(B, y2, tmp1);
  #ifdef PETSC_USE_LOG
    PetscLogEventBegin(A_solve, 0,0,0,0);
  #endif 
    ctx->lsolver_A->Solve(tmp1, tmp1, false); 
  #ifdef PETSC_USE_LOG
    PetscLogEventEnd(A_solve,0,0,0,0);
  #endif
    VecWAXPY(y1, -1.0, tmp1, z1);

    VecDestroy(&z1);
    VecDestroy(&z2);
    VecDestroy(&tmp1);
    VecDestroy(&tmp2);

    return 0;
  }
}

#endif
