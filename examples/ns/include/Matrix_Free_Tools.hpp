#ifndef MATRIX_FREE_TOOLS_HPP
#define MATRIX_FREE_TOOLS_HPP

#include "PGAssem_Block_NS_FEM_HERK.hpp"
#include "PLinear_Solver_PETSc.hpp"

namespace MF_T
{
  struct SolverContext
  {
#ifdef PETSC_USE_LOG
    PetscLogEvent A_solve, S_solve;
#endif

    PGAssem_Block_NS_FEM_HERK *const gloAssem;
    const std::unique_ptr<PLinear_Solver_PETSc> lsolver_A;
    const std::unique_ptr<PLinear_Solver_PETSc> lsolver_S;

    Vec tmp1, tmp2;    

    // Constructor
    SolverContext(PGAssem_Block_NS_FEM_HERK * in_gloAssem,
        std::unique_ptr<PLinear_Solver_PETSc> in_lsolver_A, 
        std::unique_ptr<PLinear_Solver_PETSc> in_lsolver_S)
    : gloAssem(in_gloAssem), 
      lsolver_A(std::move(in_lsolver_A)), 
      lsolver_S(std::move(in_lsolver_S))
    {
      VecDuplicate(gloAssem->subG[0], &tmp1); 
      VecDuplicate(gloAssem->subG[1], &tmp2);      
    }

    // Destructor
    ~SolverContext()
    {
      VecDestroy(&tmp1);
      VecDestroy(&tmp2);
    }
  };

  PetscErrorCode MF_MatMult(Mat shell, Vec x, Vec y)
  {
    void *ptr;
    SolverContext *ctx;
    MatShellGetContext(shell, &ptr);
    ctx = (SolverContext*) ptr;

    Vec x1, x2, y1, y2;

    const PetscScalar coef = ctx->gloAssem->Get_tangent_alpha_RK();

    // Split the VectNest into subVec: x1 is velo; x2 is pres;
    VecNestGetSubVec(x, 0, &x1);
    VecNestGetSubVec(x, 1, &x2);
    VecNestGetSubVec(y, 0, &y1);
    VecNestGetSubVec(y, 1, &y2);

    // y1 = A * x1
    MatMult(ctx->gloAssem->subK[3], x1, y1);
    // tmp1 = B * x2
    MatMult(ctx->gloAssem->subK[2], x2, ctx->tmp1);
    // y1 = A * x1 + coef * B * x2
    VecAXPY(y1, coef, ctx->tmp1);

    // y2 = B^T * x1;
    MatMultTranspose(ctx->gloAssem->subK[2], x1, y2);
    // tmp2 = D * x2
    MatMult(ctx->gloAssem->subK[0], x2, ctx->tmp2);
    // y2 = B^T * x1 - coef * D * x2
    VecAXPY(y2, -coef, ctx->tmp2);
    // y2 = coef * (B^T * x1 - coef * D * x2)
    VecScale(y2, coef);

    return 0;
  }

  PetscErrorCode SetupApproxSchur(PGAssem_Block_NS_FEM_HERK *const user, Mat &S_approx)
  {
    // Schur complement approximation: S_approx = - D - B^T inv(DIAGFORM(A)) B
    Vec diag;
    PetscInt mA_local;

    MatGetLocalSize(user->subK[3], &mA_local, NULL);

    // inverse of diagonal of A
    VecCreate(PETSC_COMM_WORLD, &diag);   
    VecSetSizes(diag, mA_local, PETSC_DETERMINE);
    VecSetType(diag, VECMPI);
    MatGetRowSum(user->subK[3], diag); // Replace MatGetDiagonal() with MatGetRowSum()
    //MatGetDiagonal(user->subK[3], diag);
    VecReciprocal(diag);

    Mat invAB;
    MatDuplicate(user->subK[2], MAT_COPY_VALUES, &invAB);
    MatDiagonalScale(invAB, diag, NULL);  // Compute diag(A)^-1 * B  
    // Compute B^T inv(DIAGFORM(A)) B
    MatTransposeMatMult(user->subK[2], invAB, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &S_approx);

    MatScale(S_approx, -1.0); 
    // S_approx = - D - S_approx
    MatAXPY(S_approx, -1.0, user->subK[0], DIFFERENT_NONZERO_PATTERN);

    VecDestroy(&diag);
    MatDestroy(&invAB);

    return 0;
  }
  
  PetscErrorCode MF_PCSchurSetup(PC pc)
  {
    void *ptr;
    SolverContext *ctx;
    PCShellGetContext(pc, &ptr);
    ctx = (SolverContext*) ptr; 

#ifdef PETSC_USE_LOG
    PetscLogEventRegister("A_solve", KSP_CLASSID, &ctx->A_solve);
    PetscLogEventRegister("S_solve", KSP_CLASSID, &ctx->S_solve);
#endif

  return 0;
  }

  PetscErrorCode MF_PCSchurApply(PC pc, Vec x, Vec y)
  {
    void *ptr;
    SolverContext *ctx;
    PCShellGetContext(pc, &ptr);
    ctx = (SolverContext*) ptr;    

    Mat B = ctx->gloAssem->subK[2];
    const PetscScalar coef = ctx->gloAssem->Get_tangent_alpha_RK();
    
    Vec x1, x2, y1, y2;

    // Split x into x1, x2
    VecNestGetSubVec(x, 0, &x1);
    VecNestGetSubVec(x, 1, &x2);
    
    // Split y into y1, y2
    VecNestGetSubVec(y, 0, &y1);
    VecNestGetSubVec(y, 1, &y2);

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(ctx->A_solve, 0,0,0,0);
#endif 
    // Step 1: Compute y1 = A^{-1} x1
    ctx->lsolver_A->Solve(x1, y1, false);
#ifdef PETSC_USE_LOG
    PetscLogEventEnd(ctx->A_solve,0,0,0,0);
#endif

    // Step 2: Compute y2 = x2 - coef * B^T * y1
    MatMultTranspose(B, y1, ctx->tmp2);
    VecWAXPY(y2, -coef, ctx->tmp2, x2);

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(ctx->S_solve, 0,0,0,0);
#endif 

    // Step 3: Compute y2 = (1/coef^2) * S_approx^{-1} y2 
    ctx->lsolver_S->Solve(y2, y2, false);
#ifdef PETSC_USE_LOG
    PetscLogEventEnd(ctx->S_solve,0,0,0,0);
#endif     

    VecScale(y2, 1.0/(coef*coef));

    // Step 4: Compute y1 = y1 - coef * A^{-1} B y2 = A^{-1} x1 - coef * A^{-1} B y2
    MatMult(B, y2, ctx->tmp1);
#ifdef PETSC_USE_LOG
    PetscLogEventBegin(ctx->A_solve, 0,0,0,0);
#endif 
    ctx->lsolver_A->Solve(ctx->tmp1, ctx->tmp1, false);
#ifdef PETSC_USE_LOG
    PetscLogEventEnd(ctx->A_solve,0,0,0,0);
#endif

    VecAXPY(y1, -coef, ctx->tmp1);
    return 0;
  }
}

#endif
