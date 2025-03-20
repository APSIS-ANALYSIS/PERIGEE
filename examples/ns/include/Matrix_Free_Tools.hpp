#ifndef MATRIX_FREE_TOOLS_HPP
#define MATRIX_FREE_TOOLS_HPP

#include "PGAssem_Block_NS_FEM_HERK.hpp"
#include "PLinear_Solver_PETSc.hpp"

namespace MF_T
{

  typedef struct 
  {
    PGAssem_Block_NS_FEM_HERK *gloAssem;
    std::unique_ptr<PLinear_Solver_PETSc> lsolver_A;
    std::unique_ptr<PLinear_Solver_PETSc> lsolver_S;
  } SolverContext;

  PetscErrorCode MF_MatMult(Mat shell, Vec x, Vec y)
  {
    void *ptr;
    PGAssem_Block_NS_FEM_HERK *user;
    MatShellGetContext(shell, &ptr);
    user = (PGAssem_Block_NS_FEM_HERK*) ptr;

    Vec x1, x2, y1, y2, tmp1, tmp2;
    // PetscInt m1_local, m2_local;

    const PetscScalar coef = user->Get_tangent_alpha_RK();

    // MatGetLocalSize(user->subK[3], &m1_local, NULL);
    // MatGetLocalSize(user->subK[0], &m2_local, NULL);

    // Create sub Vec
    // VecCreate(PETSC_COMM_WORLD, &x1);
    // VecSetSizes(x1, m1_local, PETSC_DECIDE);
    // VecCreate(PETSC_COMM_WORLD, &x2);
    // VecSetSizes(x2, m2_local, PETSC_DECIDE);
    // VecCreate(PETSC_COMM_WORLD, &y1);
    // VecSetSizes(y1, m1_local, PETSC_DECIDE);
    // VecCreate(PETSC_COMM_WORLD, &y2);
    // VecSetSizes(y2, m2_local, PETSC_DECIDE);
    
    // VecSetUp(x1);
    // VecSetUp(x2);
    // VecSetUp(y1);
    // VecSetUp(y2);

    // Split x into x1, x2  
    // PetscInt start1, start2;
    // IS is1, is2;

    // VecGetOwnershipRange(x, &start1, NULL);
    // start2 = start1 + m1_local; 

    // ISCreateStride(PETSC_COMM_WORLD, m1_local, start1, 1, &is1);
    // ISCreateStride(PETSC_COMM_WORLD, m2_local, start2, 1, &is2);

    // VecGetSubVector(x, is1, &x1);
    // VecGetSubVector(x, is2, &x2);

    // Split the VectNest into subVec
    VecNestGetSubVec(x, 0, &x1);
    VecNestGetSubVec(x, 1, &x2);
    VecNestGetSubVec(y, 0, &y1);
    VecNestGetSubVec(y, 1, &y2);

    VecDuplicate(x1, &tmp1);
    VecDuplicate(x2, &tmp2);

// //////////////
// PetscInt x1_size, x2_size, y1_size, y2_size;
// PetscInt x1_local, x2_local, y1_local, y2_local;

// VecGetSize(x1, &x1_size);
// VecGetSize(x2, &x2_size);
// VecGetSize(y1, &y1_size);
// VecGetSize(y2, &y2_size);

// VecGetLocalSize(x1, &x1_local);
// VecGetLocalSize(x2, &x2_local);
// VecGetLocalSize(y1, &y1_local);
// VecGetLocalSize(y2, &y2_local);

// PetscPrintf(PETSC_COMM_WORLD, 
//     "x1: global = %d, local = %d\n"
//     "x2: global = %d, local = %d\n"
//     "y1: global = %d, local = %d\n"
//     "y2: global = %d, local = %d\n",
//     x1_size, x1_local, 
//     x2_size, x2_local, 
//     y1_size, y1_local, 
//     y2_size, y2_local);

// for (int i = 0; i < 5; i++) {
//       PetscInt m, n, m_local, n_local;
//       MatGetSize(user->subK[i], &m, &n);
//       MatGetLocalSize(user->subK[i], &m_local, &n_local);
  
//       PetscPrintf(PETSC_COMM_WORLD, 
//           "subK[%d]: global size = %d x %d, local size = %d x %d\n",
//           i, m, n, m_local, n_local);
//   }

// //////////////

    MatMult(user->subK[3], x1, y1);      // y1 = A * x1
    MatMult(user->subK[4], x1, tmp1);    // tmp1 = A_tilde * x1
    VecAXPY(y1, coef, tmp1);             // y1 = A * x1 + coef * A_tilde * x1
    MatMult(user->subK[2], x2, tmp1);    // tmp1 = B * x2
    VecAXPY(y1, coef, tmp1);             // y1 = (A + coef * A_tilde) * x1 + coef * B * x2

    // Split y into y1, y2
    // VecGetSubVector(y, is1, &y1);
    // VecGetSubVector(y, is2, &y2);

    MatMult(user->subK[1], x1, y2);      // y2 = C * x1
    MatMult(user->subK[0], x2, tmp2);    // tmp2 = D * x2
    VecAXPY(y2, coef, tmp2);             // y2 = C * x1 + coef * D * x2
    
    // Restore
    // VecRestoreSubVector(y, is1, &y1);
    // VecRestoreSubVector(y, is2, &y2);

    // Destruction of sub vectors
    // VecDestroy(&x1);
    // VecDestroy(&x2);
    // VecDestroy(&y1);
    // VecDestroy(&y2);
    VecDestroy(&tmp1);
    VecDestroy(&tmp2);
    // ISDestroy(&is1);
    // ISDestroy(&is2);

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
    MatGetDiagonal(user->subK[3], diag);
    VecReciprocal(diag);

    MatDiagonalScale(user->subK[2], diag, NULL); // overwrites B = subK[2]) 
    MatMatMult(user->subK[1], user->subK[2], MAT_INITIAL_MATRIX, PETSC_DETERMINE, &S_approx);

    MatScale(S_approx, -1.0);
    MatAXPY(S_approx, 1.0, user->subK[0], DIFFERENT_NONZERO_PATTERN);  // S_approx = D - S_approx

    // restore B = subK[2]
    MatGetDiagonal(user->subK[3], diag);
    MatDiagonalScale(user->subK[2], diag, NULL);

    VecDestroy(&diag);

    return 0;
  }

  // PetscErrorCode MatMult_S(Mat S, Vec x, Vec y)
  // {
  //   void *ptr;
  //   PGAssem_Block_NS_FEM_HERK *user;   // 这里应该传一个包含PGA和ksp_A的结构体进去
  //   MatShellGetContext(S, &ptr);
  //   user = (PGAssem_Block_NS_FEM_HERK*) ptr;

  //   Mat B, C, D;
  //   Vec temp, Ax;

  //   B = user->subK[2];
  //   C = user->subK[1];
  //   D = user->subK[0];

  //   VecDuplicate(x, &temp);
  //   VecDuplicate(x, &Ax);

  //   // A^{-1} B x
  //   MatMult(B, x, temp);     // temp = B * x
  //   KSPSolve(user->ksp_A, temp, Ax);  // Ax = A^{-1} * (B x)

  //   // C * A^{-1} B x
  //   MatMult(C, Ax, temp);  // temp = C * Ax

  //   // y = D x - temp
  //   MatMult(D, x, y);
  //   VecAXPY(y, -1.0, temp);

  //   VecDestroy(&temp);
  //   VecDestroy(&Ax);

  //   return 0;
  // }

  // PetscErrorCode S_PCSetUp(PC pc, PGAssem_Block_NS_FEM_HERK *user)
  // {
  //   Mat S;

  //   PetscInt mD_local;

  //   MatGetLocalSize(user->subK[0], &mD_local, NULL);

  //   MatCreateShell(PETSC_COMM_WORLD, mD_local, mD_local, PETSC_DETERMINE, PETSC_DETERMINE, user, &S);
  //   MatShellSetOperation(S, MATOP_MULT, (void(*)(void))MatMult_S);

  //   lsolver_S->Setoperator(S); 

  //   return 0;
  // }

  // PetscErrorCode MyPCSetUp(PC pc)
  // {
  //   PGAssem_Block_NS_FEM_HERK *user;
  //   PetscCall(PCShellGetContext(pc, (void**)&user));

  //   Mat A = user->subK[3];
  //   lsolver_A->Setoperator(A);  

  //   Mat S;

  //   PetscInt mD_local;

  //   MatGetLocalSize(user->subK[0], &mD_local, NULL);

  //   MatCreateShell(PETSC_COMM_WORLD, mD_local, mD_local, PETSC_DETERMINE, PETSC_DETERMINE, user, &S);
  //   MatShellSetOperation(S, MATOP_MULT, (void(*)(void))MatMult_S);

  //   lsolver_S->Setoperator(S);

  //   return 0;
  // }

  
  PetscErrorCode MF_PCSchurApply(PC pc, Vec x, Vec y)
  {
    void *ptr;
    SolverContext *ctx;
    PCShellGetContext(pc, &ptr);
    ctx = (SolverContext*) ptr;    

    Mat B = ctx->gloAssem->subK[2], C = ctx->gloAssem->subK[1];
    // Mat A = ctx->gloAssem->subK[3], B = ctx->gloAssem->subK[2], C = ctx->gloAssem->subK[1], D = ctx->gloAssem->subK[0];

    Vec x1, x2, y1, y2, z1, z2, w1, w2, tmp1, tmp2;
    // IS is1, is2;
    // PetscInt m1_local, m2_local;

    // MatGetLocalSize(A, &m1_local, NULL);
    // MatGetLocalSize(D, &m2_local, NULL);

    // PetscInt start1, start2;
    // VecGetOwnershipRange(x, &start1, NULL);
    // start2 = start1 + m1_local;

    // ISCreateStride(PETSC_COMM_WORLD, m1_local, start1, 1, &is1);
    // ISCreateStride(PETSC_COMM_WORLD, m2_local, start2, 1, &is2);

    // VecGetSubVector(x, is1, &x1);
    // VecGetSubVector(x, is2, &x2);
    VecNestGetSubVec(x, 0, &x1);
    VecNestGetSubVec(x, 1, &x2);

    VecDuplicate(x1, &z1);
    VecDuplicate(x2, &z2);
    VecDuplicate(x1, &w1);
    VecDuplicate(x2, &w2);
    VecDuplicate(x1, &tmp1);
    VecDuplicate(x2, &tmp2);

    // Step 1: Compute z1 = A^{-1} x1
    ctx->lsolver_A->Solve(x1, z1, false);

    // Step 2: Compute z2 = x2 - C * z1
    MatMult(C, z1, tmp2);
    VecWAXPY(z2, -1.0, tmp2, x2);

    // Step 3: Compute w2 = S^{-1} z2
    ctx->lsolver_S->Solve(z2, w2, false);

    // Step 4: Compute w1 = z1 - A^{-1} B w2 = A^{-1} x1 - A^{-1} B w2
    MatMult(B, w2, tmp1);
    ctx->lsolver_A->Solve(tmp1, tmp1, false); 
    VecWAXPY(w1, -1.0, tmp1, z1);

    // Combine w1 and w2 into y
    // VecGetSubVector(y, is1, &y1);
    // VecGetSubVector(y, is2, &y2);
    // Split y into y1, y2
    VecNestGetSubVec(y, 0, &y1);
    VecNestGetSubVec(y, 1, &y2);

    VecCopy(w1, y1);
    VecCopy(w2, y2);

    // Restore and clean up
    // VecRestoreSubVector(x, is1, &x1);
    // VecRestoreSubVector(x, is2, &x2);
    // VecRestoreSubVector(y, is1, &y1);
    // VecRestoreSubVector(y, is2, &y2);

    // ISDestroy(&is1);
    // ISDestroy(&is2);
    VecDestroy(&z1);
    VecDestroy(&z2);
    VecDestroy(&w1);
    VecDestroy(&w2);
    VecDestroy(&tmp1);
    VecDestroy(&tmp2);

    return 0;
  }
}

#endif
