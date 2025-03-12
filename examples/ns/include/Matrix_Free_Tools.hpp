#ifndef MATRIX_FREE_TOOLS_HPP
#define MATRIX_FREE_TOOLS_HPP

#include "PGAssem_Block_NS_FEM_HERK.hpp"

namespace MF_T
{
  PetscErrorCode MF_MatMult(Mat shell, Vec x, Vec y)
  {
    void *ptr;
    PGAssem_Block_NS_FEM_HERK *user;
    MatShellGetContext(shell, &ptr);
    user = (PGAssem_Block_NS_FEM_HERK*) ptr;

    Vec x1, x2, y1, y2;
    PetscInt m1_local, m2_local;
    const PetscScalar *x_array, *y1_array, *y2_array;
    PetscScalar *y_array;

    MatGetLocalSize(user->subK[3], &m1_local, NULL);
    MatGetLocalSize(user->subK[0], &m2_local, NULL);

    // Create sub Vec
    VecCreate(PETSC_COMM_WORLD, &x1);
    VecSetSizes(x1, m1_local, PETSC_DECIDE);
    VecCreate(PETSC_COMM_WORLD, &x2);
    VecSetSizes(x2, m2_local, PETSC_DECIDE);
    VecCreate(PETSC_COMM_WORLD, &y1);
    VecSetSizes(y1, m1_local, PETSC_DECIDE);
    VecCreate(PETSC_COMM_WORLD, &y2);
    VecSetSizes(y2, m2_local, PETSC_DECIDE);
    
    VecSetUp(x1);
    VecSetUp(x2);
    VecSetUp(y1);
    VecSetUp(y2);

    // Split x into x1, x2
    VecGetArrayRead(x, &x_array);
    VecPlaceArray(x1, x_array);
    VecPlaceArray(x2, x_array + m1_local);

    // y1 = (A3 + A4) * x1 + A2 * x2
    MatMult(user->subK[3], x1, y1);
    MatMultAdd(user->subK[4], x1, y1, y1);
    MatMultAdd(user->subK[2], x2, y1, y1);

    // y2 = A1 * x1 + A0 * x2
    MatMult(user->subK[1], x1, y2);
    MatMultAdd(user->subK[0], x2, y2, y2);

    // Combine y1, y2 to y
    VecGetArrayRead(y1, &y1_array);
    VecGetArrayRead(y2, &y2_array);
    VecGetArray(y, &y_array);

    PetscMemcpy(y_array, y1_array, m1_local * sizeof(PetscScalar));
    PetscMemcpy(y_array + m1_local, y2_array, m2_local * sizeof(PetscScalar));

    // Release
    VecRestoreArrayRead(y1, &y1_array);
    VecRestoreArrayRead(y2, &y2_array);
    VecRestoreArray(y, &y_array);

    VecResetArray(x1);
    VecResetArray(x2);
    VecRestoreArrayRead(x, &x_array);

    // Destruction of sub vectors
    VecDestroy(&x1);
    VecDestroy(&x2);
    VecDestroy(&y1);
    VecDestroy(&y2);

    return 0;
  }
}

#endif