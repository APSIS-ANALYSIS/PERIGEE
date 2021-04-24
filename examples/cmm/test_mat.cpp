#include "Matrix_PETSc.hpp"

int main( int argc , char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  Mat K;

  MatCreateAIJ(PETSC_COMM_WORLD, 5, 5, PETSC_DECIDE, PETSC_DECIDE,
      5, PETSC_NULL, 3, PETSC_NULL, &K);


  for(int ii=0; ii<5; ++ii)
  {
    MatSetValue(K, ii, ii, 1.0, INSERT_VALUES);
    MatSetValue(K, ii, ii+1, 1.0, INSERT_VALUES);
  }

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

  MatView(K, PETSC_VIEWER_STDOUT_WORLD);

  MatDestroy(&K);

  PetscFinalize();
  return 0;
}
