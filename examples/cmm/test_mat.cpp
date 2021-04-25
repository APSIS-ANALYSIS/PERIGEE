#include "Matrix_PETSc.hpp"
#include "PETSc_Tools.hpp"

int main( int argc , char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  Mat K;

  MatCreateAIJ(PETSC_COMM_WORLD, 5, 5, PETSC_DECIDE, PETSC_DECIDE,
      2, PETSC_NULL, 0, PETSC_NULL, &K);

  PETSc_T::Release_nonzero_err_str(K);

  for(int ii=0; ii<4; ++ii)
  {
    MatSetValue(K, ii, ii, 1.0, INSERT_VALUES);
    MatSetValue(K, ii, ii+1, 1.0, INSERT_VALUES);
  }

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

  PETSc_T::MatInfo_Display_global(K);

  MatDestroy(&K);

  PetscFinalize();
  return 0;
}
