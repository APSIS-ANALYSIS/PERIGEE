#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "SymmMatrix_3x3.hpp"
#include "HDF5_Reader.hpp"
#include "PostVectSolution.hpp"
#include "Tensor4_3D.hpp"

int main(int argc, char *argv[])
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  
  int a = 2;

  ASSERT( a == 1, "a is not 1" );


  std::cout<<a<<'\n';

  PetscFinalize();
 
  return EXIT_SUCCESS;
}

// EOF
