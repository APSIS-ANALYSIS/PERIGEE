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
  //PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  
  //PetscAssertAbort( false, PETSC_COMM_WORLD, 60, "test\n");

  //PetscFinalize();
 
  std::vector<double> a;

  a.push_back(1.0);
  a.push_back(1.0);
  a.push_back(1.0);
  a.clear();

  //VEC_T::clean(a);

  std::cout<<a.size()<<'\t'<<a.capacity()<<'\n'; 
  
  return EXIT_SUCCESS;
}

// EOF
