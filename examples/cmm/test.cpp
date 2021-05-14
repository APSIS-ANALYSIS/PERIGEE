// Test the new TET_T functions

#include "ElemBC_3D_tet_wall.hpp"

using namespace std;

int main( int argc, char *argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  std::vector<double> a {};
  
  std::cout<<a.size()<<std::endl;

  VEC_T::print(a);

  PetscFinalize();
  return 0;
}

//EOF
