#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "SymmMatrix_3x3.hpp"
#include "HDF5_Reader.hpp"
#include "PostVectSolution.hpp"
#include "Tensor4_3D.hpp"
#include "Mesh_Tet.hpp"

int main(int argc, char *argv[])
{
  Mesh_Tet * mesh = new Mesh_Tet(100, 201, 3);

  mesh -> print_info();

  delete mesh;
  
  return EXIT_SUCCESS;
}

// EOF
