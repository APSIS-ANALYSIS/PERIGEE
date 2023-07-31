#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "SymmMatrix_3x3.hpp"
#include "HDF5_Reader.hpp"
#include "PostVectSolution.hpp"
#include "Tensor4_3D.hpp"
#include "Mesh_Tet.hpp"
#include "IEN_FEM.hpp"

int main(int argc, char *argv[])
{
  Mesh_Tet * mesh = new Mesh_Tet(100, 201, 1);

  mesh -> print_info();

  std::vector<int> inien {1,3,4,5,6};

  IIEN * ien = new IEN_FEM(2, inien);

  ien -> print_info();


  delete ien;
  delete mesh;
  
  return EXIT_SUCCESS;
}

// EOF
