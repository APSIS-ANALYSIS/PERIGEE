#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "SymmMatrix_3x3.hpp"
#include "HDF5_Reader.hpp"
#include "PostVectSolution.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"
#include "Mesh_Tet.hpp"
#include "IEN_FEM.hpp"
#include "Matrix_double_3by3_Array.hpp"
#include "Matrix_double_6by6_Array.hpp"
#include "VTK_Tools.hpp"
#include "NodalBC.hpp"
#include "Gmsh_FileIO.hpp"

int main(int argc, char *argv[])
{
  // test 1
  Gmsh_FileIO * GIO_1 = new Gmsh_FileIO( "fsi_cylinder.msh" );
  GIO_1 -> print_info();
  GIO_1 -> write_tet_h5(0, {0, 1, 2});
  GIO_1 -> write_tet_h5(1, {3, 4, 5});

  delete GIO_1;

  // test 2
  Gmsh_FileIO * GIO_2 = new Gmsh_FileIO( "fsi_beam.msh" );
  GIO_2 -> print_info();
  GIO_2 -> write_tet_h5(0, {0, 1, 2, 3, 4, 5});
  GIO_2 -> write_tet_h5(1, {7, 8, 9, 10, 11, 12});

  delete GIO_2;

  // test 3
  Gmsh_FileIO * GIO_3 = new Gmsh_FileIO( "cook_membrane.msh" );

  GIO_3 -> print_info();
  GIO_3 -> write_tri_h5(0, {0, 1, 2, 3}); // 2d problem

  delete GIO_3;


  return EXIT_SUCCESS;
}

// EOF
