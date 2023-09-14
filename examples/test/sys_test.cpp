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
#include "DataVecStr.hpp"
#include "Tet_Tools.hpp"
#include "Hex_Tools.hpp"
#include "IIEN.hpp"
#include "Gmsh_FileIO.hpp"

int main(int argc, char *argv[])
{
  for(int ii=0; ii<100; ++ii){
  Vector_3 a, b, c, d;
  a.gen_rand(0.5, 1.5);
  b.gen_rand(1.8, 2.1);
  c.gen_rand(0.3, 0.9);
  d.gen_rand(1.1, 3.0);

  std::array<Vector_3, 4> arr {{a, b, c, d}};

  double val1 = TET_T::get_aspect_ratio(arr);

  TET_T::Tet4 cell(arr);

  double val2 = cell.get_aspect_ratio();

  std::cout<<val1<<'\t'<<val2<<'\t'<<val1 - val2<<'\n';
  }
  return EXIT_SUCCESS;
}

// EOF
