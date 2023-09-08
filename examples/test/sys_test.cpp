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

int main(int argc, char *argv[])
{
  Vector_3 pt0, pt1, pt2, pt3;
  pt0.gen_rand();
  pt1.gen_rand();
  pt2.gen_rand();
  pt3.gen_rand();

  double rad;
  auto cen = MATH_T::get_tet_sphere_info(pt0, pt1, pt2, pt3, rad);

  double x, y, z, r;
  MATH_T::get_tet_sphere_info( 
      pt0.x(), pt1.x(), pt2.x(), pt3.x(),
      pt0.y(), pt1.y(), pt2.y(), pt3.y(),
      pt0.z(), pt1.z(), pt2.z(), pt3.z(),
      x, y, z, r );

  pt0.print();
  pt1.print();
  pt2.print();
  pt3.print();

  std::cout<<cen.x() - x<<'\t';
  std::cout<<cen.y() - y<<'\t';
  std::cout<<cen.z() - z<<'\t';
  std::cout<<rad - r<<'\n';

  return EXIT_SUCCESS;
}

// EOF
