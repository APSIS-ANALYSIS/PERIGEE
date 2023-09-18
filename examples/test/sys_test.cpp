#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor2_3D.hpp"
#include "SymmTensor2_3D.hpp"
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
  Vector_3 p1, p2, p3, p4;
  p1.gen_rand(-2.0, 3.3);
  p2.gen_rand(-2.1, 6.3);
  p3.gen_rand(-112.0, -23.3);
  p4.gen_rand(-2.3, 7.3);

  p1.print();
  p2.print();
  p3.print();
  p4.print();

  double cr;
  const auto cc = MATH_T::get_tet_sphere_info(p1, p2, p3, p4, cr );

  std::cout<<cr - MATH_T::get_circumradius( std::array<Vector_3, 4> {{p1, p2, p3, p4}})<<'\n';
  return EXIT_SUCCESS;
}

// EOF
