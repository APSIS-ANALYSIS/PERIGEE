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
  std::array<Vector_3, 4> pt;
  pt[0] = Vector_3(1, 1, 0);
  pt[1] = Vector_3(2, 2, 0);
  pt[2] = Vector_3(0, 1.5, 0);
  pt[3] = Vector_3(1, 1.5, 1);

  TET_T::Tet4 test_tet(pt, 2, 3, 5, 7);
  test_tet.print_info();

  std::vector<int> tet_ien {0,1,2,3};
  IIEN * IEN_v = new IEN_FEM(1, tet_ien);
  std::vector<double> ctrl_pts {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  TET_T::tetmesh_check(ctrl_pts, IEN_v, 1);

  delete IEN_v;

  std::array<Vector_3, 8> node;
  node[0] = Vector_3(0.0, 0.0, 0.0);
  node[1] = Vector_3(1.0, 0.0, 0.0);
  node[2] = Vector_3(1.0, 1.0, 0.0);
  node[3] = Vector_3(0.0, 1.0, 0.0);
  node[4] = Vector_3(0.0, 0.0, 1.0);
  node[5] = Vector_3(1.0, 0.0, 1.0);
  node[6] = Vector_3(1.0, 1.0, 1.0);
  node[7] = Vector_3(0.0, 1.0, 1.0);

  HEX_T::Hex8 test_hex(node, 11, 13, 17, 19, 23, 29, 31, 37);
  test_hex.print_info();
  std::cout << test_hex.get_aspect_ratio() << '\n' << test_hex.get_volume() << '\n' << std::endl;

  test_hex.reset(8,7,6,5,4,3,2,1);
  test_hex.print_info();
  

  return EXIT_SUCCESS;
}

// EOF
