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

int main(int argc, char *argv[])
{
  std::vector<double> vol_pts1 {0.0, 0.0, 0.0,
                                1.0, 0.0, 0.0,
                                1.0, 1.0, 0.0,
                                0.0, 1.0, 0.0,
                                0.0, 0.0, 1.0,
                                1.0, 0.0, 1.0,
                                1.0, 1.0, 1.0,
                                0.0, 1.0, 1.0};
  std::vector<int> vol_ien1 {0,1,2,3,4,5,6,7};
  IIEN * IEN_v1 = new IEN_FEM(1, vol_ien1);

  std::vector<double> vol_pts2 {0.0, 0.0, 0.0, // 0
                                2.0, 0.0, 0.0, // 1
                                2.0, 2.0, 0.0, // 2
                                0.0, 2.0, 0.0, // 3
                                0.0, 0.0, 2.0, // 4
                                2.0, 0.0, 2.0, // 5
                                2.0, 2.0, 2.0, // 6
                                0.0, 2.0, 2.0, // 7
                                1.0, 0.0, 0.0, // 8
                                0.0, 1.0, 0.0, // 9
                                0.0, 0.0, 1.0, // 10
                                2.0, 1.0, 0.0, // 11
                                2.0, 0.0, 1.0, // 12
                                1.0, 2.0, 0.0, // 13
                                2.0, 2.0, 1.0, // 14
                                0.0, 2.0, 1.0, // 15
                                1.0, 0.0, 2.0, // 16
                                0.0, 1.0, 2.0, // 17
                                2.0, 1.0, 2.0, // 18
                                1.0, 2.0, 2.0, // 19
                                1.0, 1.0, 0.0, // 20
                                1.0, 0.0, 1.0, // 21
                                0.0, 1.0, 1.0, // 22
                                2.0, 1.0, 1.0, // 23
                                1.0, 2.0, 1.0, // 24
                                1.0, 1.0, 2.0, // 25
                                1.0, 1.0, 1.0};// 26

  // test for Hex8
  HEX_T::Hex8 test_hex;
  test_hex.print_info();
  std::cout << test_hex.get_aspect_ratio() << '\n' << test_hex.get_volume() << '\n' << std::endl;

  test_hex.reset(vol_pts1, IEN_v1, 0);
  test_hex.print_info();
  std::cout << test_hex.get_aspect_ratio() << '\n' << test_hex.get_volume() << '\n' << std::endl;
  delete IEN_v1;

  test_hex.reset(2, 3, 5, 7, 9, 11, 13, 17);
  test_hex.print_info();
  std::cout << test_hex.get_aspect_ratio() << '\n' << test_hex.get_volume() << '\n' << std::endl;

  test_hex.reset(vol_pts2, 0, 1, 2, 3, 4, 5 ,6, 7);
  test_hex.print_info();
  std::cout << test_hex.get_aspect_ratio() << '\n' << test_hex.get_volume() << '\n' << std::endl;

  return EXIT_SUCCESS;
}

// EOF
