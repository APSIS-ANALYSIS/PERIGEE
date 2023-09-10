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
  std::vector<int> nid1 {0, 1, 2, 3};

  std::vector<int> eid {0};

  // these node array and ien array are designed by Gmsh PDF P358, let node0 = [0, 0, 0], length of side = 2, u = x, v = y, w = z
  std::vector<double> node1 {0.0, 0.0, 0.0,
                             2.0, 0.0, 0.0,
                             2.0, 2.0, 0.0,
                             0.0, 2.0, 0.0};
  std::vector<int> ien1 {0,1,2,3};

  std::vector<DataVecStr<int>> input1 {};
  input1.push_back({nid1, "GlobalNodeID", AssociateObject::Node});
  input1.push_back({eid, "GlobalElementID", AssociateObject::Cell});
  HEX_T::write_quadrangle_grid( "LinearQuad", 4, 1, node1, ien1, input1 );

  std::vector<double> vol_pts1 {0.0, 0.0, 0.0,
                                2.0, 0.0, 0.0,
                                2.0, 2.0, 0.0,
                                0.0, 2.0, 0.0,
                                0.0, 0.0, 2.0,
                                2.0, 0.0, 2.0,
                                2.0, 2.0, 2.0,
                                0.0, 2.0, 2.0};
  std::vector<int> vol_ien1 {0,1,2,3,4,5,6,7};
  IIEN * IEN_v1 = new IEN_FEM(1, vol_ien1);

  const Vector_3 out1 = HEX_T::get_out_normal("LinearQuad.vtp", vol_pts1, IEN_v1);
  out1.print();

  delete IEN_v1; 

  std::vector<int> nid2 {0, 1, 2, 3, 4, 5, 6, 7, 8};

  std::vector<double> node2 {0.0, 0.0, 0.0,
                             2.0, 0.0, 0.0,
                             2.0, 2.0, 0.0,
                             0.0, 2.0, 0.0,
                             1.0, 0.0, 0.0,
                             2.0, 1.0, 0.0,
                             1.0, 2.0, 0.0,
                             0.0, 1.0, 0.0,
                             1.0, 1.0, 0.0};
  std::vector<int> ien2 {0,1,2,3,4,5,6,7,8};

  std::vector<DataVecStr<int>> input2 {};
  input2.push_back({nid2, "GlobalNodeID", AssociateObject::Node});
  input2.push_back({eid, "GlobalElementID", AssociateObject::Cell});
  HEX_T::write_quadratic_quadrangle_grid( "QuadraticQuad", 9, 1, node2, ien2, input2 );

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
  std::vector<int> vol_ien2 {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
  IIEN * IEN_v2 = new IEN_FEM(1, vol_ien2);

  Vector_3 out2 = HEX_T::get_out_normal("QuadraticQuad.vtu", vol_pts2, IEN_v2);
  out2.print();

  delete IEN_v2;

  return EXIT_SUCCESS;
}

// EOF
