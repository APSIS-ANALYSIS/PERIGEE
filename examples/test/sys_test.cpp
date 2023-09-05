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
#include "Matrix_double_3by3_Array.hpp"
#include "Matrix_double_6by6_Array.hpp"
#include "VTK_Tools.hpp"
#include "NodalBC.hpp"
#include "DataVecStr.hpp"
#include "Tet_Tools.hpp"
#include "Hex_Tools.hpp"

int main(int argc, char *argv[])
{
  std::vector<int> nid1 {1,2,3,4,5,6,7,8};
  std::vector<int> nid2 {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27};

  std::vector<int> eid {11};

  // these node array and ien array are designed by Gmsh PDF P358, let node0 = [0, 0, 0], length of side = 2, u = x, v = y, w = z
  std::vector<double> node1 {0.0, 0.0, 0.0,
                             2.0, 0.0, 0.0,
                             2.0, 2.0, 0.0,
                             0.0, 2.0, 0.0,
                             0.0, 0.0, 2.0,
                             2.0, 0.0, 2.0,
                             2.0, 2.0, 2.0,
                             0.0, 2.0, 2.0};
  std::vector<int> ien1 {0,1,2,3,4,5,6,7};

  std::vector<double> node2 {0.0, 0.0, 0.0, // 0
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
  std::vector<int> ien2 {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};

  std::vector<DataVecStr<int>> input1 {};
  input1.push_back({nid1, "GlobalNodeID", AssociateObject::Node});
  input1.push_back({eid, "GlobalElementID", AssociateObject::Cell});
  HEX_T::write_hex_grid( "LinearHex", 8, 1, node1, ien1, input1 );

  std::vector<DataVecStr<int>> input2 {};
  input2.push_back({nid2, "GlobalNodeID", AssociateObject::Node});
  input2.push_back({eid, "GlobalElementID", AssociateObject::Cell});
  HEX_T::write_hex_grid( "QuadraticHex", 27, 1, node2, ien2, input2 ); 

  return EXIT_SUCCESS;
}

// EOF
