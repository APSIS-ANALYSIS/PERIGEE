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
  std::vector<int> nid1 {1,2,3,4};

  std::vector<int> eid {11};

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

  std::vector<int> nid2 {1,2,3,4,5,6,7,8,9};

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

  return EXIT_SUCCESS;
}

// EOF
