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

int main(int argc, char *argv[])
{
  std::vector<int> nid {1,2,3,4,5,6,7,8,9};
  std::vector<int> eid {11, 21};
  std::vector<double> node {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.0, 0.5, 1.0, 0.0, 1.0, 0.5, 0.0};
  std::vector<int> ien {0,1,2,4,5,6,3,2,1,7,5,8};
  // old function cannot add Physics_tag
  TET_T::write_quadratic_triangle_grid( "old", 9, 2, node, ien, nid, eid );

  std::vector<DataVecStr<int>> input {};
  input.push_back({nid, "GlobalNodeID", AssociateObject::Node});
  input.push_back({eid, "GlobalElementID", AssociateObject::Cell});
  TET_T::write_quadratic_triangle_grid( "new", 9, 2, node, ien, input ); 

  return EXIT_SUCCESS;
}

// EOF
