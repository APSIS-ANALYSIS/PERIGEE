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
  std::vector<int> nid {32, 5, 3, 22};
  std::vector<int> eid {11, 21};
  // old function cannot add Physics_tag
  TET_T::write_triangle_grid( "old", 4, 2, {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0}, {0,1,2,2,1,3}, nid, eid );

  std::vector<DataVecStr<int>> input {};
  input.push_back({nid, "GlobalNodeID", AssociateObject::Node});
  input.push_back({eid, "GlobalElementID", AssociateObject::Cell});
  input.push_back({{-2, 23}, "Physics_tag", AssociateObject::Cell});

  TET_T::write_triangle_grid( "new", 4, 2, {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0}, {0,1,2,2,1,3}, input ); 

  return EXIT_SUCCESS;
}

// EOF
