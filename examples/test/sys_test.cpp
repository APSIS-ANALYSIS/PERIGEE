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
  std::vector<int> eid {11, 21};
  std::vector<int> ptag {-2, 23};
  TET_T::write_tet_grid( "old", 5, 2, {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0}, {0,1,2,3,1, 2, 3,4}, {-1, 32, 5, 3, 22}, eid, ptag, true );

  std::vector<DataVecStr<int>> input {};
  input.push_back({{-1,32,5,3, 22}, "GlobalNID", AssociateObject::Node});
  input.push_back({eid, "GlobalEID", AssociateObject::Cell});
  input.push_back({{-2, 23}, "Physics_tag", AssociateObject::Cell});

  TET_T::write_tet_grid( "new", 5, 2, {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0}, {0,1,2,3, 1,2,3,4}, {} ); 

  return EXIT_SUCCESS;
}

// EOF
