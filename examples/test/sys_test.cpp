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

int main(int argc, char *argv[])
{
  AssociateObject ao = AssociateObject::Cell;

  std::vector<int> b {1,2};

  DataVecStr<int> a { {1,2,3,4,5,6}, "tt", AssociateObject::Node };

  std::cout<<a<<std::endl;

  std::vector<DataVecStr<int>> mm {};

  mm.push_back({ b, "first", AssociateObject::Cell});

  VEC_T::print(mm[0].get_data());

  std::cout<<mm[0].get_name()<<std::endl;
  
  std::cout<<mm[0]<<std::endl;
  
  return EXIT_SUCCESS;
}

// EOF
