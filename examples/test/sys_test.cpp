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
  std::vector<int> b {1,2};

  DataVecStr<int> a( {1,2,3,4,5,6}, "tt" );

  VEC_T::print(a.data);

  std::cout<<a.name<<std::endl;
  return EXIT_SUCCESS;
}

// EOF
