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
#include "NodalBC_3D_vtp.hpp"
#include "NodalBC_3D_vtu.hpp"

int main(int argc, char *argv[])
{
  const int nFunc = 13500;

  std::vector<std::string> fname { "wall_vol.vtp" };

  INodalBC * a = new NodalBC_3D_vtp( fname, nFunc );

  INodalBC * b = new NodalBC( fname, nFunc );

  a -> print_info();

  b -> print_info();

  delete a; delete b;

  return EXIT_SUCCESS;
}

// EOF
