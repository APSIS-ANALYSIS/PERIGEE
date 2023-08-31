#include <unistd.h>
#include "Vec_Tools.hpp"
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

int main(int argc, char *argv[])
{
  SymmMatrix_3x3 smat; smat.gen_rand();
  Matrix_3x3 mat = smat.convert_to_full();

  Tensor4_3D ten = gen_symm_id();
  SymmTensor4_3D sten = gen_ST4_symm_id();

  ten.gen_Ptilde( mat );
  sten.gen_Ptilde( smat );

  SymmMatrix_3x3 smat2; smat2.gen_rand();
  Matrix_3x3 mat2 = smat2.convert_to_full();

  ten.add_SymmOutProduct(3.14159, mat, mat2);
  sten.add_SymmOutProduct(3.14159, smat, smat2);
  

  sten.print_in_mat();

  if( sten.is_identical(ten, 1.0e-13) ) std::cout<<"passed! \n";
  else std::cout<<"error. \n";

  return EXIT_SUCCESS;
}

// EOF
