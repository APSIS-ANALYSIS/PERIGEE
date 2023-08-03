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

int main(int argc, char *argv[])
{
  Matrix_double_3by3_Array A;

  A.gen_rand();

  A.LU_fac();

  Vector_3 RHS;

  RHS.gen_rand();

  Vector_3 sol = A.LU_solve(RHS);

  auto rrhs = RHS.to_std_array();

  auto ssol = A.LU_solve(rrhs);

  double x1, x2, x3;
  A.LU_solve(RHS(0), RHS(1), RHS(2), x1, x2, x3);

  std::cout<<x1 - ssol[0]<<'\t';
  std::cout<<x2 - ssol[1]<<'\t';
  std::cout<<x3 - ssol[2]<<'\n';

  sol.print();

  return EXIT_SUCCESS;
}

// EOF
