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

int main(int argc, char *argv[])
{
  Matrix_double_6by6_Array A(-1.1, -2.3, 3.5, 4.8, 5.5, -6.0, 7.023, -18.0, 9.0);

  A.LU_fac();

  double rhs [] = {1.0, 222222.0, -332325.2, 4.0, 5.5, 6.9};

  double sol [] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  A.LU_solve(rhs, sol);

  std::array<double, 6> rrhs;

  for(int ii=0; ii<6; ++ii) rrhs[ii] = rhs[ii];

  auto ssol = A.LU_solve(rrhs);

  std::cout<<sol[0] - ssol[0]<<'\t';
  std::cout<<sol[1] - ssol[1]<<'\t';
  std::cout<<sol[2] - ssol[2]<<'\t';
  std::cout<<sol[3] - ssol[3]<<'\t';
  std::cout<<sol[4] - ssol[4]<<'\t';
  std::cout<<sol[5] - ssol[5]<<'\n';

  for(int ii=0; ii<6; ++ii) std::cout<<sol[ii]<<'\n';

  return EXIT_SUCCESS;
}

// EOF
