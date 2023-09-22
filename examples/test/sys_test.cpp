#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"
#include "Matrix_double_3by3_Array.hpp"
#include "Matrix_double_6by6_Array.hpp"

int main(int argc, char *argv[])
{
  MATH_T::Matrix_Dense<5> m1 {};
  m1.gen_rand(-1.2, 2.3);

  MATH_T::Matrix_dense m2(5);
  for(int ii=0; ii<5; ++ii)
    for(int jj=0; jj<5; ++jj)
      m2.set_value(ii,jj, m1(ii,jj));

  m1.LU_fac();
  m2.LU_fac();

  std::array<double,5> right, sol1, sol2;
  for(int ii=0; ii<5; ++ii) right[ii] = MATH_T::gen_double_rand(-28.6, 10.3);
  
  sol1 = m1.LU_solve(right);
  m2.LU_solve(&right[0], &sol2[0]);

  for(int ii=0; ii<5; ++ii) std::cout<<sol1[ii] - sol2[ii]<<'\t';

  return EXIT_SUCCESS;
}
// EOF
