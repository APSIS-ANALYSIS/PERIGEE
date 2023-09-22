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
  MATH_T::Matrix_Dense<7> aa {};
  aa.gen_rand(-5.2, 7.3);

  auto bb = aa; bb.transpose();

  auto cc = bb;
  cc.Mult(aa, bb);

  MATH_T::Matrix_SymPos_Dense<7> m1 {};
  MATH_T::Matrix_SymPos_dense m2(7);
  for(int ii=0; ii<7; ++ii)
  {
    for(int jj=0; jj<7; ++jj)
    {
      m1(ii, jj) = cc(ii, jj);
      m2.set_value(ii,jj, cc(ii,jj));
    }
  }
  
  m1.print_info();
  m2.print_info();
  
  m1.LDLt_fac();
  m2.LDLt_fac();

  std::array<double,7> right, sol1, sol2;
  for(int ii=0; ii<7; ++ii) right[ii] = MATH_T::gen_double_rand(-28.6, 70.3);

  sol1 = m1.LDLt_solve(right);
  m2.LDLt_solve(&right[0], &sol2[0]);

  for(int ii=0; ii<7; ++ii) std::cout<<sol1[ii]<<'\t'<<sol1[ii] - sol2[ii]<<'\t';

  return EXIT_SUCCESS;
}
// EOF
