#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"

int main(int argc, char *argv[])
{
  Tensor2_3D aa {};
  SymmTensor2_3D bb {};
  
  bb.gen_rand(-3.23, 2.678);

  aa = bb.convert_to_full();

  double et1, et2, et3, st1, st2, st3;
  Vector_3 av1, av2, av3, bv1, bv2, bv3;
  aa.eigen_decomp( et1, et2, et3, av1, av2, av3 );
  bb.eigen_decomp( st1, st2, st3, bv1, bv2, bv3 );

  std::cout<<et1 - st1<<'\t';
  std::cout<<et2 - st2<<'\t';
  std::cout<<et3 - st3<<'\t';

  av1 -= bv1; av1.print();
  av2 -= bv2; av2.print();
  av3 -= bv3; av3.print();

  return EXIT_SUCCESS;
}

// EOF
