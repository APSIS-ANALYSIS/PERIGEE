#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "SymmMatrix_3x3.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"
#include "Matrix_double_3by3_Array.hpp"
#include "Matrix_double_6by6_Array.hpp"

int main(int argc, char *argv[])
{

  const unsigned int size = 3;

  MATH_T::Matrix_Dense<size> mtest1 {};

  mtest1.gen_rand();

  mtest1.print_info();

  MATH_T::Matrix_Dense<size> mtest1_copy = mtest1;

  std::cout<<"================Matrix multiplication & LU test==============="<<std::endl;

  std::cout<<"----------------LU factorized----------------"<<std::endl;
  
  mtest1.LU_fac();

  mtest1.print_info();

  MATH_T::Matrix_Dense<size> l_mtest1 {};
  MATH_T::Matrix_Dense<size> u_mtest1 {};


  for (unsigned int ii=0; ii<size; ++ii)
  {
 
    for (unsigned int jj=0; jj<=ii; ++jj)
    {
      l_mtest1(ii, jj) =  mtest1(ii, jj);
    }

    l_mtest1(ii, ii) = 1.0;

    for (unsigned int jj=ii; jj<size; ++jj)
    {
      u_mtest1(ii, jj) =  mtest1(ii, jj);
    }
  }


  std::cout<<"----------------L----------------"<<std::endl;

  l_mtest1.print_info();


  std::cout<<"----------------U----------------"<<std::endl;

  u_mtest1.print_info();


  std::cout<<"----------------L*U----------------"<<std::endl;
  MATH_T::Matrix_Dense<size> lu_mtest1 = l_mtest1*u_mtest1;

  lu_mtest1.print_info();

  int pp[size]{};

  for (unsigned int ii = 0; ii<size; ++ii)
  {
    pp[ii] = mtest1.get_p(ii);
    //std::cout<<"p:"<<pp[ii]<<std::endl;
  }    

  for (unsigned int ii = 0; ii<size; ++ii)
  {
    for (unsigned int jj = 0; jj<size; ++jj)
    {
      std::cout << mtest1_copy(pp[ii], jj) - lu_mtest1(ii , jj) << std::endl;
    }
  }

  std::cout<<"----------------Lu solve----------------"<<std::endl;

  std::array<double, size> bb {};

  srand(time(NULL));

  for(unsigned int ii=0; ii<size; ++ii)
  {
    double value = rand() % 1000; 

    bb[ii] = value * 1.0e-2 - 5;
  }

  std::array<double, size> xx = mtest1.LU_solve(bb);

  std::array<double, size> mx = mtest1_copy*xx;

  for (unsigned ii=0; ii<size; ++ii)
    std::cout<<mx.at(ii)-bb.at(ii)<<std::endl;

  return EXIT_SUCCESS;
}
// EOF
