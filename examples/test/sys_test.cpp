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

  std::cout<<"================template<unsigned int N> class Matrix_Dense test==============="<<std::endl;

  std::cout<<"================Matrix multiplication & LU test==============="<<std::endl;

  std::cout<<"----------------LU factorized----------------"<<std::endl;
  
  mtest1.LU_fac();

  mtest1.print_info();

  //MATH_T::Matrix_Dense<size> mtest11 = mtest1;

  //mtest11.print_info();

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


  std::cout<<"================template<unsigned int N> class Matrix_SymPos_Dense test==============="<<std::endl;

  MATH_T::Matrix_Dense<size> mtest1_copy_T = transpose(mtest1_copy);

  //mtest1_copy_T.print_info();

  //mtest1_copy.print_info();

  for(unsigned int ii=0; ii<size; ++ii)
  {
    for(unsigned int jj=0; jj<size; ++jj) std::cout<<mtest1_copy_T (jj, ii) - mtest1_copy(ii, jj)<<std::endl;         
  }

  //(mtest1_copy_T * mtest1_copy).print_info();

  MATH_T::Matrix_SymPos_Dense<size> symtest1(mtest1_copy_T * mtest1_copy); 

  symtest1.print_info();

  symtest1.check_symm();

  MATH_T::Matrix_SymPos_Dense<size> symtest1_copy = symtest1;

/*
  MATH_T::Matrix_SymPos_Dense<size> * psymtest1 = new  MATH_T::Matrix_SymPos_Dense<size>();

  psymtest1 -> gen_rand();

  MATH_T::Matrix_SymPos_Dense<size> * psymtest2 = new  MATH_T::Matrix_SymPos_Dense<size>();

  //psymtest2-> copy(psymtest1);

  *psymtest2 = *psymtest1;

  for (unsigned int ii = 0; ii<size; ++ii)
  {
    for (unsigned int jj = 0; jj<size; ++jj)
    {
      std::cout<<psymtest1 -> operator()(ii, jj) - ( psymtest2 -> operator()(ii, jj) )<<std::endl;
    }
  }
*/ 
  std::cout<<"----------------LDLt----------------"<<std::endl;

  symtest1.LDLt_fac();

  symtest1.print_info();

  MATH_T::Matrix_Dense<size> l_symtest1 {};
  MATH_T::Matrix_Dense<size> d_symtest1 {};

  for (unsigned int ii=0; ii<size; ++ii)
  {
    for (unsigned int jj=0; jj<=ii; ++jj)
    {
      l_symtest1(ii, jj) =  symtest1(ii, jj);
    }

    l_symtest1(ii, ii) = 1.0;

    d_symtest1(ii, ii) =  symtest1(ii, ii);
  }

  std::cout<<"----------------L----------------"<<std::endl;

  l_symtest1.print_info();

  std::cout<<"----------------D----------------"<<std::endl;
  d_symtest1.print_info();    

  std::cout<<"----------------L*D*L^t----------------"<<std::endl;
    
  MATH_T::Matrix_Dense<size> ldlt_symtest1 = (l_symtest1* d_symtest1) * transpose(l_symtest1); 

  ldlt_symtest1.print_info();

  for (unsigned int ii=0; ii<size * size; ++ii)
    std::cout << ldlt_symtest1(ii) - symtest1_copy(ii) << std::endl; 

  std::cout<<"----------------LDLt solve----------------"<<std::endl;

  std::array<double, size> xx2 = symtest1.LDLt_solve(bb);
  std::array<double, size> mx2 = symtest1_copy*xx2;

  for (unsigned ii=0; ii<size; ++ii)
    std::cout<<mx2.at(ii) - bb.at(ii)<<std::endl;


  std::cout<<"================other test==============="<<std::endl;

  std::array<double, size*size> array_test1 {};

  srand(time(NULL));

  for(unsigned int ii=0; ii<size*size; ++ii)
  {
    double value = rand() % 1000; 

    array_test1[ii] = value * 1.0e-2 - 5;
  }

  MATH_T::Matrix_Dense<size> array_mtest1 {array_test1};
  MATH_T::Matrix_Dense<size> array_symtest1 {array_test1};

  //array_mtest1.print_info();
  //array_symtest1.print_info();

  for(unsigned int ii=0; ii<size*size; ++ii)
  {
    std::cout << array_mtest1(ii) - array_test1[ii] << std::endl;
    std::cout << array_symtest1(ii) - array_test1[ii] << std::endl;
  }


  return EXIT_SUCCESS;
}
// EOF
