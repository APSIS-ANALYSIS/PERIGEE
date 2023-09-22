#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "SymmMatrix_3x3.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"
#include "Matrix_double_3by3_Array.hpp"
#include "Matrix_double_6by6_Array.hpp"

int main(int argc, char *argv[])
{

  // ===== Initialization of PETSc =====
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  const int size = 9;

  MATH_T::Matrix_Dense<size> mtest1 {};

  mtest1.gen_rand(-10, 10);

  mtest1.print_info();

  auto mtest1_copy = mtest1;  


  std::cout<<"================template< int N> class Matrix_Dense test==============="<<std::endl;

  std::cout<<"================Matrix multiplication & LU test==============="<<std::endl;

  std::cout<<"----------------LU factorized----------------"<<std::endl;
  
  mtest1.LU_fac();

  MATH_T::Matrix_Dense<size> l_mtest1 {};
  MATH_T::Matrix_Dense<size> u_mtest1 {};

  for (int ii=0; ii<size; ++ii)
  {
 
    for (int jj=0; jj<=ii; ++jj)
    {
      l_mtest1(ii, jj) =  mtest1(ii, jj);
    }

    l_mtest1(ii, ii) = 1.0;

    for (int jj=ii; jj<size; ++jj)
    {
      u_mtest1(ii, jj) =  mtest1(ii, jj);
    }
  }

  std::cout<<"----------------L----------------"<<std::endl;

  l_mtest1.print_info();

  std::cout<<"----------------U----------------"<<std::endl;

  u_mtest1.print_info();


  std::cout<<"----------------L*U----------------"<<std::endl;
  MATH_T::Matrix_Dense<size> lu_mtest1 {}; 

  lu_mtest1.Mult(l_mtest1,u_mtest1);

  lu_mtest1.print_info();

  int pp[size]{};

  for(int ii = 0; ii<size; ++ii)
  {
    pp[ii] = mtest1.get_p(ii);
    //std::cout<<"p:"<<pp[ii]<<std::endl;
  }    

  std::cout<<"----------------Lu solve----------------"<<std::endl;

  std::array<double, size> bb {};

  //srand(time(NULL));

  for(int ii=0; ii<size; ++ii)
  {
    //double value = rand() % 1000; 

    //bb[ii] = value * 1.0e-2 - 5;
    bb[ii] = MATH_T::gen_double_rand(-10, 10);
  }

  const std::array<double, size> xx = mtest1.LU_solve(bb);

  //std::cout<< mtest1_copy.get_is_fac()<<std::endl;

  const std::array<double, size> mx = mtest1_copy.Mult(xx);

  for (int ii=0; ii<size; ++ii)
    std::cout<<mx.at(ii)-bb.at(ii)<<std::endl;


  std::cout<<"================template< int N> class Matrix_SymPos_Dense test==============="<<std::endl;

  MATH_T::Matrix_Dense<size> mtest1_copy_T = mtest1_copy;
  mtest1_copy_T.transpose();

  //mtest1_copy_T.print_info();

  //mtest1_copy.print_info();

  for(int ii=0; ii<size; ++ii)
  {
    for(int jj=0; jj<size; ++jj) std::cout<<mtest1_copy_T (jj, ii) - mtest1_copy(ii, jj)<<std::endl;         
  }

  for(int ii=0; ii<size; ++ii) std::cout<<mtest1_copy.get_p(ii) - mtest1_copy_T.get_p(ii)<<std::endl;
  for(int ii=0; ii<size; ++ii) std::cout<<mtest1_copy.get_is_fac()<<" "<<mtest1_copy_T.get_is_fac()<<std::endl;

  //(mtest1_copy_T.Mult(mtest1_copy)).print_info();

  MATH_T::Matrix_SymPos_Dense<size> symtest1{}; 

  symtest1.Mult(mtest1_copy_T, mtest1_copy); 

  MATH_T::Matrix_Dense<size> temp_mtest1{}; 

  temp_mtest1.Mult(mtest1_copy_T, mtest1_copy);

  const MATH_T::Matrix_SymPos_Dense<size> temp_symtest1{temp_mtest1}; 

  for(int ii=0; ii<size*size; ++ii) std::cout<<symtest1(ii) - temp_symtest1(ii)<<std::endl;
  for(int ii=0; ii<size; ++ii) std::cout<<symtest1.get_p(ii) - temp_symtest1.get_p(ii)<<std::endl;
  for(int ii=0; ii<size; ++ii) std::cout<<symtest1.get_is_fac()<<" "<<temp_symtest1.get_is_fac()<<std::endl;

  symtest1.print_info();

  symtest1.check_symm();

  const MATH_T::Matrix_SymPos_Dense<size> symtest1_copy = symtest1;

  for(int ii=0; ii<size*size; ++ii) std::cout<<symtest1_copy(ii) - symtest1(ii)<<std::endl;
  for(int ii=0; ii<size; ++ii) std::cout<<symtest1_copy.get_p(ii) - symtest1.get_p(ii)<<std::endl;
  for(int ii=0; ii<size; ++ii) std::cout<<symtest1_copy.get_is_fac()<<" "<<symtest1.get_is_fac()<<std::endl;

/*
  MATH_T::Matrix_SymPos_Dense<size> * psymtest1 = new  MATH_T::Matrix_SymPos_Dense<size>();

  psymtest1 -> gen_rand();

  MATH_T::Matrix_SymPos_Dense<size> * psymtest2 = new  MATH_T::Matrix_SymPos_Dense<size>();

  //psymtest2-> copy(psymtest1);

  *psymtest2 = *psymtest1;

  for ( int ii = 0; ii<size; ++ii)
  {
    for ( int jj = 0; jj<size; ++jj)
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

  for(int ii=0; ii<size; ++ii)
  {
    for(int jj=0; jj<=ii; ++jj)
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
  
  MATH_T::Matrix_Dense<size> ld_symtest1 {};  

  ld_symtest1.Mult(l_symtest1, d_symtest1);

  MATH_T::Matrix_SymPos_Dense<size> ldlt_symtest1 {};

  auto temp3 = l_symtest1; temp3.transpose();
  ldlt_symtest1.Mult(ld_symtest1, temp3); 

  ldlt_symtest1.print_info();

  for(int ii=0; ii<size * size; ++ii)
    std::cout << ldlt_symtest1(ii) - symtest1_copy(ii) << std::endl; 

  std::cout<<"----------------LDLt solve----------------"<<std::endl;

  const std::array<double, size> xx2 = symtest1.LDLt_solve(bb);
  const std::array<double, size> mx2 = symtest1_copy.Mult(xx2);

  for(int ii=0; ii<size; ++ii)
    std::cout<<mx2.at(ii) - bb.at(ii)<<std::endl;


  std::cout<<"================other test==============="<<std::endl;

  std::array<double, size*size> array_test1 {};

  for(int ii=0; ii<size*size; ++ii)
  {
    //double value = rand() % 1000; 

    //array_test1[ii] = value * 1.0e-2 - 5;
    array_test1[ii] = MATH_T::gen_double_rand(-10, 10);

  }

  MATH_T::Matrix_Dense<size> array_mtest1 {array_test1};
  MATH_T::Matrix_SymPos_Dense<size> array_symtest1 {array_test1};

  //array_mtest1.print_info();
  //array_symtest1.print_info();

  for(int ii=0; ii<size*size; ++ii)
  {
    std::cout << array_mtest1(ii) - array_test1[ii] << std::endl;
    std::cout << array_symtest1(ii) - array_test1[ii] << std::endl;
  }

  for(int ii=0; ii<size; ++ii) std::cout<<array_mtest1.get_p(ii) - array_symtest1.get_p(ii)<<std::endl;
  for(int ii=0; ii<size; ++ii) std::cout<<array_mtest1.get_is_fac()<<" "<<array_symtest1.get_is_fac()<<std::endl;

  std::cout<<"================Comparison 1==============="<<std::endl;
  SYS_T::Timer * mytimer = new SYS_T::Timer();

  // Matrix_double_3by3_Array vs Matrix_Dense<3>

  double totalTime_3x3 = 0;
  double totalTimeD_3x3 = 0;

  for(int nn = 0; nn<10000; ++nn)
  {
    Matrix_double_3by3_Array mtest1_3x3 {};
 
    mtest1_3x3.gen_rand(-10, 10);

    //the pivot is also required to initialize
    MATH_T::Matrix_Dense<3> mtest1d_3x3 {};

    for(int ii=0; ii<9; ++ii)
    { 
      mtest1d_3x3(ii) = mtest1_3x3(ii);
    }

    //mtest1d_3x3.print_info();
    //mtest1_3x3.print();

    std::array<double, 3> b3 {};

    for(int ii=0; ii<3; ++ii)
    {
      b3[ii] = MATH_T::gen_double_rand(-10, 10);
    }
/*
    for( int ii=0; ii<3; ++ii)
    {
      std::cout<<b3.at(ii)<< std::endl;
    }
*/
    mytimer->Reset();
    mytimer->Start();

    mtest1_3x3.LU_fac();
    //mtest1_3x3.print();

    std::array<double, 3> sol3 = mtest1_3x3.LU_solve(b3);

    mytimer -> Stop();

    double Time_3x3 = mytimer->get_sec();

    totalTime_3x3 += Time_3x3;

    mytimer->Reset();
    mytimer->Start();

    mtest1d_3x3.LU_fac();
    //mtest1d_3x3.print_info();
    std::array<double, 3> sol3d = mtest1d_3x3.LU_solve(b3);

    mytimer -> Stop();

    double TimeD_3x3 = mytimer->get_sec();

    totalTimeD_3x3 += TimeD_3x3;
    //std::cout<<"sol3 - sol3d \n"; 
    //for( int ii=0; ii<3; ++ii)
   // {
   //   std::cout<<sol3[ii] -  sol3d[ii] << '\t';
   // }
   // std::cout<<'\n';
    
    for(int ii=0; ii<3; ++ii)
    {
      if(std::abs(sol3.at(ii) - sol3d.at(ii)) > 1.0e-10)
        std::cout<<"nn="<< nn<<":diff between Matrix_double_3by3_Array and Matrix_Dense<3>:"<< sol3.at(ii) - sol3d.at(ii) << std::endl;
    }
  }

  std::cout << "Matrix_double_3by3_Array LU total runtime：" << totalTime_3x3 << " s" << std::endl;
  std::cout << "Matrix_Dense<3> LU total runtime：         " << totalTimeD_3x3 << " s" << std::endl;

  std::cout<<"---------------------------------"<<std::endl;
  // Matrix_double_6by6_Array vs Matrix_Dense<6>

  double totalTime_6x6 = 0;
  double totalTimeD_6x6 = 0;

  for(int nn = 0; nn<10000; ++nn)
  {
    MATH_T::Matrix_Dense<6> mtest1d_6x6 {};
    
    mtest1d_6x6.gen_rand(-10, 10);

    double atest1_6x6[36] = {0.0} ; 

    for(int ii=0; ii<36; ++ii)
    { 
      atest1_6x6[ii] = mtest1d_6x6(ii);
    }

    Matrix_double_6by6_Array mtest1_6x6 (atest1_6x6);

    //mtest1d_6x6.print_info();
    //mtest1_6x6.print();

    std::array<double, 6> b6 {};

    for(int ii=0; ii<6; ++ii)
    {
      //double value = rand() % 1000; 

      //b6[ii] = value * 1.0e-2 - 5;
      b6[ii] = MATH_T::gen_double_rand(-10, 10);
    }
/*
    for( int ii=0; ii<6; ++ii)
    {
      std::cout<<b6.at(ii)<< std::endl;
    }
*/
    mytimer->Reset();
    mytimer->Start();

    mtest1_6x6.LU_fac();
    //mtest1_6x6.print();
    const std::array<double, 6> sol6 = mtest1_6x6.LU_solve(b6);

    mytimer -> Stop();

    double Time_6x6 = mytimer->get_sec();

    totalTime_6x6 += Time_6x6;

    mytimer->Reset();
    mytimer->Start();

    mtest1d_6x6.LU_fac();
    //mtest1d_6x6.print_info();
    const std::array<double, 6> sol6d = mtest1d_6x6.LU_solve(b6);

    mytimer -> Stop();

    double TimeD_6x6 = mytimer->get_sec();

    totalTimeD_6x6 += TimeD_6x6;
/*
    for( int ii=0; ii<6; ++ii)
    {
      std::cout<<sol6[ii]<<" "<< sol6d[ii] << std::endl;
    }
*/
    for(int ii=0; ii<6; ++ii)
    {
      if(std::abs(sol6.at(ii) - sol6d.at(ii)) > 1.0e-10)
        std::cout<<"nn="<< nn<<":diff between Matrix_double_6by6_Array and Matrix_Dense<6>:"<< sol6.at(ii) - sol6d.at(ii) << std::endl;
    }
  }

  std::cout << "Matrix_double_6by6_Array LU total runtime：" << totalTime_6x6 << " s" << std::endl;
  std::cout << "Matrix_Dense<6> LU total runtime：         " << totalTimeD_6x6 << " s" << std::endl;


  std::cout<<"================Comparison 2==============="<<std::endl;
  // Matrix_double_3by3_Array vs Matrix_SymPos_Dense<3>

  totalTime_3x3 = 0;
  totalTimeD_3x3 = 0;

  std::cout<<"error norm diff between Matrix_double_3by3_Array and Matrix_SymPos_Dense<3>:"<<std::endl;

  for(int nn = 0; nn<10000; ++nn)
  {
    MATH_T::Matrix_Dense<3> mtest1d_3x3 {};

    mtest1d_3x3.gen_rand(-10, 10);

    MATH_T::Matrix_SymPos_Dense<3> smtest1d_3x3{};

    auto temp1 = mtest1d_3x3; temp1.transpose();

    smtest1d_3x3.Mult(temp1, mtest1d_3x3);

    const MATH_T::Matrix_SymPos_Dense<3> smtest1d_3x3_copy = smtest1d_3x3;

    smtest1d_3x3.check_symm();

    smtest1d_3x3_copy.check_symm();

    Matrix_double_3by3_Array mtest1_3x3 {};

    for(int ii=0; ii<9; ++ii)
    { 
      mtest1_3x3(ii) = smtest1d_3x3(ii);
    }

    Matrix_double_3by3_Array mtest1_3x3_copy = mtest1_3x3;

    //smtest1d_3x3.print_info();
    //mtest1_3x3.print();

    std::array<double, 3> b3 {};

    for(int ii=0; ii<3; ++ii)
    {
      //double value = rand() % 1000; 

      //b3[ii] = value * 1.0e-2 - 5;
      b3[ii] = MATH_T::gen_double_rand(-10, 10);
    }
/*
    for( int ii=0; ii<3; ++ii)
    {
      std::cout<<b3.at(ii)<< std::endl;
    }
*/
    mytimer->Reset();
    mytimer->Start();

    mtest1_3x3.LU_fac();
    //mtest1_3x3.print();

    const std::array<double, 3> sol3 = mtest1_3x3.LU_solve(b3);

    mytimer -> Stop();

    double Time_3x3 = mytimer->get_sec();

    totalTime_3x3 += Time_3x3;

    // compute norm of error
    const double sol3_arr[3] {sol3[0], sol3[1], sol3[2]};

    double act_sol3[3] = {0.0};

    mtest1_3x3_copy.VecMult(sol3_arr ,act_sol3);

    double error_3[3] {}; 

    for(int ii=0; ii<3; ++ii) error_3[ii] = act_sol3[ii]-b3[ii];

    double sum_square = 0.0;

    for(int ii=0; ii<3; ++ii)
    {
      sum_square += (act_sol3[ii]-b3[ii])*(act_sol3[ii]-b3[ii]);
    }

    const double norm2_3 = std::sqrt(sum_square);

    mytimer->Reset();
    mytimer->Start();

    smtest1d_3x3.LDLt_fac();
    //smtest1d_3x3.print_info();
    const std::array<double, 3> sol3d = smtest1d_3x3.LDLt_solve(b3);

    mytimer -> Stop();

    double TimeD_3x3 = mytimer->get_sec();

    totalTimeD_3x3 += TimeD_3x3;
/*
    for( int ii=0; ii<3; ++ii)
    {
      std::cout<<sol3[ii]<<" "<< sol3d[ii] << std::endl;
    }
*/

    // compute norm of error     
    const std::array<double, 3> act_sol3d = smtest1d_3x3_copy.Mult(sol3d);

    double error_3d[3] {}; 

    for(int ii=0; ii<3; ++ii) error_3d[ii] = act_sol3d[ii]-b3[ii];

    sum_square = 0.0;

    for(int ii=0; ii<3; ++ii)
    {
      sum_square += (act_sol3d[ii]-b3[ii])*(act_sol3d[ii]-b3[ii]);
    }

    const double norm2_3d = std::sqrt(sum_square);

    if(std::abs(norm2_3- norm2_3d) > 1.0e-10)
    {
      std::cout<<"nn="<< nn<<": "<< norm2_3<<" "<<norm2_3d<< std::endl;
        
      for(int ii=0; ii<3; ++ii) std::cout<<"             "<< error_3[ii]<<" "<<error_3d[ii]<< std::endl;
    }

  }

  std::cout << "Matrix_double_3by3_Array LU total runtime：" << totalTime_3x3 << " s" << std::endl;
  std::cout << "Matrix_SymPos_Dense<3> LDLt total runtime：" << totalTimeD_3x3 << " s" << std::endl;


  std::cout<<"---------------------------------"<<std::endl;
  // Matrix_double_6by6_Array vs Matrix_SymPos_Dense<6>

  totalTime_6x6 = 0;
  totalTimeD_6x6 = 0;
    
  std::cout<<"error norm diff between Matrix_double_6by6_Array and Matrix_SymPos_Dense<6>:"<<std::endl;

  for(int nn = 0; nn<10000; ++nn)
  {
    MATH_T::Matrix_Dense<6> mtest1d_6x6 {};
    
    mtest1d_6x6.gen_rand(-10, 10);

    MATH_T::Matrix_SymPos_Dense<6> smtest1d_6x6{};

    auto temp2 = mtest1d_6x6; temp2.transpose();

    smtest1d_6x6.Mult(temp2, mtest1d_6x6);

    const MATH_T::Matrix_SymPos_Dense<6> smtest1d_6x6_copy = smtest1d_6x6;

    smtest1d_6x6.check_symm();

    smtest1d_6x6_copy.check_symm();

    double atest1_6x6[36] = {0.0} ; 

    for(int ii=0; ii<36; ++ii)
    { 
      atest1_6x6[ii] = smtest1d_6x6(ii);
    }

    Matrix_double_6by6_Array mtest1_6x6 (atest1_6x6);

    //mtest1d_6x6.print_info();
    //mtest1_6x6.print();

    std::array<double, 6> b6 {};

    for(int ii=0; ii<6; ++ii)
    {
      //double value = rand() % 1000; 

      //b6[ii] = value * 1.0e-2 - 5;
      b6[ii] = MATH_T::gen_double_rand(-10, 10);
    }
/*
    for( int ii=0; ii<6; ++ii)
    {
      std::cout<<b6.at(ii)<< std::endl;
    }
*/
    mytimer->Reset();
    mytimer->Start();

    mtest1_6x6.LU_fac();
    //mtest1_6x6.print();

    const std::array<double, 6> sol6 = mtest1_6x6.LU_solve(b6);

    mytimer -> Stop();

    double Time_6x6 = mytimer->get_sec();

    totalTime_6x6 += Time_6x6;

    // compute norm of error
    const double sol6_arr[6] {sol6[0], sol6[1], sol6[2], sol6[3], sol6[4], sol6[5]};

    double act_sol6[6] = {0.0};

    for(int ii=0; ii<6; ++ii)
    {
      for(int jj=0; jj<6; ++jj)
      act_sol6[ii] += atest1_6x6[6*ii+jj] * sol6_arr[jj];
    }

    double error_6[6] {}; 

    for(int ii=0; ii<6; ++ii) error_6[ii] = act_sol6[ii]-b6[ii];

    double sum_square = 0.0;

    for(int ii=0; ii<6; ++ii)
    {
      sum_square += (act_sol6[ii]-b6[ii])*(act_sol6[ii]-b6[ii]);
    }

    const double norm2_6 = std::sqrt(sum_square);

    mytimer->Reset();
    mytimer->Start();

    smtest1d_6x6.LDLt_fac();
    //mtest1d_6x6.print_info();
    const std::array<double, 6> sol6d = smtest1d_6x6.LDLt_solve(b6);

    mytimer -> Stop();

    double TimeD_6x6 = mytimer->get_sec();

    totalTimeD_6x6 += TimeD_6x6;
    
    std::cout<<"sol6 - sol6d : \t";
    for( int ii=0; ii<6; ++ii)
    {
      std::cout<<sol6[ii] - sol6d[ii] << '\t';
    }
    std::cout<<'\n';

    // compute error     
    const std::array<double, 6> act_sol6d = smtest1d_6x6_copy.Mult(sol6d);
     
    double error_6d[6] {}; 

    for(int ii=0; ii<6; ++ii) error_6d[ii] = act_sol6d[ii]-b6[ii];

    sum_square = 0.0;

    for(int ii=0; ii<6; ++ii)
    {
      sum_square += (act_sol6d[ii]-b6[ii])*(act_sol6d[ii]-b6[ii]);
    }

    const double norm2_6d = std::sqrt(sum_square);

    if(std::abs(norm2_6 - norm2_6d) > 1.0e-10)
    {
      std::cout<<"nn="<< nn<<": "<< norm2_6<<" "<<norm2_6d<< std::endl;
        
      for(int ii=0; ii<6; ++ii) std::cout<<"             "<< error_6[ii]<<" "<<error_6d[ii]<< std::endl;
    }

  }

  std::cout << "Matrix_double_6by6_Array LU total runtime：" << totalTime_6x6 << " s" << std::endl;
  std::cout << "Matrix_SymPos_Dense<6> LDLt total runtime：" << totalTimeD_6x6 << " s" << std::endl;

  PetscFinalize();

  return EXIT_SUCCESS;
}
// EOF
