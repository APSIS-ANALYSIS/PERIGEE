#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor4_3D.hpp"
#include "SymmTensor4_3D.hpp"
#include "Matrix_double_6by6_Array.hpp"

int main(int argc, char *argv[])
{
  // Matrix_double_6by6_Array TEST

  std::array<double, 9> temp9_arr {};

  for (int ii=0; ii<9; ++ii) temp9_arr.at(ii) = MATH_T::gen_double_rand(-10, 10);

  Matrix_double_6by6_Array mtest6x6 {temp9_arr.at(0), temp9_arr.at(1), temp9_arr.at(2), temp9_arr.at(3), temp9_arr.at(4), temp9_arr.at(5), temp9_arr.at(6), temp9_arr.at(7), temp9_arr.at(8)};

  MATH_T::Matrix_Dense<6> mtestD6x6 {};

  const double aa = temp9_arr.at(0); const double bb= temp9_arr.at(1); const double cc = temp9_arr.at(2);  
  const double dd = temp9_arr.at(3); const double ee= temp9_arr.at(4); const double ff = temp9_arr.at(5);  
  const double gg = temp9_arr.at(6); const double hh= temp9_arr.at(7); const double ii = temp9_arr.at(8);  

  mtestD6x6(0) = aa*aa;     mtestD6x6(1) = dd*dd;     mtestD6x6(2) = gg*gg;
  mtestD6x6(3) = 2.0*aa*dd; mtestD6x6(4) = 2.0*aa*gg; mtestD6x6(5) = 2.0*dd*gg;

  mtestD6x6(6) = aa*bb;     mtestD6x6(7)  = dd*ee;     mtestD6x6(8)  = gg*hh;
  mtestD6x6(9) = aa*ee+bb*dd; mtestD6x6(10) = aa*hh+bb*gg; mtestD6x6(11) = dd*hh+ee*gg;
  
  mtestD6x6(12) = aa*cc;     mtestD6x6(13) = dd*ff;     mtestD6x6(14) = gg*ii;
  mtestD6x6(15) = aa*ff+cc*dd; mtestD6x6(16) = aa*ii+cc*gg; mtestD6x6(17) = dd*ii+ff*gg;
  
  mtestD6x6(18) = bb*bb;     mtestD6x6(19) = ee*ee;     mtestD6x6(20) = hh*hh;
  mtestD6x6(21) = 2.0*bb*ee; mtestD6x6(22) = 2.0*bb*hh; mtestD6x6(23) = 2.0*ee*hh;

  mtestD6x6(24) = bb*cc;     mtestD6x6(25) = ee*ff;     mtestD6x6(26) = hh*ii;
  mtestD6x6(27) = bb*ff+cc*ee; mtestD6x6(28) = bb*ii+cc*hh; mtestD6x6(29) = ee*ii+ff*hh;

  mtestD6x6(30) = cc*cc;     mtestD6x6(31) = ff*ff;     mtestD6x6(32) = ii*ii;
  mtestD6x6(33) = 2.0*cc*ff; mtestD6x6(34) = 2.0*cc*ii; mtestD6x6(35) = 2.0*ff*ii;

  std::cout<<"------------The 6x6 Matrix to be solved-----------"<<std::endl;

  mtest6x6.print();  
  mtestD6x6.print_info();

  std::cout<<"------------The right-hand-side vector b6 to be solved-----------"<<std::endl;

  std::array<double, 6> b6_arr {};

  for (int ii=0; ii<6; ++ii) 
  {
    b6_arr.at(ii) = MATH_T::gen_double_rand(-10, 10);
    std::cout<<b6_arr.at(ii)<<std::endl;
  }
    
  std::cout<<"------------Matrix_double_6by6_Array vs Matrix_Dense<6>-----------"<<std::endl;   

  mtest6x6.LU_fac();

  std::array<double, 6> sol_b6_arr = mtest6x6.LU_solve(b6_arr);

  mtestD6x6.LU_fac();

  std::array<double, 6> solD_b6_arr = mtestD6x6.LU_solve(b6_arr);

  for (int ii=0; ii<6; ++ii) std::cout<<solD_b6_arr.at(ii)-sol_b6_arr[ii]<<std::endl; 
  return EXIT_SUCCESS;
}
// EOF
