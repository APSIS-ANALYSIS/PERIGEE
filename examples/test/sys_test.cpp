#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "SymmMatrix_3x3.hpp"
#include "HDF5_Reader.hpp"
#include "PostVectSolution.hpp"

int main( int argc, char * argv[] )
{
  //std::vector<std::string> a;
  //a.push_back("hello a\n");
  //a.push_back("cmame \t");
  //a.push_back("july 31th \n");
  //a.push_back("aug 31th \n");

  //for(auto &out : a) out = "hello\n";

  //for(auto out : a) std::cout<<out;

  SymmMatrix_3x3 A;

  A.gen_rand();

  Matrix_3x3 B( A.xx(), A.xy(), A.xz(), A.yx(), A.yy(), A.yz(), A.zx(), A.zy(), A.zz() );

  std::cout<<A.I1() - B.I1()<<std::endl;
  std::cout<<A.I2() - B.I2()<<std::endl;
  std::cout<<A.I3() - B.I3()<<std::endl;
  std::cout<<A.det() - B.det()<<std::endl;
  std::cout<<A.tr() - B.tr()<<std::endl;

  A.print_in_row();

  A.print();

  A.print_Voigt();

  sleep(1);
  SymmMatrix_3x3 A1;
  A1.gen_rand();
  Matrix_3x3 B1( A1.xx(), A1.xy(), A1.xz(), A1.yx(), A1.yy(), A1.yz(), A1.zx(), A1.zy(), A1.zz() );

  std::cout<<(A1 + A).xx() - (B1 + B).xx()<<std::endl;
  std::cout<<(A1 + A).xy() - (B1 + B).xy()<<std::endl;
  std::cout<<(A1 + A).xz() - (B1 + B).xz()<<std::endl;
  std::cout<<(A1 + A).yy() - (B1 + B).yy()<<std::endl;
  std::cout<<(A1 + A).yz() - (B1 + B).yz()<<std::endl;
  std::cout<<(A1 + A).zz() - (B1 + B).zz()<<std::endl;
  std::cout<<(A1 + A).yx() - (B1 + B).yx()<<std::endl;
  std::cout<<(A1 + A).zx() - (B1 + B).zx()<<std::endl;
  std::cout<<(A1 + A).zy() - (B1 + B).zy()<<std::endl;

  std::cout<<(A1 - A).xx() - (B1 - B).xx()<<std::endl;
  std::cout<<(A1 - A).xy() - (B1 - B).xy()<<std::endl;
  std::cout<<(A1 - A).xz() - (B1 - B).xz()<<std::endl;
  std::cout<<(A1 - A).yy() - (B1 - B).yy()<<std::endl;
  std::cout<<(A1 - A).yz() - (B1 - B).yz()<<std::endl;
  std::cout<<(A1 - A).zz() - (B1 - B).zz()<<std::endl;
  std::cout<<(A1 - A).yx() - (B1 - B).yx()<<std::endl;
  std::cout<<(A1 - A).zx() - (B1 - B).zx()<<std::endl;
  std::cout<<(A1 - A).zy() - (B1 - B).zy()<<std::endl;

  srand(time(NULL));
  double val = ( rand() % 1000 ) * 1.0e-3 - 0.5;
  A1.scale(val); B1.scale(val);
  std::cout<<A1.xx() - B1.xx()<<std::endl;
  std::cout<<A1.xy() - B1.xy()<<std::endl;
  std::cout<<A1.xz() - B1.xz()<<std::endl;
  std::cout<<A1.yy() - B1.yy()<<std::endl;
  std::cout<<A1.yz() - B1.yz()<<std::endl;
  std::cout<<A1.zz() - B1.zz()<<std::endl;
  std::cout<<A1.yx() - B1.yx()<<std::endl;
  std::cout<<A1.zx() - B1.zx()<<std::endl;
  std::cout<<A1.zy() - B1.zy()<<std::endl;

  for (int ii=0; ii<3; ++ii) A1 += A;
  for (int ii=0; ii<3; ++ii) B1 += B;
  std::cout<<A1.xx() - B1.xx()<<std::endl;
  std::cout<<A1.xy() - B1.xy()<<std::endl;
  std::cout<<A1.xz() - B1.xz()<<std::endl;
  std::cout<<A1.yy() - B1.yy()<<std::endl;
  std::cout<<A1.yz() - B1.yz()<<std::endl;
  std::cout<<A1.zz() - B1.zz()<<std::endl;
  std::cout<<A1.yx() - B1.yx()<<std::endl;
  std::cout<<A1.zx() - B1.zx()<<std::endl;
  std::cout<<A1.zy() - B1.zy()<<std::endl;

  for (int ii=0; ii<3; ++ii) A1 -= A;
  for (int ii=0; ii<3; ++ii) B1 -= B;
  std::cout<<A1.xx() - B1.xx()<<std::endl;
  std::cout<<A1.xy() - B1.xy()<<std::endl;
  std::cout<<A1.xz() - B1.xz()<<std::endl;
  std::cout<<A1.yy() - B1.yy()<<std::endl;
  std::cout<<A1.yz() - B1.yz()<<std::endl;
  std::cout<<A1.zz() - B1.zz()<<std::endl;
  std::cout<<A1.yx() - B1.yx()<<std::endl;
  std::cout<<A1.zx() - B1.zx()<<std::endl;
  std::cout<<A1.zy() - B1.zy()<<std::endl;

  return EXIT_SUCCESS;
}

// EOF
