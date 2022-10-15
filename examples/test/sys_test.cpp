#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "SymmMatrix_3x3.hpp"
#include "HDF5_Reader.hpp"
#include "PostVectSolution.hpp"


int main(int argc, char *argv[])
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

  sleep(1);
  srand(time(NULL));
  double val = ( rand() % 1000 ) * 1.0e-3 - 0.5;

  std::cout<<"rotation test: \n";
  const double pi = 4.0 * atan(1.0);
  double thetay = pi/val;
  double thetax = 2*pi/val;
  double thetaz = 3*pi/val;
  Matrix_3x3 Ry(cos(thetay),   0,   sin(thetay),
                     0,        1,       0,
               -sin(thetay),   0,   cos(thetay));

  Matrix_3x3 Rx(1,         0,           0,
                0,   cos(thetax),  -sin(thetax),
                0,   sin(thetax),   cos(thetax));

  Matrix_3x3 Rz(cos(thetaz), -sin(thetaz),   0,
                sin(thetaz),  cos(thetaz),   0,
                     0,          0,          1);
  Matrix_3x3 Q;
  Q.MatMult(Rx, Ry);
  Q.MatMult(Q, Rz);

  return EXIT_SUCCESS;
}

// EOF
