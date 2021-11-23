#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "HDF5_Reader.hpp"
#include "PostVectSolution.hpp"

int main( int argc, char * argv[] )
{
  std::vector<std::string> a;
  a.push_back("hello a\n");
  a.push_back("cmame \t");
  a.push_back("july 31th \n");
  a.push_back("aug 31th \n");

  for(auto &out : a) out = "hello\n";

  //for(auto out : a) std::cout<<out;


  Vector_3 va, vb, vc;
  Matrix_3x3 A(
      1.0, 5.2, -0.33,
      5.2, 2.0, 3.0,
     -0.33, 3.0, -1.0);

  A.gen_rand();

  sleep(1);

  Matrix_3x3 B;
  B.gen_rand();
  
  A.print();
  B.print();

  Matrix_3x3 C,D;
  C.MatMult(A,B);

  D = A * B;

  D -= C;

  D.print();

  return EXIT_SUCCESS;
}

// EOF
