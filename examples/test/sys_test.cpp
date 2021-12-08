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

  Matrix_3x3 A;

  A.gen_rand();

  Matrix_3x3 B = cofactor(A);

  Matrix_3x3 C( A );

  C.inverse();
  C.transpose();
  C.scale(A.det());

  C -= B;

  C.print_in_row();

  return EXIT_SUCCESS;
}

// EOF
