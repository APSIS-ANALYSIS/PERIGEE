#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "HDF5_Reader.hpp"

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
  Matrix_3x3 A;
  A.gen_rand();

  va.gen_rand();

  vb = A.VecMultT(va);

  A.VecMultT(va(0), va(1), va(2), vc(0), vc(1), vc(2) );

  A.print();

  va.print();

  std::cout<<dist(vb, vc)<<std::endl;
  return EXIT_SUCCESS;
}

// EOF
