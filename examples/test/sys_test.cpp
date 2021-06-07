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

  bool flg = false;

  int aa = ( flg ? 1 : 0);
  std::cout<<aa<<'\n';

  return EXIT_SUCCESS;
}

// EOF
