#include "Vec_Tools.hpp"
#include "Matrix_3x3.hpp"

int main( int argc, char * argv[] )
{
  std::vector<std::string> a;
  a.push_back("hello a\n");
  a.push_back("cmame \t");
  a.push_back("july 31th \n");
  a.push_back("aug 31th \n");

  for(auto &out : a) out = "hello\n";

  for(auto out : a) std::cout<<out;

  return EXIT_SUCCESS;
}

// EOF
