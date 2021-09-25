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

  for(auto out : a) std::cout<<out;

  
  SYS_T::execute("rm -rf part_cmmbc_0");
  //SYS_T::execute("mkdir part_cmmbc_0");

  std::string st( "part_cmmbc_0" );

  st.append("/");

  std::cout<<st<<std::endl;

  Vector_3 va(1.0, 2.0, 0.0);
  Vector_3 vb(0.0, 1.0, 0.0);

  std::cout<<dist(va,vb)<<std::endl;

  return EXIT_SUCCESS;
}

// EOF
