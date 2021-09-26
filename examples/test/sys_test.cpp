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
  Matrix_3x3 A(
      1.0, 5.2, -0.33,
      5.2, 2.0, 3.0,
     -0.33, 3.0, -1.0);

  va.gen_rand();

  vb = A * va;

  //A.print();

  double eta1, eta2, eta3;
  Vector_3 v1, v2, v3;

  std::cout<<A.eigen_decomp( eta1, eta2, eta3, v1, v2, v3 )<<std::endl;

  std::cout<<std::setprecision(16)<<eta1<<'\t'<<eta2<<'\t'<<eta3<<'\n';

  std::cout<<std::setprecision(16)<<v1(0)<<'\t'<<v1(1)<<'\t'<<v1(2)<<'\n';

  v2.print();

  v3.print();

  return EXIT_SUCCESS;
}

// EOF
