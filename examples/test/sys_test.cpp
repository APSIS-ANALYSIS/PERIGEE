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

  Matrix_3x3 invA( A ); invA.inverse();

  Matrix_3x3 B = inverse(A);

  B -= invA;

  B.print();
  invA.print();

  double eta1, eta2, eta3;
  Vector_3 v1, v2, v3;

  std::cout<<A.eigen_decomp( eta1, eta2, eta3, v1, v2, v3 )<<std::endl;

  std::cout<<std::setprecision(16)<<eta1<<'\t'<<eta2<<'\t'<<eta3<<'\n';

  std::cout<<std::setprecision(16)<<v1(0)<<'\t'<<v1(1)<<'\t'<<v1(2)<<'\n';

  v2.print();

  v3.print();

  int nsize = 4;

  std::vector<int> b (nsize, -100), c (5, 99), d (2, 1);

  std::cout<<b.size()<<'\t'<<b.capacity()<<std::endl;
  
  VEC_T::print(b);

  b = {1, 3, 5, -1, -2, -3, -5};
  
  VEC_T::print(b);

  std::cout<<b.size()<<'\t'<<b.capacity()<<'\t'<<VEC_T::sum(b)<<std::endl;

  return EXIT_SUCCESS;
}

// EOF
