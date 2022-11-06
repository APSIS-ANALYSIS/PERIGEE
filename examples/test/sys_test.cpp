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

  std::cout<<"eigen decomp test: \n";

  SymmMatrix_3x3 A1;
  A1.gen_rand();

  double eta1 = 0.0; double eta2 = 0.0; double eta3 = 0.0;
  Vector_3 v1; Vector_3 v2; Vector_3 v3;
  int out = A1.eigen_decomp(eta1, eta2, eta3, v1, v2, v3);
  std::cout << "case:" << out << std::endl;
  
  Matrix_3x3 Basis1; Matrix_3x3 Basis2; Matrix_3x3 Basis3;
  Basis1.gen_outprod(v1);
  Basis2.gen_outprod(v2);
  Basis3.gen_outprod(v3);

  Matrix_3x3 A1_ = eta1 * Basis1 + eta2 * Basis2 + eta3 * Basis3;
  Matrix_3x3 A1full = A1.convert_to_full();
  (A1full - A1_).print();  
  
  double trA1 = A1.tr();  double frac_13= 1.0/3.0;
  SymmMatrix_3x3 A1p(A1);
  A1p(0) -= frac_13 * trA1;
  A1p(1) -= frac_13 * trA1;
  A1p(2) -= frac_13 * trA1;

  Vector_3 v1p; Vector_3 s1; Vector_3 s2;
  double eta1p = eta1 - frac_13 * trA1;
  A1p.find_eigen_vector(eta1p, v1p, s1, s2);  
  (v1p - v1).print(); 

  // Compare with Matrix_3x3
  double eta1b = 0.0; double eta2b = 0.0; double eta3b = 0.0;
  Vector_3 v1b; Vector_3 v2b; Vector_3 v3b;
  Matrix_3x3 B1( A1.xx(), A1.xy(), A1.xz(), A1.yx(), A1.yy(), A1.yz(), A1.zx(), A1.zy(), A1.zz() );
  (A1full - B1).print(); // test of convert_to_full()
  
  B1.eigen_decomp(eta1b, eta2b, eta3b, v1b, v2b, v3b);
  (v3b - v3).print();
  (v2b - v2).print();
  (v1b - v1).print();
  std::cout << (eta1b - eta1) << std::endl;
  std::cout << (eta2b - eta2) << std::endl;
  std::cout << (eta3b - eta3) << std::endl;

  // test-case: three eigenvalues are the same
  v1.gen_e1(); v2.gen_e2(); v3.gen_e3();
  Basis1.gen_outprod(v1);
  Basis2.gen_outprod(v2);
  Basis3.gen_outprod(v3);

  eta1 = val; eta2 = val; eta3 = val;
  Matrix_3x3 B2 = eta1 * Basis1 + eta2 * Basis2 + eta3 * Basis3;
  
  double eta11 = 0.0; double eta22 = 0.0; double eta33 = 0.0;
  Vector_3 v11; Vector_3 v22; Vector_3 v33;
  SymmMatrix_3x3 A2 = gen_symm_part(B2);
  out = A2.eigen_decomp(eta11, eta22, eta33, v11, v22, v33);

  std::cout << "case:" << out << std::endl;
  (v33 - v3).print();
  (v22 - v2).print();
  (v11 - v1).print();  
  std::cout << (eta1 - eta11) << std::endl;
  std::cout << (eta2 - eta22) << std::endl;
  std::cout << (eta3 - eta33) << std::endl;

  // test-case2: two eigenvalues are the same
  eta1 = 2 * val; eta2 = val; eta3 = val;
  v1.gen_rand(); 
  v1.normalize();
  v2(0) = val;
  v2(1) = 2 * val;
  v2(2) = - ( v1(0) * v2(0) + v1(1) * v2(1) ) / ( v1(2) ) ; 
  v2.normalize();
  v3 = cross_product(v1, v2);
  v3.normalize();

  Basis1.gen_outprod(v1);
  Basis2.gen_outprod(v2);
  Basis3.gen_outprod(v3);
  Matrix_3x3 B3 = eta1 * Basis1 + eta2 * Basis2 + eta3 * Basis3;
  SymmMatrix_3x3 A3 = gen_symm_part(B3);
 
  double eta111 = 0.0; double eta222 = 0.0; double eta333 = 0.0;
  Vector_3 v111; Vector_3 v222; Vector_3 v333;
  out = A3.eigen_decomp(eta111, eta222, eta333, v111, v222, v333);

  Basis1.gen_outprod(v111);
  Basis2.gen_outprod(v222);
  Basis3.gen_outprod(v333);
  Matrix_3x3 B3_ = eta111 * Basis1 + eta222 * Basis2 + eta333 * Basis3;
  std::cout << "case:" << out << std::endl;   
  std::cout << (eta1 - eta111) << std::endl;
  std::cout << (eta2 - eta222) << std::endl;
  std::cout << (eta3 - eta333) << std::endl;
  SymmMatrix_3x3 A3_ = gen_symm_part(B3_);
  (A3 - A3_).print();

  return EXIT_SUCCESS;
}

// EOF
