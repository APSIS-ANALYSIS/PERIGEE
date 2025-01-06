#include "Vector_3.hpp"
#include "SymmTensor2_3D.hpp"

int main( int argc, char * argv[] )
{
  //Test for unary minus operator of Vector_3 
  std::cout << "Test for unary minus operator of Vector_3: \n";

  Vector_3 A; A.gen_rand();

  std::cout << "current A: \n";
  A.print();

  std::cout << "unary operator -A: \n";
  (-A).print();

  A.gen_rand();
  std::cout << "current A: \n";
  A.print();  

  std::cout << "unary operator --A: \n";
  (-(-A)).print();

  
  A.gen_rand(); Vector_3 B ( A );

  std::cout << "current A: \n";
  A.print();  

  std::cout << "current B: \n";
  B.print();  

  std::cout << "unary operator -A+B: \n";
  (-A + B).print();

  A.gen_rand(); B = A;

  std::cout << "current A: \n";
  A.print();  

  std::cout << "current B: \n";
  B.print();    

  std::cout << "unary operator --A + B: \n";
  (-(-A) + B).print();  

  return EXIT_SUCCESS;
}

// EOF
