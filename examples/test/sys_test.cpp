#include "Vec_Tools.hpp"
#include "Matrix_3x3.hpp"

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  
  Matrix_3x3 A, Q;

  A = { 1,0,1, -1, -2, 0, 0, 1, -1 };

  Q = { -0.120000260381534, -0.809712281592778, 0.574426634607224, 0.901752646908814, 0.153122822484370, 0.404222172854692, -0.415261485453819, 0.566497504206538, 0.711785414592383 };

  A.MatRot(Q);

  A.print();  
  
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
