// This is a simple driver for testing the output of the membrane element
// routines.
#include "QuadPts_debug.hpp"
#include "FEAElement_Triangle3_membrane.hpp"
#include "FEAElement_Triangle6_membrane.hpp"

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  std::vector<double> in_qp{ 0.5, 0.2, 0.3 };
  std::vector<double> in_qw{ 1.0, 1.0, 1.0 };
  
  IQuadPts * quad = new QuadPts_debug(3, 1, in_qp, in_qw);

  quad -> print_info();

  FEAElement * elem = new FEAElement_Triangle3_membrane( 3 );

  delete quad;
  delete elem;

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
