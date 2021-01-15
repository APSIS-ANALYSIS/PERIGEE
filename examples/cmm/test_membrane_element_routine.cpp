// This is a simple driver for testing the output of the membrane element
// routines.
#include "QuadPts_debug.hpp"
#include "FEAElement_Triangle3_membrane.hpp"
#include "FEAElement_Triangle6_membrane.hpp"

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  const int dim   = 3;
  const int numpt = 4;

  std::vector<double> in_qp{ 1.0/3.0, 1.0/3.0, 1.0/3.0, 0.6, 0.2, 0.2, 0.2, 0.6, 0.2, 0.2, 0.2, 0.6 };
  std::vector<double> in_qw{ -0.5625, 0.520833333333333, 0.520833333333333, 0.520833333333333 };
  
  IQuadPts * quad = new QuadPts_debug(dim, numpt, in_qp, in_qw);

  quad -> print_info();

  FEAElement * elem = new FEAElement_Triangle3_membrane( 3 );
  double ctrl_x[3] = {  0.3971,  0.4969, 0.4516 };
  double ctrl_y[3] = { -1.4233, -1.2942, 1.3001 };
  double ctrl_z[3] = {  9.7337,  9.6558, 9.8612 };

  elem -> buildBasis( quad, ctrl_x, ctrl_y, ctrl_z );

  delete quad;
  delete elem;

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
