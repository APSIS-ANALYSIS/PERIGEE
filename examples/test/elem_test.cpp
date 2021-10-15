#include "FEAElement_Triangle6.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "QuadPts_debug.hpp"

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  const double r = 0.237;
  const double s = 0.15;
  IQuadPts * quadv = new QuadPts_debug( 4, 1, {r, s, 0.0, 1.0-r-s-0.0}, {1.0} );
  IQuadPts * quads = new QuadPts_debug( 3, 1, {r, s, 1.0-r-s}, {1.0});

  FEAElement * elementv = new FEAElement_Tet10_v2( 1 );
  FEAElement * elements = new FEAElement_Triangle6( 1 );

  std::vector<double> ept_x {0.0, 1.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.0};
  std::vector<double> ept_y {0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5};
  std::vector<double> ept_z {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5};

  elementv -> buildBasis( quadv, &ept_x[0], &ept_y[0], &ept_z[0] );

  std::vector<double> spt_x { ept_x[0], ept_x[1], ept_x[2], ept_x[4], ept_x[5], ept_x[6] };
  std::vector<double> spt_y { ept_y[0], ept_y[1], ept_y[2], ept_y[4], ept_y[5], ept_y[6] };

  elements -> buildBasis( quads, &spt_x[0], &spt_y[0] );

  VEC_T::print( elementv->get_d2R_dxy(0) );
  VEC_T::print( elements->get_d2R_dxy(0) );

  delete quadv; delete quads; delete elementv; delete elements;
  PetscFinalize();
  return EXIT_SUCCESS;
}
