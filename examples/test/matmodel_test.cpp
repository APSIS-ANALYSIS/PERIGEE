#include "FEAElementFactory.hpp"
#include "QuadPtsFactory.hpp"

int main( int argc, char * argv[] )
{
  const FEType type = FEType::Hex27; 
  const int nqp_vol = 8;
  const int nqp_sur = 4;
  const int nqp_vol_1D = 2;
  const int nqp_sur_1D = 2;

  FEAElement * elementv = nullptr;
  FEAElement * elements = nullptr;
  IQuadPts * quadv = nullptr;
  IQuadPts * quads = nullptr;

  if( type  == FEType::Tet4 )
  {
    if( nqp_vol > 5 ) SYS_T::commPrint("Warning: the tet element is linear and you are using more than 5 quadrature points.\n");
    if( nqp_sur > 4 ) SYS_T::commPrint("Warning: the tri element is linear and you are using more than 4 quadrature points.\n");

    elementv = new FEAElement_Tet4( nqp_vol ); // elem type Tet4
    elements = new FEAElement_Triangle3_3D_der0( nqp_sur );
    quadv = new QuadPts_Gauss_Tet( nqp_vol );
    quads = new QuadPts_Gauss_Triangle( nqp_sur );
  }
  else if( type == FEType::Tet10 )
  {
    elementv = new FEAElement_Tet10( nqp_vol ); // elem type Tet10
    elements = new FEAElement_Triangle6_3D_der0( nqp_sur );
    quadv = new QuadPts_Gauss_Tet( nqp_vol );
    quads = new QuadPts_Gauss_Triangle( nqp_sur );
  }
  else if( type == FEType::Hex8 )
  {
    elementv = new FEAElement_Hex8( nqp_vol_1D * nqp_vol_1D * nqp_vol_1D ); // elem type Hex8
    elements = new FEAElement_Quad4_3D_der0( nqp_sur_1D * nqp_sur_1D );
    quadv = new QuadPts_Gauss_Hex( nqp_vol_1D );
    quads = new QuadPts_Gauss_Quad( nqp_sur_1D );
  }
  else if( type == FEType::Hex27 )
  {
    elementv = new FEAElement_Hex27( nqp_vol_1D * nqp_vol_1D * nqp_vol_1D ); // elem type Hex27
    elements = new FEAElement_Quad9_3D_der0( nqp_sur_1D * nqp_sur_1D );
    quadv = new QuadPts_Gauss_Hex( nqp_vol_1D );
    quads = new QuadPts_Gauss_Quad( nqp_sur_1D );
  }
  else SYS_T::print_fatal("Error: Element type not supported.\n");

  // print the information of element and quadrature rule
  //elementv->print_info();
  quadv->print_info();

  // new appraoch
  auto qv = QuadPtsFactory::createVolQuadrature(type, nqp_vol);
  qv -> print_info();

  quads->print_info();
  auto qs = QuadPtsFactory::createSurQuadrature(type, nqp_sur);
  qs->print_info();

  auto qq = QuadPtsFactory::createVisQuadrature(type);
  qq->print_info();

  delete elementv; delete elements; delete quadv; delete quads;
  return EXIT_SUCCESS;
}

// EOF

