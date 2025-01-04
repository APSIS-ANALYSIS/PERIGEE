#include "FEAElementFactory.hpp"
#include "QuadPtsFactory.hpp"

int main( int argc, char * argv[] )
{
  const FEType type = FEType::Hex27; 
  const int nqp_vol = 216;
  const int nqp_sur = 36;
  const int nqp_vol_1D = 6;
  const int nqp_sur_1D = 6;

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
    SYS_T::print_fatal_if( nqp_vol < 29, "Error: not enough quadrature points for tets.\n" );
    SYS_T::print_fatal_if( nqp_sur < 13, "Error: not enough quadrature points for triangles.\n" );

    elementv = new FEAElement_Tet10( nqp_vol ); // elem type Tet10
    elements = new FEAElement_Triangle6_3D_der0( nqp_sur );
    quadv = new QuadPts_Gauss_Tet( nqp_vol );
    quads = new QuadPts_Gauss_Triangle( nqp_sur );
  }
  else if( type == FEType::Hex8 )
  {
    SYS_T::print_fatal_if( nqp_vol_1D < 2, "Error: not enough quadrature points for hex.\n" );
    SYS_T::print_fatal_if( nqp_sur_1D < 1, "Error: not enough quadrature points for quad.\n" );

    elementv = new FEAElement_Hex8( nqp_vol_1D * nqp_vol_1D * nqp_vol_1D ); // elem type Hex8
    elements = new FEAElement_Quad4_3D_der0( nqp_sur_1D * nqp_sur_1D );
    quadv = new QuadPts_Gauss_Hex( nqp_vol_1D );
    quads = new QuadPts_Gauss_Quad( nqp_sur_1D );
  }
  else if( type == FEType::Hex27 )
  {
    SYS_T::print_fatal_if( nqp_vol_1D < 4, "Error: not enough quadrature points for hex.\n" );
    SYS_T::print_fatal_if( nqp_sur_1D < 3, "Error: not enough quadrature points for quad.\n" );

    elementv = new FEAElement_Hex27( nqp_vol_1D * nqp_vol_1D * nqp_vol_1D ); // elem type Hex27
    elements = new FEAElement_Quad9_3D_der0( nqp_sur_1D * nqp_sur_1D );
    quadv = new QuadPts_Gauss_Hex( nqp_vol_1D );
    quads = new QuadPts_Gauss_Quad( nqp_sur_1D );
  }
  else SYS_T::print_fatal("Error: Element type not supported.\n");

  // print the information of element and quadrature rule
  elementv->print_info();
  //quadv->print_info();
  //quads->print_info();

  // new appraoch
  const auto ev = ElementFactory::createVolumeElement( type, nqp_vol );

  ev -> print_info();

  std::cout<<elementv->get_nLocBas()<<'\t'<<FE_T::to_nLocBas(type)<<'\n';
  std::cout<<elementv->get_numQuapts()<<'\t'<<ev->get_numQuapts()<<'\n';

  elements->print_info();
  const auto es = ElementFactory::createSurfaceElement( type, nqp_sur );
  es -> print_info();
  std::cout<<elements->get_nLocBas()<<'\t'<<es->get_nLocBas()<<'\n';
  std::cout<<elements->get_numQuapts()<<'\t'<<es->get_numQuapts()<<'\n';

  delete elementv; delete elements; delete quadv; delete quads;
  return EXIT_SUCCESS;
}

// EOF

