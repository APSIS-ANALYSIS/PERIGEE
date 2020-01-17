#include "Tet_Tools.hpp"

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  std::string geo_file("./walls_combined.vtp");

  int nFunc, nElem;
  std::vector<double> ptcoor;
  std::vector<int> ien;
  TET_T::read_vtp_grid(geo_file, nFunc, nElem, ptcoor, ien );

  std::vector<double> ptout;
  std::vector<int> ienout, ptidx, elemidx;

  ptout.push_back( 0.0 ); ptout.push_back( 0.0 ); ptout.push_back( 0.0 );
  ptout.push_back( 1.0 ); ptout.push_back( 0.0 ); ptout.push_back( 0.0 );
  ptout.push_back( 0.0 ); ptout.push_back( 1.0 ); ptout.push_back( 0.0 );
  ptout.push_back( 0.5 ); ptout.push_back( 0.0 ); ptout.push_back( 0.0 );
  ptout.push_back( 0.5 ); ptout.push_back( 0.5 ); ptout.push_back( 0.0 );
  ptout.push_back( 0.0 ); ptout.push_back( 0.5 ); ptout.push_back( 0.0 );
  ienout.push_back(0);
  ienout.push_back(1);
  ienout.push_back(2);
  ienout.push_back(3);
  ienout.push_back(4);
  ienout.push_back(5);
  ptidx.push_back(0);
  ptidx.push_back(1);
  ptidx.push_back(2);
  ptidx.push_back(3);
  ptidx.push_back(4);
  ptidx.push_back(0);

  //elemidx.push_back(-1);
  elemidx.push_back(231);

  std::string out_name("./test-surface");

  std::vector<int> ptag;
  ptag.push_back(2);
  ptag.push_back(4);

  TET_T::write_triangle_grid( out_name, 6, 1, ptout, ienout,
      ptidx, elemidx );

  PetscFinalize();
  return 0;
}
