#include "Tet_Tools.hpp"

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  std::string geo_file("./test-surface.vtu");

  int nFunc, nElem;
  std::vector<double> ptcoor;
  std::vector<int> ien, nidx, eidx;
  TET_T::read_vtu_grid(geo_file, nFunc, nElem, ptcoor, ien, nidx, eidx );

  cout<<nFunc<<'\t'<<nElem<<'\n';

  VEC_T::print(ptcoor);

  VEC_T::print(ien);
  VEC_T::print(nidx);
  VEC_T::print(eidx);

  std::vector<double> ptout;
  std::vector<int> ienout, ptidx, elemidx, elemidx2;

  ptout.push_back( 0.0 ); ptout.push_back( 0.0 ); ptout.push_back( 0.1 );
  ptout.push_back( 1.0 ); ptout.push_back( 0.0 ); ptout.push_back( 0.0 );
  ptout.push_back( 0.0 ); ptout.push_back( 1.0 ); ptout.push_back( 0.2 );
  ptout.push_back( 0.5 ); ptout.push_back( 0.0 ); ptout.push_back( 0.0 );
  ptout.push_back( 0.5 ); ptout.push_back( 0.3 ); ptout.push_back( 0.3 );
  ptout.push_back( 0.0 ); ptout.push_back( 0.7 ); ptout.push_back( 0.0 );
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
  elemidx2.push_back(31);

  std::string out_name("./test-surface");

  std::vector<int> ptag;
  ptag.push_back(2);
  ptag.push_back(4);

  TET_T::write_quadratic_triangle_grid( out_name, 6, 1, ptout, ienout,
      ptidx, elemidx );

  PetscFinalize();
  return 0;
}
