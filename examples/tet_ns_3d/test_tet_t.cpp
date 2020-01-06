#include "Tet_Tools.hpp"

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  std::string geo_file("./mesh-complete.mesh.vtu");

  int nFunc, nElem;
  std::vector<double> ptcoor;
  std::vector<int> ien;

  TET_T::read_vtu_grid(geo_file, nFunc, nElem, ptcoor, ien );
  

  std::vector<double> ptout;
  std::vector<int> ienout;
  for(int ii=0; ii<10; ++ii) 
  {
    ptout.push_back(ptcoor[3*ien[ii]+0]);
    ptout.push_back(ptcoor[3*ien[ii]+1]);
    ptout.push_back(ptcoor[3*ien[ii]+2]);
    ienout.push_back(ii);
  }

  std::string out_name("./test-output.vtu");

  TET_T::write_tet_grid( out_name, 10, 1, ptout, ienout );

  PetscFinalize();
  return 0;
}
