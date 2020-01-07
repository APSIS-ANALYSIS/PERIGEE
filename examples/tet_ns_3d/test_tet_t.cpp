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
  std::vector<int> ienout, ptidx, elemidx;
  for(int ii=0; ii<20; ++ii) 
  {
    ptout.push_back(ptcoor[3*ien[ii]+0]);
    ptout.push_back(ptcoor[3*ien[ii]+1]);
    ptout.push_back(ptcoor[3*ien[ii]+2]);
    ienout.push_back(ii);
    ptidx.push_back(ii*10);
  }

  elemidx.push_back(-1);
  elemidx.push_back(231);

  std::string out_name("./test-output");

  TET_T::write_tet_grid_node_elem_index( out_name, 20, 2, ptout, ienout,
     ptidx, elemidx );

  PetscFinalize();
  return 0;
}
