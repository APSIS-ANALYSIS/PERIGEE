// Test the new TET_T functions

#include "ElemBC_3D_tet_wall.hpp"

using namespace std;

int main( int argc, char *argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  ElemBC * wall_bc = new ElemBC_3D_tet( "wall_cyl.vtp" );
  
  std::vector<std::string> walllist; walllist.push_back("wall_cyl.vtp");

  ElemBC * wall_bc_2 = new ElemBC_3D_tet( walllist );

  wall_bc -> print_info();

  wall_bc_2 -> print_info();

  delete wall_bc; delete wall_bc_2;
  PetscFinalize();
  return 0;
}

//EOF
