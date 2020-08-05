// Test the new TET_T functions

#include "ElemBC_3D_tet_wall.hpp"

using namespace std;

int main( int argc, char *argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  ElemBC * wall_bc = new ElemBC_3D_tet_wall( "wall_cyl.vtp",
     "centerlines.vtp", 0.2 );
  
  delete wall_bc;
  PetscFinalize();
  return 0;
}

//EOF
