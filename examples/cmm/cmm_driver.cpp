// ==================================================================
// cmm_driver
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"

int main( int argc, char *argv[] )
{
  int nqp_tet = 5, nqp_tri = 4;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
