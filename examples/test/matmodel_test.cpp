#include "MaterialModel_Guccione_Incompressible_Mixed.hpp"

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  IMaterialModel * model = new MaterialModel_Guccione_Incompressible_Mixed(
      2.1, 3.3, 22.0, -11.2, 15.2, 0.5, 0.3, 0.2, -0.2, 0.3, 0.6 );
  
  model -> print_info();

  model -> write_hdf5("hello.h5");

  delete model;
  
  IMaterialModel * model2 = new MaterialModel_Guccione_Incompressible_Mixed("hello.h5");

  model2 -> print_info();

  delete model2;

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF

