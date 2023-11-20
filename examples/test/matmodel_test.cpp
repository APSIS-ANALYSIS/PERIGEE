#include "MaterialModel_Guccione_Incompressible_Mixed.hpp"

int main( int argc, char * argv[] )
{
#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  IMaterialModel * model = new MaterialModel_Guccione_Incompressible_Mixed(
      2.1, 3.3, 22.0, -11.2, 15.2, 0.5, 0.3, 0.2, -0.2, 0.3, 0.6 );
  
  model -> print_info();

  Tensor2_3D F, P, S;
  Tensor4_3D CC;

  F.xx() = 1.0; F.xy() = 0.3; F.xz() = 0.0;
  F.yx() = 0.0; F.yy() = 1.0; F.yz() = -0.05;
  F.zx() = 0.0; F.zy() = 0.1; F.zz() = 1.0;

  model -> get_PK_Stiffness(F, P, S, CC);

  P.print_in_row();
  S.print_in_row();
  CC.print();

  model -> get_PK( F, P, S );

  P.print_in_row();
  S.print_in_row();
 
  std::cout<<model -> get_strain_energy( F )<<std::endl; 
  
  delete model;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF

