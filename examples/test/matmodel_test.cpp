#include "Sys_Tools.hpp"
#include "MaterialModel_ich_StVenant_Kirchhoff.hpp"
#include "MaterialModel_StVenant_Kirchhoff_M94_Mixed.hpp"
#include "MaterialModel_ich_GOH14.hpp"
#include "MaterialModel_GOH14_ST91_Mixed.hpp"

int main( int argc, char * argv[] )
{
  Tensor2_3D F; F.gen_rand();

  if(F.det() < 0.0) F*=-1.0;

  const double mu = 1.0;

  const double nu = 0.4;

  const double E = mu * ( 2.0 + 2.0*nu );

  MaterialModel_ich_StVenant_Kirchhoff mat1( mu );

  const auto S1_ich = mat1.get_PK_2nd( F );

  Tensor2_3D P1_ich;

  const auto CC1_ich = mat1.get_PK_Stiffness( F, P1_ich );

  MaterialModel_StVenant_Kirchhoff_M94_Mixed mat2( E, nu );

  Tensor2_3D P2_ich, S2_ich;
  Tensor4_3D CC2_ich;

  mat2.get_PK_Stiffness( F, P2_ich, S2_ich, CC2_ich );

  const auto Sym_S2iso = SymmTensor2_3D( 
    S2_ich(0, 0), S2_ich(1, 1), S2_ich(2, 2),
    S2_ich(1, 2), S2_ich(0, 2), S2_ich(0, 1) );

  std::cout << "StVenant: \n";

  std::cout << "error for S_ich: \n"; 
  ( S1_ich - Sym_S2iso ).print();

  std::cout << "error for P_ich: \n"; 
  ( P1_ich - P2_ich ).print();

  std::cout << "error for CC_ich: \n"; 
  ( CC1_ich.full() - CC2_ich ).print_in_mat();



  // GOH06
  const double f1_the=49.98, f1_phi=49.98, f2_the=-49.98, f2_phi=-49.98;
  const double fk1=996.6, fk2=0.1, fkd=1.0 /3.0; 

  MaterialModel_ich_GOH14 mat3( mu, f1_the, f1_phi, f2_the, f2_phi,
    fk1, fk2, fkd );

  const auto S3_ich = mat3.get_PK_2nd( F );

  Tensor2_3D P3_ich;

  const auto CC3_ich = mat3.get_PK_Stiffness( F, P3_ich );

  MaterialModel_GOH14_ST91_Mixed mat4( 1.0, E, nu, f1_the, f1_phi, f2_the, f2_phi,
    fk1, fk2, fkd );

  Tensor2_3D P4_ich, S4_ich;
  Tensor4_3D CC4_ich;

  mat4.get_PK_Stiffness( F, P4_ich, S4_ich, CC4_ich );

  const auto Sym_S4iso = SymmTensor2_3D( 
    S4_ich(0, 0), S4_ich(1, 1), S4_ich(2, 2),
    S4_ich(1, 2), S4_ich(0, 2), S4_ich(0, 1) );

  std::cout << "\nGOH14: \n";

  std::cout << "error for S_ich: \n"; 
  ( S3_ich - Sym_S4iso ).print();

  std::cout << "error for P_ich: \n"; 
  ( P3_ich - P4_ich ).print();

  std::cout << "error for CC_ich: \n"; 
  ( CC3_ich.full() - CC4_ich ).print_in_mat();








  return EXIT_SUCCESS;
}

// EOF

