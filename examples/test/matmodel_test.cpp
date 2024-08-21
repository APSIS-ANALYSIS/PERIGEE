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

  const auto S1_iso = mat1.get_PK_2nd( F );

  Tensor2_3D P1_iso;

  const auto CC1_iso = mat1.get_PK_Stiffness( F, P1_iso );

  MaterialModel_StVenant_Kirchhoff_M94_Mixed mat2( E, nu );

  Tensor2_3D P2_iso, S2_iso;
  Tensor4_3D CC2_iso;

  mat2.get_PK_Stiffness( F, P2_iso, S2_iso, CC2_iso );

  const auto Sym_S2iso = SymmTensor2_3D( 
    S2_iso(0, 0), S2_iso(1, 1), S2_iso(2, 2),
    S2_iso(1, 2), S2_iso(0, 2), S2_iso(0, 1) );

  std::cout << "StVenant: \n";

  std::cout << "error for S_iso: \n"; 
  ( S1_iso - Sym_S2iso ).print();

  std::cout << "error for P_iso: \n"; 
  ( P1_iso - P2_iso ).print();

  std::cout << "error for CC_iso: \n"; 
  ( CC1_iso.full() - CC2_iso ).print_in_mat();



  // GOH06
  const double f1_the=49.98, f1_phi=49.98, f2_the=-49.98, f2_phi=-49.98;
  const double fk1=996.6, fk2=0.1, fkd=1.0 /3.0; 

  MaterialModel_ich_GOH14 mat3( mu, f1_the, f1_phi, f2_the, f2_phi,
    fk1, fk2, fkd );

  const auto S3_iso = mat3.get_PK_2nd( F );

  Tensor2_3D P3_iso;

  const auto CC3_iso = mat3.get_PK_Stiffness( F, P3_iso );

  MaterialModel_GOH14_ST91_Mixed mat4( 1.0, E, nu, f1_the, f1_phi, f2_the, f2_phi,
    fk1, fk2, fkd );

  Tensor2_3D P4_iso, S4_iso;
  Tensor4_3D CC4_iso;

  mat4.get_PK_Stiffness( F, P4_iso, S4_iso, CC4_iso );

  const auto Sym_S4iso = SymmTensor2_3D( 
    S4_iso(0, 0), S4_iso(1, 1), S4_iso(2, 2),
    S4_iso(1, 2), S4_iso(0, 2), S4_iso(0, 1) );

  std::cout << "\nGOH06: \n";

  std::cout << "error for S_iso: \n"; 
  ( S3_iso - Sym_S4iso ).print();

  std::cout << "error for P_iso: \n"; 
  ( P3_iso - P4_iso ).print();

  std::cout << "error for CC_iso: \n"; 
  ( CC3_iso.full() - CC4_iso ).print_in_mat();








  return EXIT_SUCCESS;
}

// EOF

