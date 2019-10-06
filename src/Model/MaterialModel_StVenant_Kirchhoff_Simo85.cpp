#include "MaterialModel_StVenant_Kirchhoff_Simo85.hpp"

MaterialModel_StVenant_Kirchhoff_Simo85::MaterialModel_StVenant_Kirchhoff_Simo85( 
    const double &in_E, const double &in_nu )
: E( in_E ), nu( in_nu ), lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{
}


MaterialModel_StVenant_Kirchhoff_Simo85::~MaterialModel_StVenant_Kirchhoff_Simo85()
{}


void MaterialModel_StVenant_Kirchhoff_Simo85::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  MaterialModel_StVenant_Kirchhoff_Simo85: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Young's Modulus E  = %e \n", E);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Possion's ratio nu = %e \n", nu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Lame coeff lambda  = %e \n", lambda);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Shear modulus mu   = %e \n", mu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Bulk modulus kappa = %e \n", kappa);
}


void MaterialModel_StVenant_Kirchhoff_Simo85::get_PK( 
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S )
{
  C.MatMultTransposeLeft(F);
  detF = F.det();
  S.copy(C);
  S.inverse();
  S.scale(kappa * std::log(detF));
  S.AXPY(mu, C);
  S.AXPY(-1.0 * mu, I);
  P.MatMult(F,S);
}


void MaterialModel_StVenant_Kirchhoff_Simo85::get_PK_Stiffness( 
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S, Tensor4_3D &CC)
{
  CC.gen_zero();
  const double val = -2.0 * kappa * std::log(detF);
  
  C.MatMultTransposeLeft(F);
  detF = F.det();
  S.copy(C);
  S.inverse();
  
  // Use S now because it is still C^-1 now. 
  CC.add_SymmProduct(val, S, S);
  CC.add_SymmProduct(2.0 * mu, I, I);
  CC.add_OutProduct( kappa, S, S);

  S.scale(kappa * std::log(detF));
  S.AXPY(mu, C);
  S.AXPY(-1.0 * mu, I);
  P.MatMult(F,S);
}


double MaterialModel_StVenant_Kirchhoff_Simo85::get_strain_energy( 
    const Matrix_3x3 &F )
{
  C.MatMultTransposeLeft(F);
  C.AXPY(-1.0, I);
  C.scale(0.5);
  C.MatMult(C,C);
  const double trE2 = C.tr();
  
  detF = F.det();

  return 0.5 * kappa * std::log(detF) * std::log(detF) + mu * trE2;
}

// EOF
