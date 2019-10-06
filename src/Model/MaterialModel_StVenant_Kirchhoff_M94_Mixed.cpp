#include "MaterialModel_StVenant_Kirchhoff_M94_Mixed.hpp"

MaterialModel_StVenant_Kirchhoff_M94_Mixed::MaterialModel_StVenant_Kirchhoff_M94_Mixed( const double &in_E, const double &in_nu )
: rho0(1.0), E( in_E ), nu( in_nu ), lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  pt33( 1.0 / 3.0 ), mpt67( -2.0 * pt33 ),
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{}


MaterialModel_StVenant_Kirchhoff_M94_Mixed::MaterialModel_StVenant_Kirchhoff_M94_Mixed( const double &in_rho0, const double &in_E, const double &in_nu )
: rho0( in_rho0 ), E( in_E ), nu( in_nu ), 
  lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  pt33( 1.0 / 3.0 ), mpt67( -2.0 * pt33 ),
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{}


MaterialModel_StVenant_Kirchhoff_M94_Mixed::~MaterialModel_StVenant_Kirchhoff_M94_Mixed()
{}


void MaterialModel_StVenant_Kirchhoff_M94_Mixed::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  MaterialModel_StVenant_Kirchhoff_M94_Mixed: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Young's Modulus E  = %e \n", E);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Possion's ratio nu = %e \n", nu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Lame coeff lambda  = %e \n", lambda);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Shear modulus mu   = %e \n", mu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Bulk modulus kappa = %e \n", kappa);
  PetscPrintf(PETSC_COMM_WORLD, "\t  rho_0              = %e \n", rho0);
}


void MaterialModel_StVenant_Kirchhoff_M94_Mixed::get_PK( 
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S )
{
  C.MatMultTransposeLeft(F);
  Cinv.copy(C); Cinv.inverse();
  trC = C.tr();
  detF = F.det();
  detFm0d67 = std::pow(detF, mpt67);

  const double mpt33CdC = (-1.0) * pt33 * C.MatContraction(C);
  PxC.copy(Cinv); PxC.scale(mpt33CdC); PxC.PY(C);
  
  PxI.copy(Cinv); PxI.scale( (-1.0) * pt33 * trC ); PxI.PY(I);

  S.copy(PxC); S.scale(detFm0d67);  S.AXPY( -1.0, PxI );
  S.scale( mu*detFm0d67 );

  P.MatMult(F,S);
}


void MaterialModel_StVenant_Kirchhoff_M94_Mixed::get_PK_Stiffness( 
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S, Tensor4_3D &CC )
{
  C.MatMultTransposeLeft(F);
  Cinv.copy(C); Cinv.inverse();
  trC = C.tr();
  detF = F.det();
  detFm0d67 = std::pow(detF, mpt67);

  const double CdC = C.MatContraction(C);

  const double mpt33CdC = (-1.0) * pt33 * CdC;
  PxC.copy(Cinv); PxC.scale(mpt33CdC); PxC.PY(C);
  
  PxI.copy(Cinv); PxI.scale( (-1.0) * pt33 * trC ); PxI.PY(I);

  S.copy(PxC); S.scale(detFm0d67);  S.AXPY( -1.0, PxI );
  S.scale( mu*detFm0d67 );

  P.MatMult(F,S);

  CC.gen_symm_id();

  CC.add_OutProduct( (-1.0) * pt33, Cinv, C );
  CC.add_OutProduct( (-1.0) * pt33, C, Cinv );

  const double val1 = pt33 * pt33 * CdC;

  CC.add_OutProduct( val1, Cinv, Cinv );

  CC.scale( 2.0 * mu * detFm0d67 * detFm0d67 );

  const double val2 = 2.0 * pt33 * detFm0d67 * mu * ( detFm0d67 * CdC - trC );

  CC.add_SymmProduct(val2, Cinv, Cinv);
  CC.add_OutProduct((-1.0)*pt33*val2, Cinv, Cinv);

  CC.add_OutProduct(mpt67, Cinv, S);
  CC.add_OutProduct(mpt67, S, Cinv);
}


double MaterialModel_StVenant_Kirchhoff_M94_Mixed::get_strain_energy( 
    const Matrix_3x3 &F )
{
  C.MatMultTransposeLeft(F);
  detF = F.det();
  detFm0d67 = std::pow(detF, mpt67);
  
  Matrix_3x3 E( C ); E.scale(0.5 * detFm0d67); E.AXPY(-0.5, I);

  E.MatMult(E,E);

  return mu * E.tr();
}


double MaterialModel_StVenant_Kirchhoff_M94_Mixed::get_rho( 
    const double &p ) const
{
  const double pk = p / kappa;

  return rho0 * (1.0 + pk);
}


double MaterialModel_StVenant_Kirchhoff_M94_Mixed::get_drho_dp( 
    const double &p ) const
{
  return rho0 / kappa;
}


double MaterialModel_StVenant_Kirchhoff_M94_Mixed::get_beta( 
    const double &p ) const
{
  return 1.0 / (p + kappa);
}


double MaterialModel_StVenant_Kirchhoff_M94_Mixed::get_dbeta_dp( 
    const double &p ) const
{
  return (-1.0) / ( (p+kappa) * (p+kappa) );
}


// EOF
