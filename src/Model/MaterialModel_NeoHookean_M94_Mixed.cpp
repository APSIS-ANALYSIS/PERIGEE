#include "MaterialModel_NeoHookean_M94_Mixed.hpp"

MaterialModel_NeoHookean_M94_Mixed::MaterialModel_NeoHookean_M94_Mixed(
    const double &in_rho0, const double &in_E, const double &in_nu )
: rho0(in_rho0), E( in_E ), nu( in_nu ), 
  lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  pt33( 1.0 / 3.0 ), mpt67( -2.0 * pt33 ),
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{}


MaterialModel_NeoHookean_M94_Mixed::~MaterialModel_NeoHookean_M94_Mixed()
{}


void MaterialModel_NeoHookean_M94_Mixed::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  MaterialModel_NeoHookean_M94_Mixed: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Young's Modulus E  = %e \n", E);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Possion's ratio nu = %e \n", nu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Lame coeff lambda  = %e \n", lambda);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Shear modulus mu   = %e \n", mu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Bulk modulus kappa = %e \n", kappa);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Ref. density rho0  = %e \n", rho0);
}

void MaterialModel_NeoHookean_M94_Mixed::write_hdf5( const char * const &fname ) const
{
  if( SYS_T::get_MPI_rank() == 0 )
  {
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * h5w = new HDF5_Writer(file_id);

    h5w -> write_string("model_name", get_model_name());
    h5w -> write_doubleScalar( "rho0", rho0 );
    h5w -> write_doubleScalar("E", E);
    h5w -> write_doubleScalar("nu", nu);
    h5w -> write_doubleScalar("lambda", lambda);
    h5w -> write_doubleScalar("mu", mu);
    h5w -> write_doubleScalar("kappa", kappa);

    delete h5w; H5Fclose(file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

void MaterialModel_NeoHookean_M94_Mixed::get_PK(
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S )
{
  C.MatMultTransposeLeft(F);
  Cinv.copy(C); Cinv.inverse();

  trC = C.tr();
  detF = F.det();
  detFm0d67 = std::pow(detF, mpt67);

  S.copy(Cinv); S.scale( (-1.0) * mu * detFm0d67 * pt33 * trC );
  S.AXPY( mu * detFm0d67, I);
  P.MatMult(F,S);
}


void MaterialModel_NeoHookean_M94_Mixed::get_PK_Stiffness(
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S, Tensor4_3D &CC )
{
  C.MatMultTransposeLeft(F);
  Cinv.copy(C); Cinv.inverse();

  trC = C.tr();
  detF = F.det();
  detFm0d67 = std::pow(detF, mpt67);

  S.copy(Cinv); S.scale( (-1.0) * mu * detFm0d67 * pt33 * trC );
  S.AXPY( mu * detFm0d67, I);
  P.MatMult(F,S);

  CC.gen_zero();
  const double val1 = 2.0 * pt33 * pt33 * mu * detFm0d67 * trC;
  const double val2 = 2.0 * pt33 * mu * detFm0d67 * trC;
  const double val3 = mpt67 * mu * detFm0d67;
  CC.add_OutProduct(val1, Cinv, Cinv);
  CC.add_SymmProduct(val2, Cinv, Cinv);
  CC.add_OutProduct(val3, I, Cinv);
  CC.add_OutProduct(val3, Cinv, I);
}


double MaterialModel_NeoHookean_M94_Mixed::get_strain_energy( 
    const Matrix_3x3 &F )
{
  C.MatMultTransposeLeft(F);
  trC = C.tr();
  detF = F.det();
  detFm0d67 = std::pow(detF, mpt67);

  return 0.5 * mu * (detFm0d67 * trC - 3.0);
}


double MaterialModel_NeoHookean_M94_Mixed::get_rho( const double &p ) const
{
  const double pk = p / kappa;

  return rho0 * (1.0 + pk);
}


double MaterialModel_NeoHookean_M94_Mixed::get_drho_dp( const double &p ) const
{
  return rho0 / kappa;
}


double MaterialModel_NeoHookean_M94_Mixed::get_beta( const double &p ) const
{
  return 1.0 / (p + kappa);
}


double MaterialModel_NeoHookean_M94_Mixed::get_dbeta_dp( const double &p ) const
{
  return (-1.0) / ( (p+kappa) * (p+kappa) );
}


// EOF
