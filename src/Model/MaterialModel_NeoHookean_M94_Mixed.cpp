#include "MaterialModel_NeoHookean_M94_Mixed.hpp"

MaterialModel_NeoHookean_M94_Mixed::MaterialModel_NeoHookean_M94_Mixed(
    const double &in_rho0, const double &in_E, const double &in_nu )
: rho0(in_rho0), E( in_E ), nu( in_nu ), 
  lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  pt33( 1.0 / 3.0 ), mpt67( -2.0 * pt33 ),
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{}

MaterialModel_NeoHookean_M94_Mixed::MaterialModel_NeoHookean_M94_Mixed(
        const char * const &fname )
: pt33( 1.0 / 3.0 ), mpt67( -2.0 * pt33 ),
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{

hid_t h5file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * h5r = new HDF5_Reader( h5file );

  SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
     "Error: MaterialModel_NeoHookean_M94_Mixed constructor does not match h5 file.\n" );

  rho0   = h5r -> read_doubleScalar("/", "rho0");
  E      = h5r -> read_doubleScalar("/", "E");
  nu     = h5r -> read_doubleScalar("/", "nu");
  lambda = h5r -> read_doubleScalar("/", "lambda");
  mu     = h5r -> read_doubleScalar("/", "mu");
  kappa  = h5r -> read_doubleScalar("/", "kappa");

  delete h5r; H5Fclose(h5file);
}

MaterialModel_NeoHookean_M94_Mixed::~MaterialModel_NeoHookean_M94_Mixed()
{}

void MaterialModel_NeoHookean_M94_Mixed::print_info() const
{
  SYS_T::commPrint("\t  MaterialModel_NeoHookean_M94_Mixed: \n");
  SYS_T::commPrint("\t  Young's Modulus E  = %e \n", E);
  SYS_T::commPrint("\t  Possion's ratio nu = %e \n", nu);
  SYS_T::commPrint("\t  Lame coeff lambda  = %e \n", lambda);
  SYS_T::commPrint("\t  Shear modulus mu   = %e \n", mu);
  SYS_T::commPrint("\t  Bulk modulus kappa = %e \n", kappa);
  SYS_T::commPrint("\t  Ref. density rho0  = %e \n", rho0);
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
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S ) const
{
  Matrix_3x3 C; C.MatMultTransposeLeft(F);
  Tensor2_3D Cinv = inverse(C);

  const double detF = F.det();
  const double detFm0d67 = std::pow(detF, mpt67);

  S.copy(Cinv); S.scale( (-1.0) * mu * detFm0d67 * pt33 * C.tr() );
  S.AXPY( mu * detFm0d67, I);
  P.MatMult(F,S);
}

void MaterialModel_NeoHookean_M94_Mixed::get_PK_Stiffness(
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S, Tensor4_3D &CC ) const
{
  Matrix_3x3 C; C.MatMultTransposeLeft(F);
  Tensor2_3D Cinv = inverse(C);

  const double detF = F.det();
  const double detFm0d67 = std::pow(detF, mpt67);

  S.copy(Cinv); S.scale( (-1.0) * mu * detFm0d67 * pt33 * C.tr() );
  S.AXPY( mu * detFm0d67, I);
  P.MatMult(F,S);

  CC.gen_zero();
  const double val1 = 2.0 * pt33 * pt33 * mu * detFm0d67 * C.tr();
  const double val2 = 2.0 * pt33 * mu * detFm0d67 * C.tr();
  const double val3 = mpt67 * mu * detFm0d67;
  CC.add_OutProduct(val1, Cinv, Cinv);
  CC.add_SymmProduct(val2, Cinv, Cinv);
  CC.add_OutProduct(val3, I, Cinv);
  CC.add_OutProduct(val3, Cinv, I);
}

double MaterialModel_NeoHookean_M94_Mixed::get_strain_energy( 
    const Matrix_3x3 &F ) const
{
  Matrix_3x3 C; C.MatMultTransposeLeft(F);
  return 0.5 * mu * ( std::pow(F.det(), mpt67) * C.tr() - 3.0);
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
