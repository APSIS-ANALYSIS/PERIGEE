#include "MaterialModel_NeoHookean_ST91.hpp"

MaterialModel_NeoHookean_ST91::MaterialModel_NeoHookean_ST91( 
    const double &in_rho, const double &in_E, const double &in_nu )
: rho0( in_rho ), E( in_E ), nu( in_nu ), 
  lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  pt33( 1.0 / 3.0 ), pt67( 2.0 / 3.0 ), mpt67( -1.0 * pt67 ),  
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{}


MaterialModel_NeoHookean_ST91::MaterialModel_NeoHookean_ST91(
        const char * const &fname )
:pt33( 1.0 / 3.0 ), pt67( 2.0 / 3.0 ), mpt67( -1.0 * pt67 ),
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{
  hid_t h5file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * h5r = new HDF5_Reader( h5file );

  SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
     "Error: MaterialModel_NeoHookean_ST91 constructor does not match h5 file.\n" );

  rho0   = h5r -> read_doubleScalar("/", "rho0");
  E      = h5r -> read_doubleScalar("/", "E");
  nu     = h5r -> read_doubleScalar("/", "nu");
  lambda = h5r -> read_doubleScalar("/", "lambda");
  mu     = h5r -> read_doubleScalar("/", "mu");
  kappa  = h5r -> read_doubleScalar("/", "kappa");

  delete h5r; H5Fclose(h5file);
}


MaterialModel_NeoHookean_ST91::~MaterialModel_NeoHookean_ST91()
{}

void MaterialModel_NeoHookean_ST91::print_info() const
{
  SYS_T::commPrint("\t  MaterialModel_NeoHookean_ST91: \n");
  SYS_T::commPrint("\t  Young's Modulus E  = %e \n", E);
  SYS_T::commPrint("\t  Possion's ratio nu = %e \n", nu);
  SYS_T::commPrint("\t  Lame coeff lambda  = %e \n", lambda);
  SYS_T::commPrint("\t  Shear modulus mu   = %e \n", mu);
  SYS_T::commPrint("\t  Bulk modulus kappa = %e \n", kappa);
  SYS_T::commPrint("\t  Ref. density rho0  = %e \n", rho0);
}

void MaterialModel_NeoHookean_ST91::write_hdf5( const char * const &fname ) const
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


void MaterialModel_NeoHookean_ST91::get_PK( const Matrix_3x3 &F, 
    Matrix_3x3 &P, Matrix_3x3 &S )
{
  C.MatMultTransposeLeft(F);
  Cinv.copy(C); Cinv.inverse();

  trC = C.tr();
  detF = F.det();
  detF2 = detF * detF;
  detFm0d67 = std::pow(detF, mpt67);

  S.copy(Cinv); S.scale( 0.5 * kappa * (detF2 - 1.0) );
  S.AXPY( mu * detFm0d67, I);
  S.AXPY( (-1.0) * mu * detFm0d67 * pt33 * trC , Cinv );
  P.MatMult(F,S);
}


void MaterialModel_NeoHookean_ST91::get_PK_Stiffness( const Matrix_3x3 &F, 
    Matrix_3x3 &P, Matrix_3x3 &S, Tensor4_3D &Stiffness ) 
{
  C.MatMultTransposeLeft(F);
  Cinv.copy(C); Cinv.inverse();

  trC = C.tr();
  detF = F.det();
  detF2 = detF * detF;
  detFm0d67 = pow(detF, mpt67);

  S.copy(Cinv); S.scale( 0.5 * kappa * (detF2 - 1.0) );
  S.AXPY( mu * detFm0d67, I);
  S.AXPY( (-1.0) * mu * detFm0d67 * pt33 * trC , Cinv );
  P.MatMult(F,S);

  Stiffness.gen_zero();
  const double val1 = 2.0 * pt33 * pt33 * mu * detFm0d67 * trC + kappa * detF2;
  const double val2 = 2.0 * pt33 * mu * detFm0d67 * trC - kappa * (detF2 - 1.0);
  const double val3 = (-2.0) * pt33 * mu * detFm0d67;
  Stiffness.add_OutProduct(val1, Cinv, Cinv);
  Stiffness.add_SymmProduct(val2, Cinv, Cinv);
  Stiffness.add_OutProduct(val3, I, Cinv);
  Stiffness.add_OutProduct(val3, Cinv, I);
}


double MaterialModel_NeoHookean_ST91::get_strain_energy( const Matrix_3x3 &F )
{
  C.MatMultTransposeLeft(F);
  trC = C.tr();
  detF = F.det();
  detF2 = detF * detF;
  detFm0d67 = std::pow(detF, mpt67);
  double lnJ = std::log(detF);
  return 0.5 * mu * (detFm0d67 * trC - 3.0) + 0.5 * kappa *  (0.5*(detF2-1.0) -lnJ);
}


// EOF
