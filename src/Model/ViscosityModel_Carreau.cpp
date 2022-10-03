#include "ViscosityModel_Carreau.hpp"

ViscosityModel_Carreau::ViscosityModel_Carreau( const double &in_mu_inf,
const double &in_mu_0, const double &in_lambda, const double &in_n )
: mu_inf( in_mu_inf ), mu_0( in_mu_0 ), lambda( in_lambda ), n( in_n )
{}

ViscosityModel_Carreau::ViscosityModel_Carreau( const char * const &fname )
{
	hid_t h5file = H5Fopen( fname, H5F_ACC_RDONLY, H5P_DEFAULT);

	HDF5_Reader * h5r = new HDF5_Reader( h5file );

	SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
			"Error: ViscosityModel_Carreau constructor does not match h5 file.\n");

	mu_inf = h5r -> read_doubleScalar("/", "mu_inf");
	mu_0   = h5r -> read_doubleScalar("/", "mu_0");
	lambda = h5r -> read_doubleScalar("/", "lambda");
	n      = h5r -> read_doubleScalar("/", "n");

	delete h5r;
	H5Fclose(h5file);
}

ViscosityModel_Carreau::~ViscosityModel_Carreau()
{}

void ViscosityModel_Carreau::print_info() const
{
  SYS_T::commPrint("\t  ViscosityModel_Carreau:: \n");
  SYS_T::commPrint("\t  Infinite Shear Viscosity mu_inf  = %e \n", mu_inf);
  SYS_T::commPrint("\t  Zero Shear Viscosity     mu_0    = %e \n", mu_0);
  SYS_T::commPrint("\t  Time Constant            lambda  = %e \n", lambda);
  SYS_T::commPrint("\t  Power Law Index          n       = %e \n", n);
}

void ViscosityModel_Carreau::write_hdf5( const char * const &fname ) const
{
  if( SYS_T::get_MPI_rank() == 0 )
  {
    hid_t file_id = H5Fcreate( fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    HDF5_Writer * h5w = new HDF5_Writer( file_id );

    h5w -> write_string("model_name", get_model_name());
    h5w -> write_doubleScalar( "mu_inf", mu_inf);
    h5w -> write_doubleScalar( "mu_0", mu_0);
    h5w -> write_doubleScalar( "lambda", lambda);
    h5w -> write_doubleScalar( "n", n);

    delete h5w;
    H5Fclose(file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

double ViscosityModel_Carreau::get_mu( const double &D_xx, const double &D_yy,
                                       const double &D_zz, const double &D_yz,
                                       const double &D_xz, const double &D_xy ) const
{
  const SymmMatrix_3x3 D( D_xx, D_yy, D_zz, D_yz, D_xz, D_xy);
  const double DII = std::abs( D.I2() );
  const double pow_base = 1.0 + lambda * lambda * 4.0 * DII;
  return mu_inf + ( mu_0 - mu_inf ) * std::pow( pow_base, (n - 1.0) / 2.0 );
}

double ViscosityModel_Carreau::get_mu( const Matrix_3x3 &grad_velo ) const
{
  const SymmMatrix_3x3 D( grad_velo );
  const double DII = std::abs( D.I2() );
  const double pow_base = 1.0 + lambda * lambda * 4.0 * DII;
  return mu_inf + ( mu_0 - mu_inf ) * std::pow( pow_base, (n - 1.0) / 2.0 );
}

double ViscosityModel_Carreau::get_dmu_dI1( const double &D_xx,
                        const double &D_yy, const double &D_zz ) const
{
  return 0.0;
}

double ViscosityModel_Carreau::get_dmu_dI1( const Matrix_3x3 &grad_velo ) const
{
  return 0.0;
}

double ViscosityModel_Carreau::get_dmu_dI2( const double &D_xx,
    const double &D_yy, const double &D_zz, const double &D_yz,
    const double &D_xz, const double &D_xy ) const
{
  const SymmMatrix_3x3 D( D_xx, D_yy, D_zz, D_yz, D_xz, D_xy);
  const double DII = std::abs( D.I2() );
  const double pow_base = 1.0 + lambda * lambda * 4.0 * DII;
  const double dmu_dvelo = ( mu_0 - mu_inf ) * ( n - 1.0 ) * lambda * 2.0 *
                           std::pow( pow_base, ( n - 3.0 ) / 2.0);
  return dmu_dvelo;
}

double ViscosityModel_Carreau::get_dmu_dI2( const Matrix_3x3 &grad_velo ) const
{
  const SymmMatrix_3x3 D( grad_velo );
  const double DII = std::abs( D.I2() );
  const double pow_base = 1.0 + lambda * lambda * 4.0 * DII;
  const double dmu_dvelo = ( mu_0 - mu_inf ) * ( n - 1.0 ) * lambda * 2.0 *
                           std::pow( pow_base, ( n - 3.0 ) / 2.0);
  return dmu_dvelo;
}

double ViscosityModel_Carreau::get_dmu_dI3( const double &D_xx,
         const double &D_yy, const double &D_zz, const double &D_yz,
         const double &D_xz, const double &D_xy ) const
{
  return 0.0;
}

double ViscosityModel_Carreau::get_dmu_dI3( const Matrix_3x3 &grad_velo ) const
{
  return 0.0;
}

