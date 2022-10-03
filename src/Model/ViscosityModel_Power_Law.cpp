#include "ViscosityModel_Power_Law.hpp"

ViscosityModel_Power_Law::ViscosityModel_Power_Law( const double &in_m,
                                                    const double &in_n )
: m( in_m ), n( in_n )
{}

ViscosityModel_Power_Law::ViscosityModel_Power_Law( const char * const &fname )
{
	hid_t h5file = H5Fopen( fname, H5F_ACC_RDONLY, H5P_DEFAULT);

	HDF5_Reader * h5r = new HDF5_Reader( h5file );

	SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
			"Error: ViscosityModel_Power_Law constructor does not match h5 file.\n");

	m      = h5r -> read_doubleScalar("/", "m");
	n      = h5r -> read_doubleScalar("/", "n");

	delete h5r;
	H5Fclose(h5file);
}

ViscosityModel_Power_Law::~ViscosityModel_Power_Law()
{}

void ViscosityModel_Power_Law::print_info() const
{
  SYS_T::commPrint("\t  ViscosityModel_Power_Law:: \n");
  SYS_T::commPrint("\t  Zero Shear Viscosity     m       = %e \n", m);
  SYS_T::commPrint("\t  Power Law Index          n       = %e \n", n);
}

void ViscosityModel_Power_Law::write_hdf5( const char * const &fname ) const
{
  if( SYS_T::get_MPI_rank() == 0 )
  {
    hid_t file_id = H5Fcreate( fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    HDF5_Writer * h5w = new HDF5_Writer( file_id );

    h5w -> write_string("model_name", get_model_name());
    h5w -> write_doubleScalar( "m", m);
    h5w -> write_doubleScalar( "n", n);

    delete h5w;
    H5Fclose(file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

double ViscosityModel_Power_Law::get_mu( const double &D_xx, const double &D_yy,
                                         const double &D_zz, const double &D_yz,
                                         const double &D_xz, const double &D_xy ) const
{
  const SymmMatrix_3x3 D( D_xx, D_yy, D_zz, D_yz, D_xz, D_xy);
  const double DII = std::abs( D.I2() );
  return m * std::pow( 2.0 * std::sqrt( DII ), n - 1.0 );
}

double ViscosityModel_Power_Law::get_mu( const Matrix_3x3 &grad_velo ) const
{
  const SymmMatrix_3x3 D( grad_velo );
  const double DII = std::abs( D.I2() );
  return m * std::pow( 2.0 * std::sqrt( DII ), n - 1.0 );
}

double ViscosityModel_Power_Law::get_dmu_dI1( const double &D_xx,
                          const double &D_yy, const double &D_zz ) const
{
  return 0.0;
}

double ViscosityModel_Power_Law::get_dmu_dI1( const Matrix_3x3 &grad_velo ) const
{
  return 0.0;
}

double ViscosityModel_Power_Law::get_dmu_dI2( const double &D_xx,
      const double &D_yy, const double &D_zz, const double &D_yz,
      const double &D_xz, const double &D_xy ) const
{
  const SymmMatrix_3x3 D( D_xx, D_yy, D_zz, D_yz, D_xz, D_xy);
  const double DII = std::abs( D.I2() );
  const double dmu_dvelo = m * ( n - 1.0) *
                           std::pow( 2.0 * std::sqrt( DII ), n - 2.0 ) /
                           std::sqrt( DII );
  return dmu_dvelo;
}

double ViscosityModel_Power_Law::get_dmu_dI2( const Matrix_3x3 &grad_velo ) const
{
  const SymmMatrix_3x3 D( grad_velo );
  const double DII = std::abs( D.I2() );
  const double dmu_dvelo = m * ( n - 1.0) *
                           std::pow( 2.0 * std::sqrt( DII ), n - 2.0 ) /
                           std::sqrt( DII );
  return dmu_dvelo;
}

double ViscosityModel_Power_Law::get_dmu_dI3( const double &D_xx,
      const double &D_yy, const double &D_zz, const double &D_yz,
      const double &D_xz, const double &D_xy ) const
{
  return 0.0;
}

double ViscosityModel_Power_Law::get_dmu_dI3( const Matrix_3x3 &grad_velo ) const
{
  return 0.0;
}

