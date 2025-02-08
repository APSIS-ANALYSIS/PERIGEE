#include "ViscosityModel_Power_Law.hpp"

ViscosityModel_Power_Law::ViscosityModel_Power_Law( const double &in_m_cons,
                                                    const double &in_n_pli )
: m_cons( in_m_cons ), n_pli( in_n_pli )
{}

ViscosityModel_Power_Law::ViscosityModel_Power_Law( const char * const &fname )
{
	hid_t h5file = H5Fopen( fname, H5F_ACC_RDONLY, H5P_DEFAULT);

	HDF5_Reader * h5r = new HDF5_Reader( h5file );

	SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
			"Error: ViscosityModel_Power_Law constructor does not match h5 file.\n");

	m_cons = h5r -> read_doubleScalar("/", "m_cons");
	n_pli  = h5r -> read_doubleScalar("/", "n_pli");

	delete h5r;
	H5Fclose(h5file);
}

ViscosityModel_Power_Law::~ViscosityModel_Power_Law()
{}

void ViscosityModel_Power_Law::print_info() const
{
  SYS_T::commPrint("\t  ViscosityModel_Power_Law:: \n");
  SYS_T::commPrint("\t  Consistency       m_cons  = %e \n", m_cons);
  SYS_T::commPrint("\t  Power Law Index   n_pli   = %e \n", n_pli);
}

double ViscosityModel_Power_Law::get_mu( const SymmTensor2_3D &strain_rate ) const
{
  const double strain_rate_II = strain_rate.MatContraction( strain_rate );
  return m_cons * std::pow( std::sqrt( 2.0 * strain_rate_II ), n_pli - 1.0 ); 
}

double ViscosityModel_Power_Law::get_dmu_dI1( const SymmTensor2_3D &strain_rate ) const
{
  return 0.0;
}

double ViscosityModel_Power_Law::get_dmu_dI2( const SymmTensor2_3D &strain_rate ) const
{
  const double strain_rate_II = strain_rate.MatContraction( strain_rate );
  const double dmu_dvelo = m_cons * ( n_pli - 1.0) * 
                           std::pow( std::sqrt( 2.0 * strain_rate_II ), n_pli - 2.0 ) /
                           std::sqrt( 2.0 * strain_rate_II );
  return dmu_dvelo;   
}

double ViscosityModel_Power_Law::get_dmu_dI3( const SymmTensor2_3D &strain_rate ) const
{
  return 0.0;    
}

// EOF
