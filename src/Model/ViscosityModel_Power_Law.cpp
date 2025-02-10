#include "ViscosityModel_Power_Law.hpp"

ViscosityModel_Power_Law::ViscosityModel_Power_Law( const double &in_m_cons,
                                                    const double &in_n_pli,
                                                    const double &in_mu_max,
                                                    const double &in_mu_min )
: m_cons( in_m_cons ), n_pli( in_n_pli ), mu_max( in_mu_max ), mu_min( in_mu_min )
{}

void ViscosityModel_Power_Law::print_info() const
{
  SYS_T::commPrint("\t  ViscosityModel_Power_Law:: \n");
  SYS_T::commPrint("\t  Consistency       m_cons    = %e \n", m_cons);
  SYS_T::commPrint("\t  Power Law Index   n_pli     = %e \n", n_pli);
  SYS_T::commPrint("\t  Power Law maximum viscosity = %e \n", mu_max);
  SYS_T::commPrint("\t  Power Law minimum viscosity = %e \n", mu_min);
}

double ViscosityModel_Power_Law::get_mu( const SymmTensor2_3D &strain_rate ) const
{
  const double strain_rate_II = strain_rate.MatContraction( strain_rate );
  const double temp_mu = m_cons * std::pow( std::sqrt( 2.0 * strain_rate_II ), n_pli - 1.0 );
  
  return (temp_mu <= mu_min) ? mu_min : (temp_mu > mu_max ? mu_max : temp_mu);
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

  const double temp_mu = m_cons * std::pow( std::sqrt( 2.0 * strain_rate_II ), n_pli - 1.0 );

  return (temp_mu <= mu_min || temp_mu > mu_max) ? 0.0 : dmu_dvelo; 
}

double ViscosityModel_Power_Law::get_dmu_dI3( const SymmTensor2_3D &strain_rate ) const
{
  return 0.0;    
}

// EOF
