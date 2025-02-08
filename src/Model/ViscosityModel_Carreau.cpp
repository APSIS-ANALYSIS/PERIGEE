#include "ViscosityModel_Carreau.hpp"

ViscosityModel_Carreau::ViscosityModel_Carreau( const double &in_mu_inf,
const double &in_mu_0, const double &in_lambda, const double &in_n_pli )
: mu_inf( in_mu_inf ), mu_0( in_mu_0 ), lambda( in_lambda ), n_pli( in_n_pli )
{}

void ViscosityModel_Carreau::print_info() const
{
  SYS_T::commPrint("\t  ViscosityModel_Carreau:: \n");
  SYS_T::commPrint("\t  Infinite Shear Viscosity mu_inf  = %e \n", mu_inf);
  SYS_T::commPrint("\t  Zero Shear Viscosity     mu_0    = %e \n", mu_0);
  SYS_T::commPrint("\t  Time Constant            lambda  = %e \n", lambda);
  SYS_T::commPrint("\t  Power Law Index          n_pli   = %e \n", n_pli);
}

double ViscosityModel_Carreau::get_mu( const SymmTensor2_3D &strain_rate ) const
{
  const double strain_rate_II = strain_rate.MatContraction( strain_rate );
  const double pow_base = 1.0 + lambda * lambda * 2.0 * strain_rate_II;
  return mu_inf + ( mu_0 - mu_inf ) * std::pow( pow_base, (n_pli - 1.0) * 0.5 );
}

double ViscosityModel_Carreau::get_dmu_dI1( const SymmTensor2_3D &strain_rate ) const
{
  return 0.0;
}

double ViscosityModel_Carreau::get_dmu_dI2( const SymmTensor2_3D &strain_rate ) const
{
  const double strain_rate_II = strain_rate.MatContraction( strain_rate );
  const double pow_base = 1.0 + lambda * lambda * 2.0 * strain_rate_II;
  const double dmu_dvelo = ( mu_0 - mu_inf ) * ( n_pli - 1.0 ) * lambda * lambda 
                            * std::pow( pow_base, ( n_pli - 3.0 ) * 0.5 );
  return dmu_dvelo;
}

double ViscosityModel_Carreau::get_dmu_dI3( const SymmTensor2_3D &strain_rate ) const
{
  return 0.0;
}

// EOF
