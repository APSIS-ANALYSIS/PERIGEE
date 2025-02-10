#ifndef VISCOSITYMODEL_POWER_LAW_HPP
#define VISCOSITYMODEL_POWER_LAW_HPP
// ============================================================================
// ViscosityModel_Power_Law.hpp
//
// Interface for Power Law non-Newtonian model
// ============================================================================
#include "IViscosityModel.hpp"

class ViscosityModel_Power_Law final : public IViscosityModel
{
  public:
    ViscosityModel_Power_Law( const double &in_m_cons, const double &in_n_pli,
        const double &in_mu_max, const double &in_mu_min ) : m_cons( in_m_cons ),
    n_pli( in_n_pli ), mu_max( in_mu_max ), mu_min( in_mu_min ) {};

    ~ViscosityModel_Power_Law() override = default;

    void print_info() const override
    {
      SYS_T::commPrint("\t  ViscosityModel_Power_Law:: \n");
      SYS_T::commPrint("\t  Consistency       m_cons    = %e \n", m_cons);
      SYS_T::commPrint("\t  Power Law Index   n_pli     = %e \n", n_pli);
      SYS_T::commPrint("\t  Power Law maximum viscosity = %e \n", mu_max);
      SYS_T::commPrint("\t  Power Law minimum viscosity = %e \n", mu_min);
    }

    std::string get_model_name() const override {return "Power-Law";}

    double get_mu( const SymmTensor2_3D &strain_rate ) const override
    {
      const double strain_rate_II = strain_rate.MatContraction( strain_rate );
      const double temp_mu = m_cons * std::pow( std::sqrt( 2.0 * strain_rate_II ), n_pli - 1.0 );

      return (temp_mu <= mu_min) ? mu_min : (temp_mu > mu_max ? mu_max : temp_mu);
    }

    double get_dmu_dI1( const SymmTensor2_3D &strain_rate ) const override
    {
      return 0.0;
    }

    double get_dmu_dI2( const SymmTensor2_3D &strain_rate ) const override
    {
      const double strain_rate_II = strain_rate.MatContraction( strain_rate );

      const double dmu_dvelo = m_cons * ( n_pli - 1.0) *
        std::pow( std::sqrt( 2.0 * strain_rate_II ), n_pli - 2.0 ) /
        std::sqrt( 2.0 * strain_rate_II );

      const double temp_mu = m_cons * std::pow( std::sqrt( 2.0 * strain_rate_II ), n_pli - 1.0 );

      return (temp_mu <= mu_min || temp_mu > mu_max) ? 0.0 : dmu_dvelo;
    }

    double get_dmu_dI3( const SymmTensor2_3D &strain_rate ) const override
    {
      return 0.0;
    }

  private:
    // ------------------------------------------------------------------------
    // m_cons : consistency ( Pa * s ^ n ) when n = 1, it is the same as viscosity
    // n_pli  : an index shows the dependency of viscosity on shear rate
    // ------------------------------------------------------------------------
    const double m_cons, n_pli;

    // ------------------------------------------------------------------------
    // mu_max : maximum viscosity
    // mu_min : minimum viscosity
    // ------------------------------------------------------------------------
    const double mu_max, mu_min;
    
    ViscosityModel_Power_Law() = delete;
};

#endif
