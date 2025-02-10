#ifndef VISCOSITYMODEL_POWER_LAW_HPP
#define VISCOSITYMODEL_POWER_LAW_HPP
// ============================================================================
// ViscosityModel_Power_Law.hpp
//
// Interface for Power Law non-Newtonian model
// ============================================================================
#include "IViscosityModel.hpp"

class ViscosityModel_Power_Law : public IViscosityModel
{
  public:
    ViscosityModel_Power_Law() = delete;
    
    ViscosityModel_Power_Law( const double &in_m_cons, const double &in_n_pli, const double & in_mu_max, const double & in_mu_min );

    virtual ~ViscosityModel_Power_Law() = default;

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "Power Law";
      return mname;
    }
    
    virtual double get_mu( const SymmTensor2_3D &strain_rate ) const;

    virtual double get_dmu_dI1( const SymmTensor2_3D &strain_rate ) const;

    virtual double get_dmu_dI2( const SymmTensor2_3D &strain_rate ) const;

    virtual double get_dmu_dI3( const SymmTensor2_3D &strain_rate ) const;

    private:
    // ------------------------------------------------------------------------
    // m_cons : consistency ( Pa * s ^ n ) when n = 1, it is the same as viscosity
    // n_pli  : an index shows the dependency of viscosity on shear rate
    // ------------------------------------------------------------------------
    double m_cons, n_pli;

    // ------------------------------------------------------------------------
    // mu_max : maximum viscosity
    // mu_min : minimum viscosity
    // ------------------------------------------------------------------------
    const double mu_max, mu_min;
};

#endif
