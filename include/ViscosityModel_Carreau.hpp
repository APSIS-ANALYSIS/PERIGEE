#ifndef VISCOSITYMODEL_CARREAU_HPP
#define VISCOSITYMODEL_CARREAU_HPP
// ============================================================================
// ViscosityModel_Carreau.hpp
//
// Interface for Carreau non-Newtonian model
// Ref. Cho and Kensey, Biorheology, 28:241-262, 1991
//
// Auther: Xinhai Yue
// Email: seavegetableyxh@outlook.com
// ============================================================================
#include "IViscosityModel.hpp"

class ViscosityModel_Carreau final : public IViscosityModel
{
  public:
    ViscosityModel_Carreau( const double &in_mu_inf, const double &in_mu_0,
        const double &in_lambda, const double &in_n_pli ) : mu_inf( in_mu_inf ),
    mu_0( in_mu_0 ), lambda( in_lambda ), n_pli( in_n_pli ) {};

    ~ViscosityModel_Carreau() override = default;

    void print_info() const override
    {
      SYS_T::commPrint("\t  ViscosityModel_Carreau:: \n");
      SYS_T::commPrint("\t  Infinite Shear Viscosity mu_inf  = %e \n", mu_inf);
      SYS_T::commPrint("\t  Zero Shear Viscosity     mu_0    = %e \n", mu_0);
      SYS_T::commPrint("\t  Time Constant            lambda  = %e \n", lambda);
      SYS_T::commPrint("\t  Power Law Index          n_pli   = %e \n", n_pli);
    }

    std::string get_model_name() const override {return "Carreau";}

    double get_mu( const SymmTensor2_3D &strain_rate ) const override
    {
      const double strain_rate_II = strain_rate.MatContraction();
      const double pow_base = 1.0 + lambda * lambda * 2.0 * strain_rate_II;
      return mu_inf + ( mu_0 - mu_inf ) * std::pow( pow_base, (n_pli - 1.0) * 0.5 );
    }

    double get_dmu_dI1( const SymmTensor2_3D &strain_rate ) const override 
    {return 0.0;}

    double get_dmu_dI2( const SymmTensor2_3D &strain_rate ) const override
    {
      const double strain_rate_II = strain_rate.MatContraction();
      const double pow_base = 1.0 + lambda * lambda * 2.0 * strain_rate_II;
      const double dmu_dvelo = ( mu_0 - mu_inf ) * ( n_pli - 1.0 ) * lambda * lambda
        * std::pow( pow_base, ( n_pli - 3.0 ) * 0.5 );
      return dmu_dvelo;
    }

    double get_dmu_dI3( const SymmTensor2_3D &strain_rate ) const override
    {return 0.0;}

  private:
    // ----------------------------------------------------------------------------
    // mu_inf : viscosity as shear rate tends to infinity
    // mu_0   : viscosity as shear rate tends to 0
    // lambda : a time constant
    // n_pli  : Power law index, an index shows the dependency of viscosity on 
    //          shear rate
    // ----------------------------------------------------------------------------
    const double mu_inf, mu_0, lambda, n_pli;
    
    ViscosityModel_Carreau() = delete;
};

#endif
