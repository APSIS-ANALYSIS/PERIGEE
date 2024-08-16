#ifndef MATERIALMODEL_VOL_ST91_HPP
#define MATERIALMODEL_VOL_ST91_HPP
// ============================================================================
// MaterialModel_vol_ST91.hpp
// 
// Volumetric model from J.C.Simo & R.L.Taylor. CMAME 85:273-310, 1991. 
//
// Date: Aug. 16 2024
// Author: Ju Liu, Jiawei Luo
// Contact: liujuy@gmail.com
// ============================================================================
#include "IMaterialModel_vol.hpp"

class MaterialModel_vol_ST91 : public IMaterialModel_vol
{
  public:
    MaterialModel_vol_ST91(const double &in_rho_0, const double &in_kappa)
      : rho_0(in_rho_0), kappa( in_kappa ) {};

    virtual ~MaterialModel_vol_ST91() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\t  MaterialModel_vol_ST91: \n");
      SYS_T::commPrint("\t  parameter rho_0 = %e \n", rho_0);
      SYS_T::commPrint("\t  parameter kappa = %e \n", kappa);
    }

    virtual std::string get_model_name() const
    {
      return std::string("Simo-Taylor-1991");
    }

    virtual bool is_Gibbs_supported() const {return true;}

    virtual bool is_Helmholtz_supported() const {return true;}

    virtual double get_elastic_kappa() const {return kappa;}

    virtual double get_Gibbs_energy( const double &p ) const 
    {return (p * std::sqrt(p * p + kappa * kappa) - p*p) / (2.0*kappa) - kappa / 2.0 * std::log((std::sqrt(p*p + kappa*kappa) - p) / kappa);}

    virtual double get_Helmholtz_energy( const double &J ) const
    {return kappa / 4.0 * (J*J - 2.0 * std::log(J) - 1.0);}

    virtual double get_rho_0() const 
    {return rho_0;}

    virtual double get_rho( const double &p ) const 
    {return rho_0 * (std::sqrt(p*p + kappa*kappa) + p ) / kappa;}

    virtual double get_drho_dp( const double &p ) const 
    {return (rho_0 * p) / (kappa * std::sqrt(p*p + kappa*kappa)) + rho_0 / kappa;}

    virtual double get_beta( const double &p ) const 
    {return 1.0/std::sqrt(p*p+kappa*kappa);}

    virtual double get_dbeta_dp( const double &p ) const 
    {return (-1.0 * p) / ( std::pow(p*p+kappa*kappa, 1.5) );}

    virtual double get_vol_stress( const double &J ) const
    {return kappa / 2.0 *( J - (1.0/J) );}

    virtual double get_dvol_stress_dJ( const double &J ) const
    {return kappa/2.0 * (1.0 + 1.0 / (J*J) );}

  private:
    const double rho_0, kappa;
};

#endif
