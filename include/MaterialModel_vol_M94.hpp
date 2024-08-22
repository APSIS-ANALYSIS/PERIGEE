#ifndef MATERIALMODEL_VOL_M94_HPP
#define MATERIALMODEL_VOL_M94_HPP
// ============================================================================
// MaterialModel_vol_M94.hpp
// 
// Volumetric model from C.Miehe. IJNME 37:1981-2004, 1994. 
//
// Date: Aug. 14 2024
// Author: Ju Liu
// Contact: liujuy@gmail.com
// ============================================================================
#include "IMaterialModel_vol.hpp"

class MaterialModel_vol_M94 : public IMaterialModel_vol
{
  public:
    MaterialModel_vol_M94(const double &in_rho_0, const double &in_kappa)
      : rho_0( in_rho_0 ), kappa( in_kappa ) {};

    virtual ~MaterialModel_vol_M94() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\t  MaterialModel_vol_M94: \n");
      SYS_T::commPrint("\t  parameter rho_0 = %e \n", rho_0);
      SYS_T::commPrint("\t  parameter kappa = %e \n", kappa);
    }

    virtual std::string get_model_name() const
    {
      return std::string("Miehe-1994");
    }

    virtual bool is_Gibbs_supported() const {return true;}

    virtual bool is_Helmholtz_supported() const {return true;}

    virtual double get_elastic_kappa() const {return kappa;}

    virtual double get_Gibbs_energy( const double &p ) const 
    {return kappa * ( std::log(p+kappa)-std::log(kappa) );}

    virtual double get_Helmholtz_energy( const double &J ) const
    {return kappa * (J - std::log(J) - 1.0);}

    virtual double get_rho_0() const 
    {return rho_0;}

    virtual double get_rho( const double &p ) const 
    {return rho_0 * ( 1.0 + (p/kappa) );}

    virtual double get_drho_dp( const double &p ) const 
    {return rho_0 / kappa;}

    virtual double get_beta( const double &p ) const 
    {return 1.0/(p+kappa);}

    virtual double get_dbeta_dp( const double &p ) const 
    {return (-1.0) / ( (p+kappa) * (p+kappa) );}

    virtual double get_vol_stress( const double &J ) const
    {return kappa*( 1.0 - (1.0/J) );}

    virtual double get_dvol_stress_dJ( const double &J ) const
    {return kappa/(J*J);}

  private:
    const double rho_0, kappa;
};

#endif
