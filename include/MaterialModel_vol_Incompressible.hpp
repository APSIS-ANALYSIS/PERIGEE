#ifndef MATERIALMODEL_VOL_INCOMPRESSIBLE_HPP
#define MATERIALMODEL_VOL_INCOMPRESSIBLE_HPP
// ============================================================================
// MaterialModel_vol_Incompressible.hpp
// 
// Volumetric model for fully incompressible materials.
//
// Note: The Helmholtz free energy is undefined.
//
// Date: Aug. 12 2024
// Author: Ju Liu
// Contact: liujuy@gmail.com
// ============================================================================
#include "IMaterialModel_vol.hpp"

class MaterialModel_vol_Incompressible : public IMaterialModel_vol
{
  public:
    MaterialModel_vol_Incompressible(const double &in_rho_0) : rho_0(in_rho_0) {};

    virtual ~MaterialModel_vol_Incompressible() = default;

    virtual void print_info() const
    {
      SYS_T::commPrint("\t  MaterialModel_vol_Incompressible: \n");
      SYS_T::commPrint("\t  parameter rho_0 = %e \n", rho_0);
    }

    virtual std::string get_model_name() const
    {
      return std::string("Incompressible");
    }

    virtual double get_Gibbs_energy( const double &p ) const {return p;}

    virtual bool is_Gibbs_supported() const {return true;}

    virtual bool is_Helmholtz_supported() const {return false;}

    virtual double get_rho_0() const {return rho_0;}

    virtual double get_rho( const double &p ) const {return rho_0;}

    virtual double get_drho_dp( const double &p ) const {return 0.0;}

    virtual double get_beta( const double &p ) const {return 0.0;}

    virtual double get_dbeta_dp( const double &p ) const {return 0.0;}

  private:
    const double rho_0;
};

#endif
