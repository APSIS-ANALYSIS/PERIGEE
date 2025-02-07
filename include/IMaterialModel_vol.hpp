#ifndef IMATERIALMODEL_VOL_HPP
#define IMATERIALMODEL_VOL_HPP
// ============================================================================
// IMaterialModel_vol.hpp
//
// Interface for material model of the volumetric behavior.
// 
// This is an interface for different volumetric models that are represented in
// terms of Gibbs and/or Helmholtz free energies. (The Helmholtz free energy is
// not supported only for the fully incompressible model; there are a few
// energies that we cannot obtain an analytic form for the Gibss free energy and
// their Gibbs representation is not supported therefore.)
// 
// This class supports the material model that are described using the
// Isochoric-Dialatational split.
//
// Date: Aug. 12 2024
// Author: Ju Liu
// Contact: liujuy@gmail.com
// ============================================================================

class IMaterialModel_vol
{
  public:
    IMaterialModel_vol() = default;

    virtual ~IMaterialModel_vol() = default;

    virtual void print_info() const = 0;

    virtual std::string get_model_name() const
    {
      SYS_T::commPrint("Warning: IMaterialModel_vol::get_model_name() is not implemented. \n");
      return std::string( "unknown" );
    }

    // ------------------------------------------------------------------------
    // Output: The bulk modulus used for defining the volumetric energy
    // ------------------------------------------------------------------------
    virtual double get_elastic_kappa() const
    {
      SYS_T::commPrint("Warning: IMaterialModel_vol::get_elastic_kappa() is not implemented. \n");
      return 0.0;
    }

    // ------------------------------------------------------------------------
    // Input: Pressure p
    // Output: The Gibbs free energy (energy per unit volume), denoted by G
    // ------------------------------------------------------------------------
    virtual double get_Gibbs_energy( const double &p ) const
    {
      SYS_T::commPrint("Warning: IMaterialModel_vol::get_Gibbs_energy(p) is not implemented. \n");
      return 0.0;
    }

    virtual bool is_Gibbs_supported() const = 0;

    // ------------------------------------------------------------------------
    // Input: Volume ration J
    // Output: The Helmholtz free energy (energy per unit volume), denoted by W
    // or by Psi in Holzapfel book
    // ------------------------------------------------------------------------
    virtual double get_Helmholtz_energy( const double &J ) const
    {
      SYS_T::commPrint("Warning: IMaterialModel_vol::get_Helmholtz_energy(J) is not implemented. \n");
      return 0.0;
    }

    virtual bool is_Helmholtz_supported() const = 0;
    
    // ------------------------------------------------------------------------
    // Input: Pressure p
    // Output: density rho = rho_0 (dG / dp)^(-1)
    // ------------------------------------------------------------------------
    virtual double get_rho( const double &p ) const
    {
      SYS_T::commPrint("Warning: IMaterialModel_vol::get_rho(p) is not implemented. \n");
      return 0.0;
    }

    // ------------------------------------------------------------------------
    // Output: The referencial density rho_0
    // ------------------------------------------------------------------------
    virtual double get_rho_0() const = 0;

    // ------------------------------------------------------------------------
    // Input: Pressure p
    // Output: drho / dp
    // ------------------------------------------------------------------------
    virtual double get_drho_dp( const double &p ) const
    {
      SYS_T::commPrint("Warning: IMaterialModel_vol::get_drho_dp(p) is not implemented. \n");
      return 0.0;
    }

    // ------------------------------------------------------------------------
    // Input: Pressure p
    // Output: isothermal compressibility factor beta := (drho / dp) / rho
    // ------------------------------------------------------------------------
    virtual double get_beta( const double &p ) const
    {
      SYS_T::commPrint("Warning: IMaterialModel_vol::get_beta(p) is not implemented. \n");
      return 0.0;
    }

    // ------------------------------------------------------------------------
    // Input: Pressure p
    // Output: dbeta / dp
    // ------------------------------------------------------------------------
    virtual double get_dbeta_dp( const double &p ) const
    {
      SYS_T::commPrint("Warning: IMaterialModel_vol::get_dbeta_dp(p) is not implemented. \n");
      return 0.0;
    }

    // ------------------------------------------------------------------------
    // Input: Volume ratio J
    // Output: the volumetric force sigma_vol := dW/dJ
    // Note: sigma_vol = -p
    // ------------------------------------------------------------------------
    virtual double get_vol_stress( const double &J ) const
    {
      SYS_T::commPrint("Warning: IMaterialModel_vol::get_vol_stress(J) is not implemented. \n");
      return 0.0;
    }

    // ------------------------------------------------------------------------
    // Input: Volume ratio J
    // Output: d sigma_vol / dp
    // Note: In Holzaplfe book pp. 254-255, the volumetric part of the
    // elasticity tensor is defined as J tilde_p invC otimes invC - 2Jp invC O
    // invC, in which tilde_p = sigma_vol + J dsigma_vol / dJ
    // ------------------------------------------------------------------------
    virtual double get_dvol_stress_dJ( const double &J ) const
    {
      SYS_T::commPrint("Warning: IMaterialModel_vol::get_dvol_stress_dJ(J) is not implemented. \n");
      return 0.0;
    }
};

#endif
