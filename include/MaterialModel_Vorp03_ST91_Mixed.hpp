#ifndef MATERIALMODEL_VORP03_ST91_MIXED_HPP
#define MATERIALMODEL_VORP03_ST91_MIXED_HPP
// ==================================================================
// MaterialModel_Vorp03_ST91_Mixed.hpp
// Compressible Vorp03 model with the volumetric penalization
// form given by Simo-Taylor91. This is a mixed formulation, where
// the pressure p is treated as an independent variable.
// 
// Dialatational part:
// Density : rho = rho_0  /  ( p/k + sqrt(p^2/k^2 + 1) )
// Isothermal Compressibility : beta = -1 / sqrt( p^2 + k^2 )
// Isothermal Compressibility derivative : 
//                       d beta / dp = p (p^2 + k^2)^-1.5
//
// Date: Aug 25 2024
// Author: Xinhai Yue
// ==================================================================
#include <cmath>
#include "IMaterialModel.hpp"

class MaterialModel_Vorp03_ST91_Mixed : public IMaterialModel
{
  public:
    MaterialModel_Vorp03_ST91_Mixed( const double &in_rho0, 
        const double &in_E, const double &in_nu,
        const double &in_c1, const double &in_c2 );

    MaterialModel_Vorp03_ST91_Mixed(
        const char * const &fname = "material_model.h5" );


    virtual ~MaterialModel_Vorp03_ST91_Mixed();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "Vorp-ST91-Mixed";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5" ) const;

    // Returns the S_iso and P_iso
    virtual void get_PK( const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S ) const;

    // Returns the S_iso, P_iso, and CC_iso
    virtual void get_PK_Stiffness( const Tensor2_3D &F, Tensor2_3D &P,
        Tensor2_3D &S, Tensor4_3D &CC ) const;

    // Return the Psi_iso
    virtual double get_strain_energy( const Tensor2_3D &F ) const;

    // Elastic material properties
    virtual double get_elastic_rho0() const {return rho0;}

    virtual double get_elastic_E() const {return E;}

    virtual double get_elastic_nu() const {return nu;}

    virtual double get_elastic_lambda() const {return lambda;}

    virtual double get_elastic_mu() const {return mu;}

    virtual double get_elastic_kappa() const {return kappa;}

    // Dialatational properties
    // rho 
    virtual double get_rho( const double &p ) const;

    // drho = drho/dp
    virtual double get_drho_dp( const double &p ) const;

    // beta
    virtual double get_beta( const double &p ) const;

    // dbeta
    virtual double get_dbeta_dp( const double &p ) const;

  private:
    const double pt33, mpt67;
    double rho0, E, nu, lambda, mu, kappa, c1, c2;
    const Tensor2_3D I;
};

#endif
