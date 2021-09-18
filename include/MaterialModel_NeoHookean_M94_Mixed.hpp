#ifndef MATERIALMODEL_NEOHOOKEAN_M94_MIXED_HPP
#define MATERIALMODEL_NEOHOOKEAN_M94_MIXED_HPP
// ==================================================================
// MaterialModel_NeoHookean_M94_Mixed.hpp
// Compressible Neo-Hookean model with the volumetric penalization
// form given by Miehe 94. This is a mixed formulation, where
// the pressure p is treated as an independent variable.
// 
// Density : rho = rho_0 ( 1 + p/k) 
// Isothermal Compressibility : beta = 1 / ( p + k )
// Isothermal Compressibility derivative : 
//                       d beta / dp = -1/ (p + k)^2
//
// Isochoric strain energy Psi_iso = (mu/2) (tr bar(C) - 3)
//                                 = (mu/2) (bar(I_1) - 3).
// See, Holzapfel book, p. 247.
//
// Date: Aug. 28 2017
// Author: Ju Liu
// Contact: liujuy@gmail.com
// ==================================================================
#include <cmath>
#include "IMaterialModel.hpp"

class MaterialModel_NeoHookean_M94_Mixed : public IMaterialModel
{
  public:
    MaterialModel_NeoHookean_M94_Mixed( const double &in_rho0,
        const double &in_E, const double &in_nu );

    MaterialModel_NeoHookean_M94_Mixed(
        const char * const &fname = "material_model.h5" );

    virtual ~MaterialModel_NeoHookean_M94_Mixed();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "NeoHookean-M94-Mixed";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5" ) const;

    // Returns the S_iso and P_iso
    virtual void get_PK( const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S ) const;

    // Returns the S_iso, P_iso, and CC_iso
    virtual void get_PK_Stiffness( const Matrix_3x3 &F, Matrix_3x3 &P,
        Matrix_3x3 &S, Tensor4_3D &CC ) const;

    // Return the Psi_iso
    virtual double get_strain_energy( const Matrix_3x3 &F ) const;

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
    double rho0, E, nu, lambda, mu, kappa;
    const double pt33, mpt67;
    const Matrix_3x3 I;
};

#endif
