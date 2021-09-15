#ifndef MATERIALMODEL_NEOHOOKEAN_ST91_MIXED_HPP
#define MATERIALMODEL_NEOHOOKEAN_ST91_MIXED_HPP
// ==================================================================
// MaterialModel_NeoHookean_ST91_Mixed.hpp
// Compressible Neo-Hookean model with the volumetric penalization
// form given by Simo-Taylor91. This is a mixed formulation, where
// the pressure p is treated as an independent variable.
// 
// Dialatational part:
// Density : rho = rho_0  /  ( p/k + sqrt(p^2/k^2 + 1) )
// Isothermal Compressibility : beta = -1 / sqrt( p^2 + k^2 )
// Isothermal Compressibility derivative : 
//                       d beta / dp = p (p^2 + k^2)^-1.5
//
// Isochoric strain energy Psi_iso = (mu/2) (tr bar(C) - 3)
//                                 = (mu/2) (bar(I_1) - 3).
// See, Holzapfel book, p. 247.
//
// Date: Dec. 2 2016
// Author: Ju Liu
// Contact: liujuy@gmail.com
// ==================================================================
#include <cmath>
#include "IMaterialModel.hpp"

class MaterialModel_NeoHookean_ST91_Mixed : public IMaterialModel
{
  public:
    MaterialModel_NeoHookean_ST91_Mixed( const double &in_E, 
        const double &in_nu );

    MaterialModel_NeoHookean_ST91_Mixed( const double &in_rho0, 
        const double &in_E, const double &in_nu );

    MaterialModel_NeoHookean_ST91_Mixed(
        const char * const &fname = "material_model.h5" );


    virtual ~MaterialModel_NeoHookean_ST91_Mixed();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "NeoHookean-ST91-Mixed";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5") const;

    // Returns the S_iso and P_iso
    virtual void get_PK( const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S );

    // Returns the S_iso, P_iso, and CC_iso
    virtual void get_PK_Stiffness( const Matrix_3x3 &F, Matrix_3x3 &P,
        Matrix_3x3 &S, Tensor4_3D &CC );

    // Return the Psi_iso
    virtual double get_strain_energy( const Matrix_3x3 &F );

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
    double rho0;
    double E, nu, lambda, mu, kappa;
    const double pt33, mpt67;
    double trC, detF, detFm0d67;

    const Matrix_3x3 I;
    Matrix_3x3 C, Cinv;
};

#endif
