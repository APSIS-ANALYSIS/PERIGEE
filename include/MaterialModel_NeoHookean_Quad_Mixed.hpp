#ifndef MATERIALMODEL_NEOHOOKEAN_QUAD_MIXED_HPP
#define MATERIALMODEL_NEOHOOKEAN_QUAD_MIXED_HPP
// ==================================================================
// MaterialModel_NeoHookean_Quad_Mixed.hpp
// Compressible Neo-Hookean model with the volumetric penalization
// form given by the quadratic form. This is a mixed formulation, 
// where the pressure p is treated as an independent variable.
// 
// Dialatational part:
// Density : rho = rho_0  /  ( 1 - p/k )
// Isothermal Compressibility : beta = 1 / ( k - p )
//
// Isochoric strain energy Psi_iso = (mu/2) (tr bar(C) - 3)
//                                 = (mu/2) (bar(I_1) - 3).
//
// Date: Aug. 22 2016
// Author: Ju Liu
// Contact: liujuy@gmail.com
// ==================================================================
#include <cmath>
#include "IMaterialModel.hpp"

class MaterialModel_NeoHookean_Quad_Mixed : public IMaterialModel
{
  public:
    MaterialModel_NeoHookean_Quad_Mixed( const double &in_rho,
        const double &in_E, const double &in_nu );

    virtual ~MaterialModel_NeoHookean_Quad_Mixed();

    virtual void print_info() const;

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

    // A empty fibre_dir function
    virtual void get_fibre_dir( const int &dir,
        double &fa1, double &fa2, double &fa3 ) const {}

  private:
    const double rho0;
    const double E, nu, lambda, mu, kappa;
    const double pt33, mpt67;
    double trC, detF, detFm0d67;

    const Matrix_3x3 I;
    Matrix_3x3 C, Cinv;
};

#endif
