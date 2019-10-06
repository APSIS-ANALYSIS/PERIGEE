#ifndef MATERIALMODEL_NEOHOOKEAN_ST91_HPP
#define MATERIALMODEL_NEOHOOKEAN_ST91_HPP
// ==================================================================
// MaterialModel_NeoHookean_ST91.hpp
// Generalized compressible Neo-Hookean model with the volumetric 
// penalization form given by Simo-Taylor91.
//
// Reference: Y. Bazilevs, et al. Comput. Mech. 2008 43:3-37
//            J. Simo & R.L. Taylor, CMAME 1991 85:273-310 
// Date: Sept. 15 2016
// Author: Ju Liu
// Contact: liujuy@gmail.com
// ==================================================================
#include <cmath>
#include "IMaterialModel.hpp"

class MaterialModel_NeoHookean_ST91 : public IMaterialModel
{
  public:
    MaterialModel_NeoHookean_ST91( const double &in_rho,
        const double &in_E, const double &in_nu );

    virtual ~MaterialModel_NeoHookean_ST91();

    virtual void print_info() const;

    virtual void get_PK( const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S ); 

    virtual void get_PK_Stiffness( const Matrix_3x3 &F, Matrix_3x3 &P,
                Matrix_3x3 &S, Tensor4_3D &CC);

    virtual double get_strain_energy( const Matrix_3x3 &F );

    virtual double get_elastic_E() const {return E;}

    virtual double get_elastic_nu() const {return nu;}
    
    virtual double get_elastic_lambda() const {return lambda;}
    
    virtual double get_elastic_mu() const {return mu;}
    
    virtual double get_elastic_kappa() const {return kappa;}

    virtual double get_elastic_rho0() const {return rho0;}

  private:
    const double rho0, E, nu, lambda, mu, kappa;
    const double pt33, pt67, mpt67;
    double trC, detF, detF2, detFm0d67;
    
    const Matrix_3x3 I;
    Matrix_3x3 C, Cinv;
};

#endif
