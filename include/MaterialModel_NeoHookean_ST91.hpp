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

    MaterialModel_NeoHookean_ST91(
        const char * const &fname = "material_model.h5" );

    virtual ~MaterialModel_NeoHookean_ST91();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "NeoHookean-ST91";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5" ) const;

    virtual void get_PK( const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S ) const; 

    virtual void get_PK_Stiffness( const Tensor2_3D &F, Tensor2_3D &P,
                Tensor2_3D &S, Tensor4_3D &CC) const;

    virtual double get_strain_energy( const Tensor2_3D &F ) const;

    virtual double get_elastic_E() const {return E;}

    virtual double get_elastic_nu() const {return nu;}
    
    virtual double get_elastic_lambda() const {return lambda;}
    
    virtual double get_elastic_mu() const {return mu;}
    
    virtual double get_elastic_kappa() const {return kappa;}

    virtual double get_elastic_rho0() const {return rho0;}

  private:
    double rho0, E, nu, lambda, mu, kappa;
    const double pt33, mpt67;
};

#endif
