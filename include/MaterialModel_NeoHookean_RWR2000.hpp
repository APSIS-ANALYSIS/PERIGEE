#ifndef MATERIALMODEL_NEOHOOKEAN_RWR2000_HPP
#define MATERIALMODEL_NEOHOOKEAN_RWR2000_HPP
// ==================================================================
// MaterialModel_NeoHookean_RWR2000.hpp
// Generalized compressible Neo-Hookean model with volumetric 
// penalization. This model is used in the 3D block compression bench
// -mark problem for (nearly) incompressible solids and numerical
// schemes.
//
// Reference: S. Resse, M. Kussner, B.D. Reddy, Comput. Struct. 75
//            2000: 291-304.
//            T. Elguedj, Y. Bazilevs, V.M. Calo, & T.J.R. Hughes,
//            CMAME 197 2008: pp 2752.
// ==================================================================
#include <cmath>
#include "IMaterialModel.hpp"

class MaterialModel_NeoHookean_RWR2000 : public IMaterialModel
{
  public:
    MaterialModel_NeoHookean_RWR2000( const double &in_E, const double &in_nu );

    MaterialModel_NeoHookean_RWR2000(
        const char * const &fname = "material_model.h5" );

    virtual ~MaterialModel_NeoHookean_RWR2000();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "NeoHookean-RWR2000";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5" ) const;
    
    virtual void get_PK( const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S ) const;

    virtual void get_PK_Stiffness( const Tensor2_3D &F, Tensor2_3D &P,
        Tensor2_3D &S, Tensor4_3D &CC ) const;

    virtual double get_strain_energy( const Tensor2_3D &F ) const;

    virtual double get_elastic_E() const {return E;}

    virtual double get_elastic_nu() const {return nu;}

    virtual double get_elastic_lambda() const {return lambda;}

    virtual double get_elastic_mu() const {return mu;}

    virtual double get_elastic_kappa() const {return kappa;}

  private:
    double E, nu, lambda, mu, kappa;
};

#endif
