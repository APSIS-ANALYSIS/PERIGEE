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

    virtual ~MaterialModel_NeoHookean_RWR2000();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "NeoHookean-RWR2000";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5") const;
    
    virtual void get_PK( const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S );

    virtual void get_PK_Stiffness( const Matrix_3x3 &F, Matrix_3x3 &P,
        Matrix_3x3 &S, Tensor4_3D &CC);

    virtual double get_strain_energy( const Matrix_3x3 &F );

    virtual double get_elastic_E() const {return E;}

    virtual double get_elastic_nu() const {return nu;}

    virtual double get_elastic_lambda() const {return lambda;}

    virtual double get_elastic_mu() const {return mu;}

    virtual double get_elastic_kappa() const {return kappa;}

  private:
    const double E, nu, lambda, mu, kappa;
    Matrix_3x3 Cinv;
};

#endif
