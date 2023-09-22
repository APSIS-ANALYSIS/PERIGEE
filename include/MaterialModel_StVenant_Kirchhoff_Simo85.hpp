#ifndef MATERIALMODEL_STVENANT_KIRCHHOFF_SIMO85_HPP
#define MATERIALMODEL_STVENANT_KIRCHHOFF_SIMO85_HPP
// ==================================================================
// MaterialModel_StVenant_Kirchhoff.hpp
// This is the modified Saint-Venant Kirchhoff model with a Simo85
// stabilizing term for the compression.
//
// Reference: G.A. Holzapfel Nonlinear solid mechanics, pp. 251
//
// Date: Dec. 15 2016
// Author: Ju Liu
// ==================================================================
#include <cmath>
#include "IMaterialModel.hpp"

class MaterialModel_StVenant_Kirchhoff_Simo85 : public IMaterialModel
{
  public:
    MaterialModel_StVenant_Kirchhoff_Simo85( const double &in_E, 
        const double &in_nu);

    MaterialModel_StVenant_Kirchhoff_Simo85(
        const char * const &fname = "material_model.h5" );

    ~MaterialModel_StVenant_Kirchhoff_Simo85();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "StVenant-Kirchhoff-Simo85";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5") const;

    virtual void get_PK( const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S ) const;

    virtual void get_PK_Stiffness( const Tensor2_3D &F, Tensor2_3D &P,
        Tensor2_3D &S, Tensor4_3D &CC) const;

    virtual double get_strain_energy( const Tensor2_3D &F ) const;

    virtual double get_elastic_E() const {return E;}

    virtual double get_elastic_nu() const {return nu;}

    virtual double get_elastic_lambda() const {return lambda;}

    virtual double get_elastic_mu() const {return mu;}

    virtual double get_elastic_kappa() const {return kappa;}

  private:
    double E, nu, lambda, mu, kappa;
    const Tensor2_3D I;
};

#endif
