#ifndef MATERIALMODEL_STVENANT_KIRCHHOFF_HPP
#define MATERIALMODEL_STVENANT_KIRCHHOFF_HPP
// ==================================================================
// MaterialModel_StVenant_Kirchhoff.hpp
// Saint-Venant Kirchhoff model.
//
// Reference: G.A. Holzapfel Nonlinear Solid Mechanics, pp 250.
//
// Date: Sept 15 2016
// Author: Ju Liu
// Contact: liujuy@gmail.com
// ==================================================================
#include <cmath>
#include "IMaterialModel.hpp"

class MaterialModel_StVenant_Kirchhoff : public IMaterialModel
{
  public:
    MaterialModel_StVenant_Kirchhoff( const double &in_E, const double &in_nu );

    MaterialModel_StVenant_Kirchhoff( const char * const &fname = "material_model.h5" );

    ~MaterialModel_StVenant_Kirchhoff();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "StVenant-Kirchhoff";
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
    double E, nu, lambda, mu, kappa;
    const Matrix_3x3 I;
};

#endif
