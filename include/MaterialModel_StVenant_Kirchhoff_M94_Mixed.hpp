#ifndef MATERIALMODEL_STVENANT_KIRCHHOFF_M94_MIXED_HPP
#define MATERIALMODEL_STVENANT_KIRCHHOFF_M94_MIXED_HPP
// ==================================================================
// MaterialModel_StVenant_Kirchhoff_M94_Mixed.hpp
// Compressible StVenant-Kirchhoff model with volumetric part given
// by the M94 formula. This is a mixed formulation, where 
// the pressure p is treated as an independent variable.
//
// Date: Aug. 14 2017
// Author: Ju Liu
// ==================================================================
#include "IMaterialModel.hpp"
#include "Math_Tools.hpp"

class MaterialModel_StVenant_Kirchhoff_M94_Mixed : public IMaterialModel
{
  public:
    MaterialModel_StVenant_Kirchhoff_M94_Mixed( const double &in_E,
        const double &in_nu );

    MaterialModel_StVenant_Kirchhoff_M94_Mixed( const double &in_rho0, 
        const double &in_E, const double &in_nu );

    MaterialModel_StVenant_Kirchhoff_M94_Mixed(
        const char * const &fname = "material_model.h5" );

    virtual ~MaterialModel_StVenant_Kirchhoff_M94_Mixed();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "StVenant-Kirchhoff-M94-Mixed";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5" ) const;

    virtual void get_PK( const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S ) const;

    virtual void get_PK_Stiffness( const Matrix_3x3 &F, Matrix_3x3 &P,
                Matrix_3x3 &S, Tensor4_3D &CC ) const;

    virtual double get_strain_energy( const Matrix_3x3 &F ) const;

    virtual double get_elastic_rho0() const {return rho0;}

    virtual double get_elastic_E() const {return E;}

    virtual double get_elastic_nu() const {return nu;}

    virtual double get_elastic_lambda() const {return lambda;}

    virtual double get_elastic_mu() const {return mu;}

    virtual double get_elastic_kappa() const {return kappa;}

    // Dialatational properties
    virtual double get_rho( const double &p ) const;

    virtual double get_drho_dp( const double &p ) const;

    virtual double get_beta( const double &p ) const;

    virtual double get_dbeta_dp( const double &p ) const;

  private:
    double rho0, E, nu, lambda, mu, kappa;
    const double pt33, mpt67;
    const Matrix_3x3 I;
};

#endif
