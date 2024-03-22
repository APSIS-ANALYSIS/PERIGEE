#ifndef MATERIALMODEL_LINEAR_ELASTICITY_HPP
#define MATERIALMODEL_LINEAR_ELASTICITY_HPP
// ==================================================================
// MaterialModel_Linear_Elasticity.hpp
// 
// Linear elasticity model. This material model is prepared for
// stress smoothing.
//
// Date: Feb 4 2024
// Author: Xinhai Yue
// ==================================================================
#include "IMaterialModel.hpp"
#include "Math_Tools.hpp"

class MaterialModel_Linear_Elasticity : public IMaterialModel
{
  public:
    MaterialModel_Linear_Elasticity(
       const double &in_modulus_E, const double &in_nu );

    MaterialModel_Linear_Elasticity(
        const char * const &fname = "material_model.h5" );
    
    virtual ~MaterialModel_Linear_Elasticity();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      return std::string("Linear_Elasticity");
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5") const;

    virtual void get_PK( const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S ) const
    {
      SYS_T::commPrint("Warning: MaterialModel::get_PK() is not implemented. \n");
    }

    virtual void get_PK_Stiffness( const Tensor2_3D &F, Tensor2_3D &P,
        Tensor2_3D &S, Tensor4_3D &CC ) const
    {
      SYS_T::commPrint("Warning: MaterialModel::get_PK_Stiffness() is not implemented. \n");
    }

    // Read access to parameters
    virtual double get_elastic_E() const {return modulus_E;}

    virtual double get_elastic_nu() const {return nu;}

    virtual double get_elastic_mu() const {return 0.5 * modulus_E / (1.0 + nu);}

    virtual double get_elastic_lambda() const {return nu * modulus_E / ((1.0 + nu) * (1.0 - 2.0 * nu));}

    virtual Tensor2_3D get_Cauchy_stress( const Tensor2_3D &F ) const;

  private:
    double modulus_E, nu;
};

#endif