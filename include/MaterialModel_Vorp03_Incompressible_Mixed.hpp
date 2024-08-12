#ifndef MATERIALMODEL_VORP03_INCOMPRESSIBLE_MIXED_HPP
#define MATERIALMODEL_VORP03_INCOMPRESSIBLE_MIXED_HPP
// ==================================================================
// MaterialModel_Vorp03_Incompressible_Mixed.hpp
// Incompressible Vorp03 model with pressure rate mixed
// formulation.
//
// Reference: Effect of variation in intraluminal thrombus 
// constitutive properties on abdominal aortic aneurysm wall stress. 
//
// Date: Aug 12 2024
// Author: Xinhai Yue
// ==================================================================
#include "IMaterialModel.hpp"
#include "Math_Tools.hpp"

class MaterialModel_Vorp03_Incompressible_Mixed : public IMaterialModel
{
  public:
    MaterialModel_Vorp03_Incompressible_Mixed( const double &in_rho,
        const double &in_c1, const double &in_c2 );
  
    ~MaterialModel_Vorp03_Incompressible_Mixed();

    MaterialModel_Vorp03_Incompressible_Mixed(
        const char * const &fname = "material_model.h5" );

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "Vorp03-Incompressible-Mixed";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5" ) const;

    virtual void get_PK( const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S ) const;

    virtual void get_PK_Stiffness( const Tensor2_3D &F, Tensor2_3D &P,
        Tensor2_3D &S, Tensor4_3D &CC ) const;

    virtual double get_strain_energy( const Tensor2_3D &F ) const;

    virtual double get_elastic_rho0() const {return rho0;}

    virtual double get_elastic_nu() const {return 0.5;}

    virtual double get_rho( const double &p ) const {return rho0;}

    virtual double get_drho_dp( const double &p ) const {return 0.0;}

    virtual double get_beta( const double &p ) const {return 0.0;}

    virtual double get_dbeta_dp( const double &p ) const {return 0.0;}

  private:
    const double pt33, mpt67;

    double rho0, c1, c2;
};

#endif
