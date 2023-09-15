#ifndef MATERIALMODEL_NEOHOOKEAN_INCOMPRESSIBLE_MIXED_HPP
#define MATERIALMODEL_NEOHOOKEAN_INCOMPRESSIBLE_MIXED_HPP
// ==================================================================
// MaterialModel_NeoHookean_Incompressible_Mixed.hpp
// Incompressible Neo-Hookean model with pressure rate mixed
// formulation.
//
// Date: March 22 2017
// Author: Ju Liu, liujuy@gmail.com
// ==================================================================
#include "IMaterialModel.hpp"
#include "Math_Tools.hpp"

class MaterialModel_NeoHookean_Incompressible_Mixed : public IMaterialModel
{
  public:
    MaterialModel_NeoHookean_Incompressible_Mixed( const double &in_mu );
    
    MaterialModel_NeoHookean_Incompressible_Mixed( 
        const double &in_rho, const double &in_E );
  
    ~MaterialModel_NeoHookean_Incompressible_Mixed();

    MaterialModel_NeoHookean_Incompressible_Mixed(
        const char * const &fname = "material_model.h5" );

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "NeoHookean-Incompressible-Mixed";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5" ) const;

    virtual void get_PK( const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S ) const;

    virtual void get_PK_Stiffness( const Tensor2_3D &F, Tensor2_3D &P,
        Tensor2_3D &S, Tensor4_3D &CC ) const;

    virtual double get_strain_energy( const Tensor2_3D &F ) const;

    virtual double get_elastic_rho0() const {return rho0;}

    virtual double get_elastic_E() const {return E;}

    virtual double get_elastic_nu() const {return nu;}

    virtual double get_elastic_mu() const {return mu;}

    virtual double get_rho( const double &p ) const {return rho0;}

    virtual double get_drho_dp( const double &p ) const {return 0.0;}

    virtual double get_beta( const double &p ) const {return 0.0;}

    virtual double get_dbeta_dp( const double &p ) const {return 0.0;}

  private:
    const double pt33;

    double rho0, E, nu, mu;

    const Tensor2_3D I;
};

#endif
