#ifndef MATERIALMODEL_GOH06_ST91_MIXED_HPP
#define MATERIALMODEL_GOH06_ST91_MIXED_HPP
// ==================================================================
// MaterialModel_GOH06_ST91_Mixed.hpp
// Quasi-incompressible Gasser-Odgen-Holzapfel model. This material
// model's isochoric part is the same as the GOH06_Incompressible 
// routine. The volumetric part is based on the Simo-Taylor91 
// penalization.
//
// Date: March 24 2017
// Author: Ju Liu, liujuy@gmail.com
// ==================================================================
#include "IMaterialModel.hpp"
#include "Math_Tools.hpp"

class MaterialModel_GOH06_ST91_Mixed : public IMaterialModel
{
  public:
    MaterialModel_GOH06_ST91_Mixed( 
        const double &in_rho, const double &in_E, const double &in_nu,
        const double &in_f1the, const double &in_f1phi,
        const double &in_f2the, const double &in_f2phi,
        const double &in_fk1, const double &in_fk2,
        const double &in_fkd );

    MaterialModel_GOH06_ST91_Mixed(
        const char * const &fname = "material_model.h5" );

    virtual ~MaterialModel_GOH06_ST91_Mixed();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "GOH06-ST91-Mixed";
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

    virtual double get_elastic_lambda() const {return lambda;}

    virtual double get_elastic_mu() const {return mu;}

    virtual double get_elastic_kappa() const {return kappa;}

    virtual double get_rho( const double &p ) const;

    virtual double get_drho_dp( const double &p ) const;

    virtual double get_beta( const double &p ) const;

    virtual double get_dbeta_dp( const double &p ) const;

    virtual Vector_3 get_fibre_dir( const int &dir ) const;

    // Update fibre direction by input direction basis for each point.
    // Notice that the original fibre direction a1 and a2 are difined on the basis vector
    // e1 = [1, 0, 0], e2 = [0, 1, 0], and e3 = [0, 0, 1].
    // And the input vector basis_r, basis_c, and basis_l correspond to e1, e2, and e3,
    // respectively.
    virtual void update_fibre_dir( const Vector_3 &basis_r, const Vector_3 &basis_l,
        const Vector_3 &basis_c );

  private:
    // useful constants
    const double pt33, mpt67, pi;

    // density
    double rho0;

    // groud matrix parameters 
    double E, nu, lambda, mu, kappa;

    // fibre direction: f_theta, f_phi in radians
    double f1_the, f1_phi, f2_the, f2_phi;

    // fibre elastic parameter
    // fkd is the dispresion parameter ranging from 0 to 1/3
    double fk1, fk2, fkd;

    // unit vector for fibre direction in the ref domain    
    Vector_3 a1, a2;

    // unit vector for fibre direction difined on the basis vector 
    // e1 = [1, 0, 0], e2 = [0, 1, 0], and e3 = [0, 0, 1]
    Vector_3 a01, a02;
};

#endif
