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

    virtual void write_hdf5( const char * const &fname = "material_model.h5") const;

    virtual void get_PK( const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S);

    virtual void get_PK_Stiffness( const Matrix_3x3 &F, Matrix_3x3 &P,
        Matrix_3x3 &S, Tensor4_3D &CC );

    virtual double get_strain_energy(const Matrix_3x3 &F );

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

    virtual void get_fibre_dir( const int &dir,
        double &fa1, double &fa2, double &fa3 ) const;

  private:
    const double pt33, mpt67, pi;
    double rho0, E, nu, lambda, mu, kappa;
    double f1_the, f1_phi, f2_the, f2_phi;
    double fk1, fk2, fkd;
    Vector_3 a1, a2;
    const Matrix_3x3 I;
};

#endif
