#ifndef MATERIALMODEL_GOH14_ST91_Mixed_HPP
#define MATERIALMODEL_GOH14_ST91_Mixed_HPP
// ==================================================================
// MaterialModel_GOH14_ST91_Mixed.hpp
// Quasi-incompressible Gasser-Odgen-Holzapfel model. The model's
// volumetric part is based on the Simo-Taylor91 penalization. The
// isochoric isotropic part is the same as the GOH06_Incompressible.
// The anisotropic part is modified by replacing the modified right
// Cauchy-Green tensor by the right Cauchy-Green tensor. This idea 
// is from Nolan et al. J Mech Behav Biomed Mater 39 (2014) 48-60
//
// Auther: Xinhai Yue
// =================================================================
#include "IMaterialModel.hpp"
#include "Math_Tools.hpp"

class MaterialModel_GOH14_ST91_Mixed : public IMaterialModel
{
  public:
    MaterialModel_GOH14_ST91_Mixed(
        const double &in_rho, const double &in_E, const double &in_nu,
        const double &in_f1the, const double &in_f1phi,
        const double &in_f2the, const double &in_f2phi,
        const double &in_fk1, const double &in_fk2,
        const double &in_fkd );
    
    MaterialModel_GOH14_ST91_Mixed(
        const char * const &fname = "material_model.h5" );

    virtual ~MaterialModel_GOH14_ST91_Mixed();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "GOH14-ST91-Mixed";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5") const;

    virtual void get_PK( const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S) const;

    virtual void get_PK_Stiffness( const Tensor2_3D &F, Tensor2_3D &P,
        Tensor2_3D &S, Tensor4_3D &CC ) const;

    virtual double get_strain_energy(const Tensor2_3D &F ) const;

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

  private:
    double rho0, E, nu, lambda, mu, kappa;
    double f1_the, f1_phi, f2_the, f2_phi;
    double fk1, fk2, fkd;
    Vector_3 a1, a2;
};

#endif
