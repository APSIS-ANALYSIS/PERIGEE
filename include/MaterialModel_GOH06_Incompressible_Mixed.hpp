#ifndef MATERIALMODEL_GOH06_INCOMPRESSIBLE_MIXED_HPP
#define MATERIALMODEL_GOH06_INCOMPRESSIBLE_MIXED_HPP
// ==================================================================
// MaterialModel_GOH06_Incompressible_Mixed.hpp
// Incompressible Gasser-Odgen-Holzapfel model. This material model
// is prepared for a mixed formulation. Since it is incompressible,
// rho = rho_0, beta = 0, and dbeta/dp = 0.
//
// The isochoric part contains the groundmatrix, a (soft) Neo-Hookean
// material and two families of fibres in a prescribed direction, using
// Fung-type model.
//
// The fibre direction is given be two angles: f_theta and f_phi:
// ( sin(f_theta) cos(f_phi), sin(f_theta) sin(f_phi), cos(f_theta) ).
// See Figure 2 in the following reference paper. Their units are degrees.
// We will transform them into radian and the vectors.
// 
// E.G. In the adventitial strip tensile test, the fibres are in the 
// x-z plane, and the angles are f_phi = 0.0 degree.
//
// There are two families of fibres in this model.
//
// Reference: Hyperelastic modeling of arterial layers with distributed
// collagen fibre orientations.
//
// Date: March 17 2017
// Author: Ju Liu, liujuy@gmail.com
// ==================================================================
#include "IMaterialModel.hpp"
#include "Math_Tools.hpp"
#include "HDF5_Reader.hpp"

class MaterialModel_GOH06_Incompressible_Mixed : public IMaterialModel
{
  public:
    MaterialModel_GOH06_Incompressible_Mixed( 
       const double &in_rho, const double &in_mu,
       const double &in_f1the, const double &in_f1phi,
       const double &in_f2the, const double &in_f2phi,
       const double &in_fk1, const double &in_fk2,
       const double &in_fkd );

    MaterialModel_GOH06_Incompressible_Mixed(
        const char * const &fname = "material_model.h5" );
    
    virtual ~MaterialModel_GOH06_Incompressible_Mixed();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "GOH06-Incompressible-Mixed";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5") const;

    virtual void get_PK( const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S);

    virtual void get_PK_Stiffness( const Matrix_3x3 &F, Matrix_3x3 &P,
                Matrix_3x3 &S, Tensor4_3D &CC );

    virtual double get_strain_energy(const Matrix_3x3 &F );

    // Read access to parameters
    virtual double get_elastic_rho0() const {return rho0;}

    virtual double get_elastic_E() const {return E;}

    virtual double get_elastic_nu() const {return nu;}

    virtual double get_elastic_mu() const {return mu;}

    // Dialatational part
    virtual double get_rho( const double &p ) const {return rho0;}

    virtual double get_drho_dp( const double &p ) const {return 0.0;}

    virtual double get_beta(const double &p) const {return 0.0;}

    virtual double get_dbeta_dp(const double &p) const {return 0.0;}

    virtual void get_fibre_dir( const int &dir,
        double &fa1, double &fa2, double &fa3 ) const;

  private:
    // useful constants
    const double pt33, mpt67, pi;

    // density
    double rho0;

    // groud matrix parameters 
    double E, nu, mu; // kappa = infty in this model

    // fibre direction: f_theta, f_phi in radians
    double f1_the, f1_phi, f2_the, f2_phi; 

    // fibre elastic parameter
    // fkd is the dispresion parameter ranging from 0 to 1/3
    double fk1, fk2, fkd;

    // unit vector for fibre direction in the ref domain    
    double a1[3], a2[3];

    const Matrix_3x3 I;
    Matrix_3x3 C, Cinv;

    // outproduct matrix for anisotropy 
    Matrix_3x3 a1xa1, a2xa2;

    // PxI  = I - 1/3 trC invC
    // PxH1 = PxI + (1-3kd)( aixa1 - 0.33 (ai dot C ai) invC )
    Matrix_3x3 PxI, PxH1, PxH2;

    double trC, detF, detFm0d67;

    double fE1, fE2;
};

#endif
