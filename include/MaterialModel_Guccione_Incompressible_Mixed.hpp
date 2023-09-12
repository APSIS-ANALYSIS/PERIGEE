#ifndef MATERIALMODEL_GUCCIONE_INCOMPRESSIBLE_MIXED_HPP
#define MATERIALMODEL_GUCCIONE_INCOMPRESSIBLE_MIXED_HPP
// ==================================================================
// MaterialModel_Guccione_Incompressible_Mixed.hpp
// Incompressible Guccione model for Cardiac mechanics.
// 
// The model is an anisotropic model with the isochoric part of the
// free energy written as
//
// W( E' ) = 0.5 C [exp(Q) - 1].
// Q = bf E'ff^2 + bt[E'ss^2 + E'nn^2 + E'sn^2 + E'ns^2]
//     + bft [E'fs^2 + E'sf^2 + E'fn^2 + E'nf^2].
// Here E := 0.5 (C-I) , E' = 0.5 (C' - I), C' = J^-2/3 C,
//      E'_{fsn} = R^T E' R,
//      R = [ f^T;
//            s^T;
//            n^T ].
//      f is the unit vector in the fibre direction
//      s is the unit vector normal to the collagen sheets
//      n = f x s / ||fxs||.
//      R is the matrix that transform the local coordinates to the 
//      cartisian coordinates.  
//
// Reference: E. Garcia-Blanco, et al. CMAME 348 (2019) 796-845.
// Author: Ju Liu
// Date: March 12 2019
// ==================================================================
#include "IMaterialModel.hpp"
#include "Math_Tools.hpp"

class MaterialModel_Guccione_Incompressible_Mixed : public IMaterialModel
{
  public:
    MaterialModel_Guccione_Incompressible_Mixed(
        const double &in_rho, const double &in_C,
        const double &in_bf, const double &in_bt, const double &in_bft,
        const double &fx, const double &fy, const double &fz,
        const double &sx, const double &sy, const double &sz );

    MaterialModel_Guccione_Incompressible_Mixed( 
        const char * const &fname = "material_model.h5" );

    virtual ~MaterialModel_Guccione_Incompressible_Mixed();

    virtual void print_info() const;

    virtual std::string get_model_name() const
    {
      const std::string mname = "Guccione-Incompressible-Mixed";
      return mname;
    }

    virtual void write_hdf5( const char * const &fname = "material_model.h5") const;

    virtual void get_PK( const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S ) const;

    virtual void get_PK_Stiffness( const Matrix_3x3 &F, Matrix_3x3 &P, 
        Matrix_3x3 &S, Tensor4_3D &CC ) const;

    virtual double get_strain_energy( const Matrix_3x3 &F ) const;

    // Read access to material parameters
    virtual double get_elastic_rho0() const {return rho0;}

    virtual double get_elastic_nu() const {return 0.5;}

    // Dialatational part
    virtual double get_rho( const double &p ) const {return rho0;}

    virtual double get_drho_dp( const double &p ) const {return 0.0;}

    virtual double get_beta(const double &p) const {return 0.0;}

    virtual double get_dbeta_dp(const double &p) const {return 0.0;}

  private:
    // useful constants
    const double pt33, mpt67;

    // Material density
    double rho0;

    // Parameters appearing in front of expQ, and inside Q definition.
    double Cq, b_f, b_t, b_ft;

    // unit vector for fibre and firbre sheet orientation
    Vector_3 f, s, n;

    const Matrix_3x3 I;

    // The rotation matrix
    Matrix_3x3 R;
};

#endif
