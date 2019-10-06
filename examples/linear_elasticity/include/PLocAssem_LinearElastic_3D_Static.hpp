#ifndef PLOCASSEM_LINEARELASTIC_3D_STATIC_HPP
#define PLOCASSEM_LINEARELASTIC_3D_STATIC_HPP
// ==================================================================
// PLocAssem_LinearElastic_3D_Static.hpp
//
// This is the local assembly routien for the linear elastostatic
// problems in 3D.
//
// Reference:
// T.J.R. Hughes, The Finite Element Method: Linear Static and Dynamic
// Finite Element Analysis.
//
// Author: Ju Liu
// Date: May 9 2017
// ==================================================================
#include "IPLocAssem.hpp"

class PLocAssem_LinearElastic_3D_Static : public IPLocAssem
{
  public:
    PLocAssem_LinearElastic_3D_Static( const double &in_mat_E,
        const double &in_mat_nu, const int &in_nlocbas, 
        const int &in_nqp, const int &in_snlocbas );
    
    virtual ~PLocAssem_LinearElastic_3D_Static();

    virtual int get_dof() const {return dof_per_node;}

    virtual int get_num_ebc_fun() const {return num_ebc_fun;}

    virtual void Zero_Tangent_Residual();

    virtual void Zero_Residual();

    virtual void Assem_Estimate();

    // Assembly routines
    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad ){};

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Mass_Residual(
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad ){};

    // Elemental BC assembly
    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    // Functions to pass material parameter out
    virtual double get_model_para_1() const {return lambda;}
    virtual double get_model_para_2() const {return mu;}

  private:
    const int num_ebc_fun;

    // 1. Material properties
    const double E; // Young's modulus
    const double nu; // Possion's ratio
    // Lame coeff 1 for linear solidi lambda = nu E / (1+nu)(1-2nu) 
    const double lambda; 
    const double mu;   // shear modulus mu = E / 2(1+nu)
    const double kappa; // bulk modulus kappa = lambda + 2mu / 3

    // 2. Memory layout
    const int nLocBas, dof_per_node, vec_size;
    const int nqp, snLocBas;

    // 5. Basis functions allocations
    double * R;
    double * dR_dx;
    double * dR_dy;
    double * dR_dz;

    double ** Sub_Tan;

    // -------------------------------------------------------------- 
    // PRIVATE FUNCTIONS
    // -------------------------------------------------------------- 
    void print_info() const;

    // --------------------------------------------------------------
    // External Forces
    // --------------------------------------------------------------
    // External volumetric body force vector
    void get_f(const double &x, const double &y, const double &z,
        double &fx, double &fy, double &fz ) const
    {
      //const double pi = MATH_T::PI;
      //const double lambda = 0.5769231; 
      //const double mu = 0.3846154;

      fx = 0.0; 
      
      fy = 0.0;
      
      fz = 0.0;
    }


    // External boundary traction
    // get_lef_H : return left surface's traction 
    void get_lef_H(const double &x, const double &y, const double &z,
        const double &nx, const double &ny, const double &nz,
        double &gx, double &gy, double &gz ) const
    {
      //const double pi = MATH_T::PI;
      //const double lambda = 0.5769231; 
      //const double mu = 0.3846154;
      
      gx = 0.0; 
      
      gy = 0.0;

      gz = 0.0;
    }

    // get_rig_H : return right surface's traction 
    void get_rig_H(const double &x, const double &y, const double &z,
        const double &nx, const double &ny, const double &nz,
        double &gx, double &gy, double &gz ) const
    {
      //const double pi = MATH_T::PI;
      //const double lambda = 0.5769231; 
      //const double mu = 0.3846154;
      
      gx = 0.0; 
      
      gy = 0.0;

      gz = 0.0;
    }

    // get_top_H : return right surface's traction 
    void get_top_H(const double &x, const double &y, const double &z,
        const double &nx, const double &ny, const double &nz,
        double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }

    // get_bot_H : return right surface's traction 
    void get_bot_H(const double &x, const double &y, const double &z,
        const double &nx, const double &ny, const double &nz,
        double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }
    
    // get_fro_H : return front surface's traction 
    void get_fro_H(const double &x, const double &y, const double &z,
        const double &nx, const double &ny, const double &nz,
        double &gx, double &gy, double &gz ) const
    {
      //const double pi = MATH_T::PI;
      //const double lambda = 0.5769231; 
      //const double mu = 0.3846154;
      
      gx = 0.0; 
      
      gy = 0.0;

      gz = 0.0;
    }

    // get_bac_H : return back surface's traction 
    void get_bac_H(const double &x, const double &y, const double &z,
        const double &nx, const double &ny, const double &nz,
        double &gx, double &gy, double &gz ) const
    {
      //const double pi = MATH_T::PI;
      //const double lambda = 0.5769231; 
      //const double mu = 0.3846154;
      
      gx = 0.0; 
      
      gy = 0.0;
      
      gz = 0.0; 
    }

    typedef void ( PLocAssem_LinearElastic_3D_Static::*locassem_le3d_fem_funs )( const double &x, const double &y, const double &z,
        const double &nx, const double &ny, const double &nz, 
        double &gx, double &gy, double &gz ) const;

    locassem_le3d_fem_funs * flist;

    void get_ebc_fun( const int &ebc_id,
        const double &x, const double &y, const double &z,
        const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      return ((*this).*(flist[ebc_id]))(x,y,z,nx,ny,nz,gx,gy,gz);
    }

    void Zero_Sub_Tan()
    {
      for(int ii=0; ii<9; ++ii)
      {
        for(int jj=0; jj<nLocBas * nLocBas; ++jj) Sub_Tan[ii][jj] = 0.0;
      }
    }
};


#endif
