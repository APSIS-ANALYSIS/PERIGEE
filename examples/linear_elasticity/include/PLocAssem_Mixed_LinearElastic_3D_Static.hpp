#ifndef PLOCASSEM_MIXED_LINEARELASTIC_3D_STATIC_HPP
#define PLOCASSEM_MIXED_LINEARELASTIC_3D_STATIC_HPP
// ==================================================================
// PLocAssem_Mixed_LinearElastic_3D_Static.hpp
//
// This is the local assembly routine for the linear elastostatic
// problem in 3D displacement-pressure mixed formulation.
//
// Reference:
// T.J.R. Hughes, The Finite Element Method: Linear Static and Dynamic
// Finite Element Analysis.
//
// Author: Ju Liu
// Date: Dec. 22 2017
// ==================================================================
#include "IPLocAssem.hpp"

class PLocAssem_Mixed_LinearElastic_3D_Static : public IPLocAssem
{
  public:
    PLocAssem_Mixed_LinearElastic_3D_Static( const double &in_mat_E,
        const double &in_mat_nu, const int &in_nlocbas,
        const int &in_nqp, const int &in_snlocbas );

    virtual ~PLocAssem_Mixed_LinearElastic_3D_Static();

    virtual int get_dof() const {return dof_per_node;}

    virtual int get_num_ebc_fun() const {return num_ebc_fun;}

    virtual void Zero_Tangent_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
      for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 0.0;
    }

    virtual void Zero_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
    }

    virtual void Assem_Estimate()
    {
      for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 1.0;
    }

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

    virtual double get_model_para_1() const {return lambda;}
    virtual double get_model_para_2() const {return mu;}

  private:
    const int num_ebc_fun;

    const double E; // Young's modulus
    const double nu; // Possion's ratio
    const double lambda;
    const double mu;   // shear modulus mu = E / 2(1+nu)
    const double kappa; // bulk modulus kappa = lambda + 2mu / 3

    const int nLocBas, dof_per_node, vec_size, nqp, snLocBas;

    double * R;
    double * dR_dx;
    double * dR_dy;
    double * dR_dz;

    double ** Sub_Tan;
    
    void print_info() const;

    void get_f(const double &x, const double &y, const double &z,
        double &fx, double &fy, double &fz ) const
    {
      fx = 0.0; fy = 0.0; fz = 0.0;
    }
    
    void get_lef_H(const double &x, const double &y, const double &z,
        const double &nx, const double &ny, const double &nz,
        double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }

    void get_rig_H(const double &x, const double &y, const double &z,
        const double &nx, const double &ny, const double &nz,
        double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }

    void get_top_H(const double &x, const double &y, const double &z,
        const double &nx, const double &ny, const double &nz,
        double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }

    void get_bot_H(const double &x, const double &y, const double &z,
        const double &nx, const double &ny, const double &nz,
        double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }
    
    void get_fro_H(const double &x, const double &y, const double &z,
        const double &nx, const double &ny, const double &nz,
        double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }
    
    void get_bac_H(const double &x, const double &y, const double &z,
        const double &nx, const double &ny, const double &nz,
        double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }

    typedef void ( PLocAssem_Mixed_LinearElastic_3D_Static::*locassem_mle3d_funs )( const double &x, const double &y, const double &z,
        const double &nx, const double &ny, const double &nz,
        double &gx, double &gy, double &gz ) const;

    locassem_mle3d_funs * flist;

    void get_ebc_fun( const int &ebc_id,
        const double &x, const double &y, const double &z,
        const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      return ((*this).*(flist[ebc_id]))(x,y,z,nx,ny,nz,gx,gy,gz);
    }

    void Zero_Sub_Tan()
    {
      for(int ii=0; ii<16; ++ii)
      {
        for(int jj=0; jj<nLocBas * nLocBas; ++jj) Sub_Tan[ii][jj] = 0.0;
      }
    }
};

#endif
