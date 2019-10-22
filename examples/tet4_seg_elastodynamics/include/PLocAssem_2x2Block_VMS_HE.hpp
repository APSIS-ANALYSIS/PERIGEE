#ifndef PLOCASSEM_2X2BLOCK_VMS_HE_HPP
#define PLOCASSEM_2X2BLOCK_VMS_HE_HPP
// ==================================================================
// PLocAssem_2x2Block_VMS_HE.hpp
//
// This is the local assembly routine for finite elasticity in 3D
// using the block solver approach. Consistent VMS is added to
// provide pressure stabilization.
//
// Author: Ju Liu
// Date: Feb. 20 2018
// ==================================================================
#include "IPLocAssem_2x2Block.hpp"
#include "IMaterialModel.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_2x2Block_VMS_HE : public IPLocAssem_2x2Block
{
  public:
    PLocAssem_2x2Block_VMS_HE( IMaterialModel * const &in_matmodel,
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas );

    virtual ~PLocAssem_2x2Block_VMS_HE();

    virtual int get_dof() const {return 7;}

    virtual int get_dof_mat() const {return 4;}

    virtual int get_dof_mat_0() const {return 1;}

    virtual int get_dof_mat_1() const {return 3;}

    virtual int get_num_ebc_fun() const {return num_ebc_fun;}

    virtual void Zero_Tangent_Residual();

    virtual void Zero_Residual();

    virtual void Assem_Estimate();

    virtual void Assem_Mass_Residual(
        const double * const &disp,
        const double * const &pres,
        const double * const &velo,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &dot_disp,
        const double * const &dot_pres,
        const double * const &dot_velo,
        const double * const &disp,
        const double * const &pres,
        const double * const &velo,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_disp,
        const double * const &dot_pres,
        const double * const &dot_velo,
        const double * const &disp,
        const double * const &pres,
        const double * const &velo,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double &in_x, const double &in_y, const double &in_z, 
        const double * const &dot_disp,
        const double * const &dot_pres,
        const double * const &dot_velo,
        const double * const &disp,
        const double * const &pres,
        const double * const &velo,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

  private:
    const double rho0, alpha_f, alpha_m, gamma;

    double tau_m, tau_c;

    // memory layout
    const int nLocBas, dof_per_node, vec_size_0, vec_size_1;
    const int nqp, snLocBas, num_ebc_fun;

    // Useful tensor objects for material model
    IMaterialModel * matmodel;
    double rho, drho, mbeta, dmbeta, detF;
    Matrix_3x3 F, invF, P_iso, S_iso;

    Tensor4_3D AA_iso;

    // basis function containers
    double * R;
    double * dR_dx, * dR_dy, * dR_dz;
    double * d2R_dxx, * d2R_dxy, * d2R_dxz;
    double * d2R_dyy, * d2R_dyz, * d2R_dzz;
    double * dxi_dx;

    double ** Sub_Tan;

    void print_info() const;

    void get_tau( double &tau_m_qua, double &tau_c_qua,
        const double &dt, const double &Jin, const double &dx ) const;

    void get_f(const double &x, const double &y, const double &z,
        const double &t, double &fx, double &fy, double &fz ) const
    {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
    }

    void get_top_H(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
      if( (x>=0.0) && (x<=0.5) && (y>=0.0) && (y<=0.5) ) gz = -3.2e8 * t;
      else gz = 0.0;
    }

    void get_bot_H(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }

    void get_lef_H(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;

    }

    void get_rig_H(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;

      gx *= -1.0; gy *= -1.0; gz *= -1.0;
    }

    void get_fro_H(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }

    void get_bac_H(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;

      gx *= -1.0; gy *= -1.0; gz *= -1.0;
    }

    typedef void ( PLocAssem_2x2Block_VMS_HE::*locassem_2x2block_vms_he_funs )( const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const;

    locassem_2x2block_vms_he_funs * flist;

    void get_ebc_fun( const int &ebc_id,
        const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      return ((*this).*(flist[ebc_id]))(x,y,z,t,nx,ny,nz,gx,gy,gz);
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
