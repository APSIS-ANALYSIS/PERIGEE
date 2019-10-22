#ifndef PLOCASSEM_VMS_SEG_HYPERELASTIC_3D_FEM_GENALPHA_HPP
#define PLOCASSEM_VMS_SEG_HYPERELASTIC_3D_FEM_GENALPHA_HPP
// ==================================================================
// PLocAssem_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha.hpp
//
// This is the local assembly routine for the stabilized mixed 
// formulation for the p-v implicit solver using the segregated
// algorithm. This assembly routine does NOT assume linear elements.
// Second order derivatives are calculated in the assembly.
//
// Time integration is the Generalized-alpha method.
//
// Author: Ju Liu
// Date Created: May 15 2017
// ==================================================================
#include "IPLocAssem.hpp"
#include "IMaterialModel.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha(
        IMaterialModel * const &in_matmodel,
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas );

    virtual ~PLocAssem_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha();

    virtual int get_dof() const {return dof_per_node;}

    virtual int get_dof_mat() const {return 4;}

    virtual int get_num_ebc_fun() const {return num_ebc_fun;}

    virtual void Zero_Tangent_Residual();

    virtual void Zero_Residual();

    virtual void Assem_Estimate();

    // Assembly routines.
    // Note: vec_a / vec_b has 3 + 1 + 3 = 7 d.o.f.
    //       vec_a : dot(variable), i.e., velo,
    //       vec_b : variable, i.e., disp.    
    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    
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
        const IQuadPts * const &quad );


    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double &in_x, const double &in_y, const double &in_z,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );


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

  private:
    const double rho0, alpha_f, alpha_m, gamma;

    const int num_ebc_fun;

    double tau_m, tau_c;

    // memory layout
    // dof_per_node = 7 to make it compatible with the problem setting
    // vec_size = 4 * nLocBas, which defines the local matrix/vector length
    const int nLocBas, dof_per_node, vec_size;
    const int nqp;
    const int snLocBas;

    // useful tensors for the material model
    IMaterialModel * matmodel;
    double rho, drho, mbeta, dmbeta, detF;
    Matrix_3x3 F, invF, P_iso, S_iso;

    Tensor4_3D AA_iso;

    // basis function allocations
    double * R;
    double * dR_dx;
    double * dR_dy;
    double * dR_dz;
    double * d2R_dxx;
    double * d2R_dxy;
    double * d2R_dxz;
    double * d2R_dyy;
    double * d2R_dyz;
    double * d2R_dzz;

    double * dxi_dx;

    double ** Sub_Tan;

    void print_info() const;

    void get_tau( double &tau_m_qua, double &tau_c_qua,
        const double &dt, const double &Jin,
        const double &dx ) const;

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
      
      double val;
      if( t<= 1.0 ) val = -3.2e8 * t;
      else val = -3.2e8;

      if( (x>=0.0) && (x<=0.5) && (y>=0.0) && (y<=0.5) ) gz = val;
      else gz = 0.0;
      
      //gx = 6.66666667e4 * t; gy = 0.0; gz = 0.0;
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
    }


    // Use pointers to the member functions to facilitate the automatic
    // treatment of ebc surface integration.
    typedef void ( PLocAssem_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha::*locassem_vms_seg_ela_fem_funs )( const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const;

    locassem_vms_seg_ela_fem_funs * flist;

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
