#ifndef PLOCASSEM_LPS0_SEG_GENALPHA_HPP
#define PLOCASSEM_LPS0_SEG_GENALPHA_HPP
// ==================================================================
// PLocAssem_LPS0_Seg_GenAlpha.hpp
//
// This is the local assembly routine for the local projection
// stabilized mixed formulation for the p-v implicit solver using
// the segregated algorithm.
// 
// LPS0 means the projection is projected to zero-th order polynomial
//
// The time implementation is the Generalized-alpha method.
//
// Author: Ju Liu
// Date created: Dec. 11 2017
// ==================================================================
#include "IPLocAssem.hpp"
#include "IMaterialModel.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_LPS0_Seg_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_LPS0_Seg_GenAlpha(
        IMaterialModel * const &in_matmodel,
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas );

    virtual ~PLocAssem_LPS0_Seg_GenAlpha();

    virtual int get_dof() const {return dof_per_node;}

    virtual int get_dof_mat() const {return 4;}

    virtual int get_num_ebc_fun() const {return num_ebc_fun;}

    virtual void Zero_Tangent_Residual();

    virtual void Zero_Residual();

    virtual void Assem_Estimate();

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
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

  private:
    const double alpha_f, alpha_m, gamma;
    
    const int num_ebc_fun;
    
    // memory layout
    const int nLocBas, dof_per_node, vec_size, nqp, snLocBas;

    // useful tensors for the material model
    IMaterialModel * matmodel;
    
    double rho, mbeta, drho, dmbeta, detF;

    Matrix_3x3 F, invF, P_iso, S_iso;

    Tensor4_3D AA_iso; // stiffness in the current configuration

    // basis function allocations
    double * R, * dR_dx, * dR_dy, * dR_dz;

    double ** Sub_Tan;

    void print_info() const;

    double get_tau( const double &dt, const double &dx ) const
    {
      return 1.0 * dx * dx;
    } 

    void get_f( const double &x, const double &y, const double &z,
        const double &t, double &fx, double &fy, double &fz ) const
    {
      fx = 0.0; fy = 0.0; fz = 0.0;
    }

    void get_rig_h( const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gz = 0.0;
      if(t<1.0) gy = 6.25e5 * t;
      else gy = 6.25e5;
    }

    typedef void ( PLocAssem_LPS0_Seg_GenAlpha::*locassem_lps0_funs )( 
        const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const;

    locassem_lps0_funs * flist;

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
