#ifndef PLOCASSEM_TET4_VMS_NS_3D_GENALPHA_HPP
#define PLOCASSEM_TET4_VMS_NS_3D_GENALPHA_HPP
// ==================================================================
// PLocAssem_Tet4_VMS_NS_3D_GenAlpha.hpp
// 
// This is the local assembly routine for VMS formulation of the 3D
// Navier-Stokes equations with Generalized-alpha for time stepping.
//
// This code assumes the usage of four-node tetradedral element, which
// simplifies the implementation w.r.t. second-order derivatives.
//
// Ref: CMAME 197 (2007) 173-201
// Author: Ju Liu
// Date Created: June 16 2017
// ==================================================================
#include "IPLocAssem.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_Tet4_VMS_NS_3D_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_Tet4_VMS_NS_3D_GenAlpha(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas );

    virtual ~PLocAssem_Tet4_VMS_NS_3D_GenAlpha();

    virtual int get_dof() const {return dof_per_node;}

    virtual int get_dof_mat() const {return 4;}

    virtual int get_num_ebc_fun() const {return num_ebc_fun;}

    virtual void Zero_Tangent_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
      for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 0.0;
    }

    virtual void Zero_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
    }

    virtual void Assem_Estimate()
    {
      for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 1.0;
    }


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
    // density, kinematic viscosity, Gen-alpha parameters
    const double rho0, nu, alpha_f, alpha_m, gamma;

    // number of element boundary contiditon
    const int num_ebc_fun;

    // memory layout
    // nLocBas = 4, dof_per_node = 4, vec_size = nLocBas * dof_per_node
    const int nLocBas, dof_per_node, vec_size;
    const int nqp;
    const int snLocBas;

    // stabilization parameters
    const double CI;
    const double CT;

    // Basis function allocations
    double R[4];
    double dR_dx[4];
    double dR_dy[4];
    double dR_dz[4];
    double dxi_dx[9]; // 3x3 matrix

    // Local Assembly local matrix
    // dof_per_node^2 x nLocBas^2
    double Sub_Tan[16][16];


    // functions
    void print_info() const;

    // The metric tensor for tetrahedron needs to be modified
    void get_metric( const double * const &dxi_dx,
        double &G11, double &G12, double &G13,
        double &G22, double &G23, double &G33 ) const;

    void get_tau( double &tau_m_qua, double &tau_c_qua,
        const double &dt, const double * const &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    void get_DC( double &dc_tau, const double * const &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    void get_f(const double &x, const double &y, const double &z,
        const double &t, double &fx, double &fy, double &fz ) const
    {
      /*
         const double pi = MATH_T::PI;

         const double x2 = x * x;
         const double y2 = y * y;
         const double z2 = z * z;

         const double x3 = x2 * x;
         const double y3 = y2 * y;
         const double z3 = z2 * z;

         const double aa = x - 1.0;
         const double bb = y - 1.0;
         const double cc = z - 1.0;
         const double dd = 2*y2 - 3*y + 1.0;

         const double aa2 = aa * aa;
         const double bb2 = bb * bb;
         const double cc2 = cc * cc;

         const double aa3 = aa2 * aa;
         const double cc3 = cc2 * cc;

         const double t2 = t * t;
         */

      fx = 0.0; fy = 0.0; fz = 0.0;

      //fx = t2*pi*cos(pi*x)*sin(pi*y)*sin(pi*z) - nu*(4*t2*x2*y*aa2*dd + 8*t2*x2*y*z*aa2*cc + 4*t2*x2*y*z*cc*dd + 4*t2*y*z*aa2*cc*dd + 4*t2*x2*z*(4*y - 3)*aa2*cc + 8*t2*x*y*z*(2*x - 2)*cc*dd) + 4*t*x2*y*z*aa2*cc*dd - 4*t2*t2*x3*y2*z2*(2*z - 1)*aa3*cc2*dd*dd + 8*t2*t2*x3*y2*z2*aa2*cc2*(2*x2 - 3*x + 1)*dd*dd - 8*t2*t2*x3*y2*z2*(x - z)*aa3*bb2*cc2*(6*y2 - 6*y + 1);

      //fy = nu*(8*t2*x*y2*z*bb2*cc - 8*t2*x*y2*z*aa*bb2 + 8*t2*x*y2*(x - z)*aa*bb2 + 8*t2*y2*z*(x - z)*bb2*cc - 8*t2*x*y2*aa*bb2*cc + 8*t2*y2*z*aa*bb2*cc + 8*t2*x*y2*z*(x - z)*aa*cc + 8*t2*x*z*(x - z)*aa*bb2*cc + 16*t2*x*y*z*(2*y - 2)*(x - z)*aa*cc) + t2*pi*cos(pi*y)*sin(pi*x)*sin(pi*z) - 8*t*x*y2*z*(x - z)*aa*bb2*cc - 8*t2*t2*x2*y3*z2*aa2*bb2*cc2*dd*(x - 2*z - 2*x*z + 3*z2) + 8*t2*t2*x2*y3*z2*aa2*bb2*cc2*dd*(2*x - z + 2*x*z - 3*x2) + 32*t2*t2*x2*y3*z2*(x - z)*(x - z)*aa2*bb2*cc2*dd;

      //fz = nu*(4*t2*y*z2*cc2*dd + 8*t2*x*y*z2*aa*cc2 + 4*t2*x*y*z2*aa*dd + 4*t2*x*y*aa*cc2*dd + 4*t2*x*z2*(4*y - 3)*aa*cc2 + 8*t2*x*y*z*(2*z - 2)*aa*dd) + t2*pi*cos(pi*z)*sin(pi*x)*sin(pi*y) - 4*t*x*y*z2*aa*cc2*dd - 4*t2*t2*x2*y2*z3*(2*x - 1)*aa2*cc3*dd*dd + 8*t2*t2*x2*y2*z3*aa2*cc2*dd*dd*(2*z2 - 3*z + 1) + 8*t2*t2*x2*y2*z3*(x - z)*aa2*bb2*cc3*(6*y2 - 6*y + 1);
    }

    void get_top_H(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
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

    void get_H1(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }

    void get_H2(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }

    typedef void ( PLocAssem_Tet4_VMS_NS_3D_GenAlpha::*locassem_tet4_vms_ns_funs )( const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const;

    locassem_tet4_vms_ns_funs * flist;

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
