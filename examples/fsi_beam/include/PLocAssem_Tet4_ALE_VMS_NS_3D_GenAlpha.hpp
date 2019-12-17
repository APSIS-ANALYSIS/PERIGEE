#ifndef PLOCASSEM_TET4_ALE_VMS_NS_3D_GENALPHA_HPP
#define PLOCASSEM_TET4_ALE_VMS_NS_3D_GENALPHA_HPP
// ==================================================================
// PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha.hpp
// 
// This is the local assembly routine for ALE-VMS formulation of the
// 3D Navier-Stokes equations with Generalized-alpha for time stepping.
//
// This code assumes the usage of four-node tetradedral element, which
// simplifies the implementation w.r.t. second-order derivatives.
//
// Ref: CM vol 43 3-37, 2008
// Author: Ju Liu
// Date Created: Aug. 2 2017
// 
// The additional shock-capturing term is added to enhance stability.
// Date Modified: Aug. 8 2017
// ==================================================================
#include "IPLocAssem.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas );

    virtual ~PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha();

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

    // vec_a : dot_sol at time step n+alpha_m
    // vec_b : sol at time step n+alpha_f
    // In ALE formulations, the vec_a, dot_sol contains
    //  [mesh velocity; dot_pressure; dot_velocity].
    // the vec_b, sol contains
    //  [mesh displacement; pressure; velocity].
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
    // nLocBas = 4, dof_per_node = 7, vec_size = nLocBas * dofMat
    const int nLocBas, dof_per_node, vec_size;
    const int nqp;
    const int snLocBas;

    // stabilization parameters
    // CI is diffusive correction, for linear elements, it is 36,
    // for quadratic, it is 60.0, for cubic, it is 128.
    const double CI;
    const double CT;

    // Basis function allocations
    double R[4];
    double dR_dx[4];
    double dR_dy[4];
    double dR_dz[4];
    double dxi_dx[9]; // 3x3 matrix

    // current mesh coordinates with length nLocBas (=4 here.)
    double curPt_x[4];
    double curPt_y[4];
    double curPt_z[4];

    // Local Assembly local matrix
    // dofMat^2 x nLocBas^2
    double Sub_Tan[16][16];


    // functions
    void print_info() const;

    // The metric tensor modification routine
    void get_metric( const double * const &dxi_dx,
        double &G11, double &G12, double &G13,
        double &G22, double &G23, double &G33 ) const;

    void get_tau( double &tau_m_qua, double &tau_c_qua,
        const double &dt, const double * const &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    // Calculate the \bar{tau} := (v' G v')^-0.5
    void get_DC( double &dc_tau, const double * const &dxi_dx,
        const double &u, const double &v, const double &w ) const;


    void get_f(const double &x, const double &y, const double &z,
        const double &t, double &fx, double &fy, double &fz ) const
    {
      fx = 0.0; fy = 0.0; fz = 0.0;
    }

    void get_H1(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      const double p0 = -0.0e0;
      gx = p0*nx; gy = p0*ny; gz = p0*nz;
    }

    // Define the Natural boundary condition functions
    typedef void ( PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha::*locassem_tet4_ale_vms_ns_funs )( const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const;

    locassem_tet4_ale_vms_ns_funs * flist;

    void get_ebc_fun( const int &ebc_id,
        const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      return ((*this).*(flist[ebc_id]))(x,y,z,t,nx,ny,nz,gx,gy,gz);
    }


    // Zero the container for tangent matrix
    void Zero_Sub_Tan()
    {
      for(int ii=0; ii<16; ++ii)
      {
        for(int jj=0; jj<nLocBas * nLocBas; ++jj) Sub_Tan[ii][jj] = 0.0;
      }
    }


    // get_currPts : use the solution vector and the ref coordinates
    //               to get the current configuration coordinates
    //               x = X + disp
    //               ept_x/y/z : length is nLocBas (4),
    //               sol : length is nLocBas x dofNum (4*7),
    //               where the first three slots are the mesh displacement.
    void get_currPts( const double * const &ept_x,
        const double * const &ept_y,
        const double * const &ept_z,
        const double * const &sol );


    // get_currBCPts : same function as get_currPts, but the length
    //                 of ept_x/y/z is snLocBas (3) 
    //                 of sol is snLocBas x dofNum (3*7)
    void get_currBCPts( const double * const &ept_x,
        const double * const &ept_y,
        const double * const &ept_z,
        const double * const &sol );
};

#endif
