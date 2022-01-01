#ifndef PLOCASSEM_2x2BLOCK_TET4_ALE_VMS_NS_GENALPHA_HPP
#define PLOCASSEM_2x2BLOCK_TET4_ALE_VMS_NS_GENALPHA_HPP
// ============================================================================
// PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha.hpp
//
// This is a local assembly routine for ALE-VMS formulation of the 3D
// Navier-Stokes equations with Generalized-alpha for time stepping.
//
// Date: Jan 1 2022
// Author: Ju Liu
// ============================================================================
#include "IPLocAssem_2x2Block.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha : public IPLocAssem_2x2Block
{
  public:
    PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nqp,
        const double &in_rho, const double &in_vis_mu,
        const double &in_beta );

    virtual ~PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha();

    virtual int get_dof_mat_0() const {return 3;}

    virtual int get_dof_mat_1() const {return 1;}

    virtual double get_model_para_1() const {return alpha_f;}

    virtual double get_model_para_2() const {return gamma;}

    virtual void Zero_Tangent_Residual()
    {
      for(int ii=0; ii<vec_size_0; ++ii) Residual0[ii] = 0.0;
      for(int ii=0; ii<vec_size_1; ++ii) Residual1[ii] = 0.0;
    
      for(int ii=0; ii<vec_size_0 * vec_size_0; ++ii) Tangent00[ii] = 0.0;
      for(int ii=0; ii<vec_size_0 * vec_size_1; ++ii) Tangent01[ii] = 0.0;
      for(int ii=0; ii<vec_size_1 * vec_size_0; ++ii) Tangent10[ii] = 0.0;
      for(int ii=0; ii<vec_size_1 * vec_size_1; ++ii) Tangent11[ii] = 0.0;
    }

    virtual void Zero_sur_Tangent_Residual()
    {
      for(int ii=0; ii<sur_size_0; ++ii) sur_Residual0[ii] = 0.0;
      for(int ii=0; ii<sur_size_0 * sur_size_0; ++ii) sur_Tangent00[ii] = 0.0;
    }

    virtual void Zero_Residual()
    {
      for(int ii=0; ii<vec_size_0; ++ii) Residual0[ii] = 0.0;
      for(int ii=0; ii<vec_size_1; ++ii) Residual1[ii] = 0.0;
    }

    virtual void Zero_sur_Residual()
    {
      for(int ii=0; ii<sur_size_0; ++ii) sur_Residual0[ii] = 0.0;
    }

    virtual void Assem_Estimate()
    {
      for(int ii=0; ii<vec_size_0 * vec_size_0; ++ii) Tangent00[ii] = 1.0;
      for(int ii=0; ii<vec_size_0 * vec_size_1; ++ii) Tangent01[ii] = 1.0;
      for(int ii=0; ii<vec_size_1 * vec_size_0; ++ii) Tangent10[ii] = 1.0;
      for(int ii=0; ii<vec_size_1 * vec_size_1; ++ii) Tangent11[ii] = 1.0;
    }

    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &dot_disp,
        const double * const &dot_velo,
        const double * const &dot_pres,
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_disp,
        const double * const &dot_velo,
        const double * const &dot_pres,
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Mass_Residual(
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double * const &disp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    // Calculate the flow rate Q := int_Omega^e v dot n dA
    virtual double get_flowrate( 
        const double * const &disp,
        const double * const &velo,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    // Calculate the pressure integrated over the element
    // as well as the area of the element.
    virtual void get_pressure_area( 
        const double * const &disp,
        const double * const &pres,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad,
        double &pressure, double &area );

    // Assembly the elemental boundary condition
    // val is the value scaling the int w . n dA, and for resistance
    // boundary condition, it is flow_rate x C_resis + p_resis
    virtual void Assem_Residual_EBC_Resistance(
        const double &val,
        const double * const &disp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    // Assembly the residual due to the back flow stabilization
    virtual void Assem_Residual_BackFlowStab(
        const double * const &dot_disp,
        const double * const &disp,
        const double * const &velo,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    // Assembly the residual and tangent due to the back flow stabilization
    virtual void Assem_Tangent_Residual_BackFlowStab(
        const double &dt,
        const double * const &dot_disp,
        const double * const &disp,
        const double * const &velo,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

  private:
    const double rho0, vis_mu, alpha_f, alpha_m, gamma, beta, CI, CT;

    const int nLocBas, snLocBas;
    const int vec_size_0, vec_size_1, sur_size_0, sur_size_1, nqp;

    void print_info() const;

    // The metric tensor for tetrahedron needs to be modified.
    // See Pauli dissertation and Whiting, C.H. RPI dissertation
    void get_metric( const double * const &dxi_dx,
        double &G11, double &G12, double &G13,
        double &G22, double &G23, double &G33 ) const;

    // Tau is different from the kinematic tau with rho come
    // into the definition of tau's
    void get_tau( double &tau_m_qua, double &tau_c_qua,
        const double &dt, const double * const &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    // Tau_DC is different from the kinematic definition with
    // a rho in the definition. It scales like Time * Density
    double get_DC( const double * const &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    Vector_3 get_f(const double &x, const double &y, const double &z, const double &t ) const
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    void get_H1(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      const double p0 = 0.0;
      gx = p0*nx; gy = p0*ny; gz = p0*nz;
    }

    // Define Natural BC functions
    typedef void ( PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha::*locassem_2x2block_tet4_ale_vms_ns_funs )( const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const;

    locassem_2x2block_tet4_ale_vms_ns_funs * flist;

    void get_ebc_fun( const int &ebc_id,
        const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      return ((*this).*(flist[ebc_id]))(x,y,z,t,nx,ny,nz,gx,gy,gz);
    }

    // Get the current point coordinates
    void get_currPts( const double * const &ept_x,
        const double * const &ept_y,
        const double * const &ept_z,
        const double * const &sol,
        const int &len,
        double * const &currPt_x,
        double * const &currPt_y,
        double * const &currPt_z ) const
    {
      for(int ii=0; ii<len; ++ii)
      {
        currPt_x[ii] = ept_x[ii] + sol[3*ii];
        currPt_y[ii] = ept_y[ii] + sol[3*ii+1];
        currPt_z[ii] = ept_z[ii] + sol[3*ii+2];
      }
    }
};

#endif
