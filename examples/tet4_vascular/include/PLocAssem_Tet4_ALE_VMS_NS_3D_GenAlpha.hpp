#ifndef PLOCASSEM_TET4_ALE_VMS_NS_3D_GENALPHA_HPP
#define PLOCASSEM_TET4_ALE_VMS_NS_3D_GENALPHA_HPP
// ==================================================================
// PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha.hpp
// 
// This is the local assembly routine for ALE-VMS formulation of the
// 3D Navier-Stokes equations with Generalized-alpha for time stepping.
//
// This routine uses the momentum equation (_mom_), meaning that we 
// have rho in front of u_{,t} term. Hence, the viscosity should be 
// the dynamic viscosity.
//
// Date: Mar. 20 2019
// Author: Ju Liu
// ==================================================================
#include "IPLocAssem.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nqp,
        const double &in_rho, const double &in_vis_mu,
        const double &in_beta );

    virtual ~PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha();

    virtual int get_dof() const {return dof_per_node;}

    virtual int get_dof_mat() const {return 4;}

    virtual double get_model_para_1() const {return alpha_f;}

    virtual double get_model_para_2() const {return gamma;}

    virtual void Zero_Tangent_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
      for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 0.0;
    }

    virtual void Zero_sur_Tangent_Residual()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Residual[ii] = 0.0;
      for(int ii=0; ii<sur_size*sur_size; ++ii) sur_Tangent[ii] = 0.0;
    }

    virtual void Zero_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
    }

    virtual void Zero_sur_Residual()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Residual[ii] = 0.0;
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

    // Calculate the flow rate Q := int_Omega^e v dot n dA
    virtual double get_flowrate( const double * const &vec,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    // Calculate the pressure integrated over the element
    // as well as the area of the element.
    virtual void get_pressure_area( const double * const &vec,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad,
        double &pres, double &area );

    // Assembly the elemental boundary condition
    // val is the value scaling the int w . n dA, and for resistance
    // boundary condition, it is flow_rate x C_resis + p_resis
    virtual void Assem_Residual_EBC_Resistance(
        const int &ebc_id,
        const double &val,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    // Assembly the residual due to the back flow stabilization
    virtual void Assem_Residual_BackFlowStab(
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );
    
    // Assembly the residual and tangent due to the back flow stabilization
    virtual void Assem_Tangent_Residual_BackFlowStab(
        const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

  private:
    // beta is the back flow stabilization parameter
    const double rho0, vis_mu, alpha_f, alpha_m, gamma, beta;

    const int nLocBas, dof_per_node, vec_size, sur_size;
    const int nqp;
    const int snLocBas;

    const double CI, CT;

    // Functions
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
        currPt_x[ii] = ept_x[ii] + sol[7*ii];
        currPt_y[ii] = ept_y[ii] + sol[7*ii+1];
        currPt_z[ii] = ept_z[ii] + sol[7*ii+2];
      }
    }

};

#endif
