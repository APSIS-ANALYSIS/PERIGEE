#ifndef PLOCASSEM_2x2BLOCK_ALE_VMS_NS_GENALPHA_HPP
#define PLOCASSEM_2x2BLOCK_ALE_VMS_NS_GENALPHA_HPP
// ============================================================================
// PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha.hpp
//
// This is a local assembly routine for ALE-VMS formulation of the 3D
// Navier-Stokes equations with Generalized-alpha for time stepping.
//
// Date: Jan 1 2022
// Author: Ju Liu
// ============================================================================
#include "IPLocAssem_2x2Block.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "SymmTensor2_3D.hpp"
#include "FEAElementFactory.hpp"
#include "QuadPtsFactory.hpp"

class PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha : public IPLocAssem_2x2Block
{
  public:
    PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha(
        const FEType &in_type, const int &in_nqp_v, const int &in_nqp_s,
        const TimeMethod_GenAlpha * const &tm_gAlpha, const double &in_rho,
        const double &in_vis_mu, const double &in_beta );

    virtual ~PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha();

    virtual int get_dof_0() const {return 3;}

    virtual int get_dof_1() const {return 1;}

    virtual int get_nLocBas_0() const {return nLocBas;}
    
    virtual int get_nLocBas_1() const {return nLocBas;}
    
    virtual int get_snLocBas_0() const {return snLocBas;}
    
    virtual int get_snLocBas_1() const {return snLocBas;}

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
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_disp,
        const double * const &dot_velo,
        const double * const &dot_pres,
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    virtual void Assem_Mass_Residual(
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double * const &disp,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    // Calculate the flow rate Q := int_Omega^e v dot n dA
    virtual double get_flowrate( 
        const double * const &disp,
        const double * const &velo,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    // Calculate the pressure integrated over the element
    // as well as the area of the element.
    virtual void get_pressure_area( 
        const double * const &disp,
        const double * const &pres,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        double &pressure, double &area );

    // Assembly the elemental boundary condition
    // val is the value scaling the int w . n dA, and for resistance
    // boundary condition, it is flow_rate x C_resis + p_resis
    virtual void Assem_Residual_EBC_Resistance(
        const double &val,
        const double * const &disp,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    // Assembly the residual due to the back flow stabilization
    virtual void Assem_Residual_BackFlowStab(
        const double * const &dot_disp,
        const double * const &disp,
        const double * const &velo,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    // Assembly the residual and tangent due to the back flow stabilization
    virtual void Assem_Tangent_Residual_BackFlowStab(
        const double &dt,
        const double * const &dot_disp,
        const double * const &disp,
        const double * const &velo,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

  private:
    // Private data
    const FEType elemType;

    const int nqpv, nqps;

    const std::unique_ptr<FEAElement> elementv, elements;

    const std::unique_ptr<IQuadPts> quadv, quads;

    const double rho0, vis_mu, alpha_f, alpha_m, gamma, beta, CI, CT;

    const int nLocBas, snLocBas;
    const int vec_size_0, vec_size_1, sur_size_0;

    // M matrix for tau_m
    //             mm[0], mm[1], mm[2]
    // M = coef *  mm[3], mm[4], mm[5]
    //             mm[6], mm[7], mm[8]
    const double coef;
    const std::array<double, 9> mm; 

    void print_info() const;

    // The metric tensor for tetrahedron needs to be modified.
    // See Pauli dissertation and Whiting, C.H. RPI dissertation
    SymmTensor2_3D get_metric( const std::array<double, 9> &dxi_dx ) const;

    // Tau is different from the kinematic tau with rho come
    // into the definition of tau's
    std::array<double, 2> get_tau(
        const double &dt, const std::array<double, 9> &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    // Tau_DC is different from the kinematic definition with
    // a rho in the definition. It scales like Time * Density
    double get_DC( const std::array<double, 9> &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    Vector_3 get_f(const Vector_3 &pt, const double &tt ) const
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    Vector_3 get_H1(const Vector_3 &pt, const double &tt, const Vector_3 &n_out) const
    {
      const double p0 = 0.0;
      return Vector_3( p0*n_out.x(), p0*n_out.y(), p0*n_out.z() );
    }

    // Define Natural BC functions
    typedef Vector_3 ( PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::*locassem_2x2block_ale_vms_ns_funs )( 
        const Vector_3 &pt, const double &tt, const Vector_3 &n_out ) const;

    locassem_2x2block_ale_vms_ns_funs * flist;

    Vector_3 get_ebc_fun( const int &ebc_id,
        const Vector_3 &pt, const double &tt, const Vector_3 &n_out ) const
    {
      return ((*this).*(flist[ebc_id]))(pt, tt, n_out);
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
