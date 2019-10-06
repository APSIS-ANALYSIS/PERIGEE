#ifndef PLOCASSEM_NS_3D_VMS_GENALPHA_EXACTJAC_NOCACHE_ADVFORM_HPP
#define PLOCASSEM_NS_3D_VMS_GENALPHA_EXACTJAC_NOCACHE_ADVFORM_HPP
// ============================================================================
// PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm.hpp
// This is the local assembly routine for the incompressible Navier-
// Stokes equations using the VMS formulation for spatial discretization
// and generalized-alpha method for temporal integraiton.
// 
// In this routine, we use EXACT Jacobian for the Newton iteration. In 
// the reference paper, they used an inexact Jacobian. Everything else
// should be the same.
//
// This is a NOCACHE version, which does NOT require precomputed 
// quadrature info. The quadrature info will be evaluated within the
// local assembly routines.
//
// AdvForm represents the formulation of the nonlinear advection term.
// In the reference paper, the authors used conservation form for the 
// nonlinear term. In this routine, we use the advection term 
//                    (u . grad) u
// for the nonlinear term. This one is simpler and allows Open BC
// implementation. 
//
// In the derivation of the VMS formulation, there is one assumption that 
// u'=0 on the boundary. Hence, the bilinear form in the stabilization terms
// involving u' (aka cross stress and Reynolds stress) can be written either
// in conservative form or advective form. They are identical due to the u'=0 on
// the boundary assumption. This is merely an assumption made in the reference 
// paper.
// 
// Reference: CMAME 197 (2007) pp. 173-201
//
// Author: Ju Liu
// Date: July 22 2015
// ============================================================================

#include <cmath>
#include "QuadPts_Gauss.hpp"
#include "QuadPts_bc0.hpp"
#include "QuadPts_bc1.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "IPLocAssem.hpp"
#include "FEAElement_NURBS_3D_der1_lap_vms.hpp"
#include "FEAElement_NURBS_3D_der1_lap_vms_jac.hpp"

class PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm : public IPLocAssem
{
  public:
    PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm(
        const class TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_sdeg, const int &in_tdeg, const int &in_udeg,
        const int &in_nqpx, const int &in_nqpy, const int &in_nqpz,
        const double &phy_len_x, const double &phy_len_y,
        const double &phy_len_z, const double &max_hx,
        const double &max_hy, const double &max_hz );

    virtual ~PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm();

    virtual int get_dof() const {return dof_per_node;}

    virtual void Zero_Tangent_Residual();
    virtual void Zero_Residual();
    virtual void Assem_Estimate();

    // ! Assembly Routines
    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const int &eindex,
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z,
        const class AInt_Weight * const &weight );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const int &eindex,
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z,
        const class AInt_Weight * const &weight );


    virtual void Assem_Residual_TopFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const int &eindex,
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z );


    virtual void Assem_Residual_BotFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const int &eindex,
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z );


    virtual void Assem_Mass_Residual(
        const double * const &vec_a,
        const int &eindex,
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z,
        const class AInt_Weight * const &weight );

  private:
    // Data group 1: Material properties
    const double nu; // kinematic viscosity

    // Data group 2: Numerical parameters
    const double alpha_f, alpha_m, gamma; // Generalized-alpha parameter

    // Data group 3: Memory layout
    const int nLocBas, dof_per_node, vec_size;
    const int sdeg, tdeg, udeg;

    const int nqp; // number of quadrature points
    const int nqpx, nqpy, nqpz; // number of quadrature points in each direction

    // Data group 4: Stabilization parameters
    const double CI;
    const double CT;

    // Basis function allocations
    double * R;
    double * dR_dx;
    double * dR_dy;
    double * dR_dz;
    double * d2R_dxx;
    double * d2R_dyy;
    double * d2R_dzz;

    double * dxi_dx;

    double ** Sub_Tan_m;
    double ** Sub_Tan_f;
    double ** Sub_Tan_p;

    // ------------------------------------------------------------------------
    // Private functions
    // ------------------------------------------------------------------------
    // ----- External Body Force Vector
    void get_f(const double &x, const double &y, const double &z, const double &t,
        double &fx, double &fy, double &fz) const
    {
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

      fx = t2*pi*cos(pi*x)*sin(pi*y)*sin(pi*z) - nu*(4*t2*x2*y*aa2*dd + 8*t2*x2*y*z*aa2*cc + 4*t2*x2*y*z*cc*dd + 4*t2*y*z*aa2*cc*dd + 4*t2*x2*z*(4*y - 3)*aa2*cc + 8*t2*x*y*z*(2*x - 2)*cc*dd) + 4*t*x2*y*z*aa2*cc*dd - 4*t2*t2*x3*y2*z2*(2*z - 1)*aa3*cc2*dd*dd + 8*t2*t2*x3*y2*z2*aa2*cc2*(2*x2 - 3*x + 1)*dd*dd - 8*t2*t2*x3*y2*z2*(x - z)*aa3*bb2*cc2*(6*y2 - 6*y + 1);

      fy = nu*(8*t2*x*y2*z*bb2*cc - 8*t2*x*y2*z*aa*bb2 + 8*t2*x*y2*(x - z)*aa*bb2 + 8*t2*y2*z*(x - z)*bb2*cc - 8*t2*x*y2*aa*bb2*cc + 8*t2*y2*z*aa*bb2*cc + 8*t2*x*y2*z*(x - z)*aa*cc + 8*t2*x*z*(x - z)*aa*bb2*cc + 16*t2*x*y*z*(2*y - 2)*(x - z)*aa*cc) + t2*pi*cos(pi*y)*sin(pi*x)*sin(pi*z) - 8*t*x*y2*z*(x - z)*aa*bb2*cc - 8*t2*t2*x2*y3*z2*aa2*bb2*cc2*dd*(x - 2*z - 2*x*z + 3*z2) + 8*t2*t2*x2*y3*z2*aa2*bb2*cc2*dd*(2*x - z + 2*x*z - 3*x2) + 32*t2*t2*x2*y3*z2*(x - z)*(x - z)*aa2*bb2*cc2*dd;

      fz = nu*(4*t2*y*z2*cc2*dd + 8*t2*x*y*z2*aa*cc2 + 4*t2*x*y*z2*aa*dd + 4*t2*x*y*aa*cc2*dd + 4*t2*x*z2*(4*y - 3)*aa*cc2 + 8*t2*x*y*z*(2*z - 2)*aa*dd) + t2*pi*cos(pi*z)*sin(pi*x)*sin(pi*y) - 4*t*x*y*z2*aa*cc2*dd - 4*t2*t2*x2*y2*z3*(2*x - 1)*aa2*cc3*dd*dd + 8*t2*t2*x2*y2*z3*aa2*cc2*dd*dd*(2*z2 - 3*z + 1) + 8*t2*t2*x2*y2*z3*(x - z)*aa2*bb2*cc3*(6*y2 - 6*y + 1);
    
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
    }


    void get_mf(const double &x, const double &y, const double &z, const double &t,
        double &f1, double &f2, double &f3, double &f4) const
    {
      f1 = 0.0;
      f2 = 0.0;
      f3 = 0.0;
      f4 = 0.0;
    }


    // This function will return the Cartesian component of the Neumann boundary
    // Input : Spatial coordinates x-y-z, time t, and normal direction of the
    //         boundary surface:  (normal_x; normal_y; normal_z).
    // Output : gnx, gny, gnz on the TOP surface
    // Note: assuming ex = (1;0;0), ey = (0;1;0), ez = (0;0;1).
    //       then we have (sigma_ij) n_j = (-p delta_ij + tau_ij) n_j = g_N,
    //                    n = (0;0;1),
    //                and g_N = gnx ex + gny ey + gnz ez.
    void get_top_gn( const double &x, const double &y, const double &z, 
        const double &t, const double &normal_x, const double &normal_y,
        const double &normal_z, double &gnx, double &gny, double &gnz ) const
    {
      /*
         const double pi = 3.14159266535898;

         const double x2 = x * x;
         const double y2 = y * y;

         const double x3 = x2 * x;

         const double aa = x - 1.0;
         const double bb = y - 1.0;

         const double aa2 = aa * aa;
         const double bb2 = bb * bb;

         gnx = 2.0*nu*y*sin(pi*t)*(2*y2 - 3*y + 1)*( x2*x2 - 2*x3 + x2 );

         gny = -nu*( 4*x*y2*sin(pi*t)*aa2*bb2 );

         gnz = 0.0;
         */

      // const double pi = 3.14159266535898;
      // const double top_pres = 1.0 + sin(pi * t);

      const double top_pres = 0.0;

      gnx = -top_pres * normal_x;
      gny = -top_pres * normal_y;
      gnz = -top_pres * normal_z;
    }

    void get_bot_gn( const double &x, const double &y, const double &z, 
        const double &t, const double &normal_x, const double &normal_y,
        const double &normal_z, double &gnx, double &gny, double &gnz ) const
    {
      /*
         const double pi = 3.14159266535898;

         const double x2 = x * x;
         const double y2 = y * y;

         const double x3 = x2 * x;
         const double x4 = x2 * x2;

         gnx = 2*nu*y*sin(pi*t)*(2*y2 - 3*y + 1)*( x4 - 2*x3 + x2 );

         gny = -nu*(4*x2*y2*sin(pi*t)*(x - 1)*(y - 1)*(y - 1));

         gnz = 0.0;
         */

      double bot_pres;

      /*
         if( z<-5.0 )
         {
         bot_pres = 7.99932e2;
         }
         else
         {
         bot_pres = 7.99932e2;
         }
         */

      bot_pres = 0.0;

      gnx = -bot_pres * normal_x;
      gny = -bot_pres * normal_y;
      gnz = -bot_pres * normal_z; 
    }

    // Get stabilization parameter
    void get_tau(double &tau_m_qua, double &tau_c_qua, const double &dt,
        const double * const &dxi_dx, const double &u, 
        const double &v, const double &w) const;


    // Print the local assembly info on screen
    void print_info() const;

};

#endif
