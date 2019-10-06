#ifndef PLOCASSEM_NS_3D_VMS_GENALPHA_EXACTJACOBIAN_NOCACHE_HPP
#define PLOCASSEM_NS_3D_VMS_GENALPHA_EXACTJACOBIAN_NOCACHE_HPP
// ==================================================================
// PLocAssem_NS_3D_VMS_GenAlpha_ExactJacobian_noCache.hpp
// This is the local assembly routine for the incompressible Navier-
// Stokes equations using the VMS formulation for spatial discretization
// and generalized-alpha method for temporal integraiton.
// 
// In this routine, we use exact Jacobian for the Newton iteration. In 
// the reference paper, they used an inexact Jacobian. Everything else
// should be the same.
//
// This is a NOCACHE version, which does NOT require precomputed 
// quadrature info. The quadrature info will be evaluated within the
// local assembly routines.
// 
// Reference: CMAME 197 (2007) pp. 173-201
//
// Author: Ju Liu
// Date: June 16 2015
// ==================================================================

#include <cmath>
#include "Sys_Tools.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "AInt_Weight.hpp"
#include "IPLocAssem.hpp"
#include "FEAElement.hpp"
#include "FEAElement_NURBS_3D_der1_lap_vms.hpp"
#include "BernsteinBasis_Array.hpp"
#include "IAExtractor.hpp"
#include "IALocal_meshSize.hpp"

class PLocAssem_NS_3D_VMS_GenAlpha_ExactJacobian_noCache : public IPLocAssem
{
  public:
    PLocAssem_NS_3D_VMS_GenAlpha_ExactJacobian_noCache(
        const class TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const double &phy_len_x, const double &phy_len_y,
        const double &phy_len_z, const double &max_hx,
        const double &max_hy, const double &max_hz );

    virtual ~PLocAssem_NS_3D_VMS_GenAlpha_ExactJacobian_noCache();

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
    double nu; // kinematic viscosity

    // Data group 2: Numerical parameters
    double alpha_f, alpha_m, gamma; // Generalized-alpha parameter

    // Data group 3: Memory layout
    int vec_size, nLocBas, dof_per_node;

    int nqp; // number of quadrature points

    // Data group 4: Stabilization parameters
    double CI;
    double CT;

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

    // ----- External Body Force Vector
    void get_f(const double &x, const double &y, const double &z, const double &t,
        double &fx, double &fy, double &fz) const
    {
      const double pi = 3.14159266535898;

      const double x2 = x * x;
      const double y2 = y * y;
      const double z2 = z * z;

      const double x3 = x2 * x;
      const double y3 = y2 * y;
      const double z3 = z2 * z;

      const double aa = x - 1.0;
      const double bb = y - 1.0;
      const double cc = z - 1.0;

      const double aa2 = aa * aa;
      const double bb2 = bb * bb;
      const double cc2 = cc * cc;

      const double aa3 = aa2 * aa;
      const double cc3 = cc2 * cc;

      fx = pi*cos(pi*x)*sin(pi*t)*sin(pi*y)*sin(pi*z) - nu*(4*x2*y*sin(pi*t)*aa2*(2*y2 - 3*y + 1) + 4*x2*y*z*sin(pi*t)*cc*(2*y2 - 3*y + 1) + 4*y*z*sin(pi*t)*aa2*cc*(2*y2 - 3*y + 1) + 4*x2*z*sin(pi*t)*(4*y - 3)*aa2*cc + 8*x2*y*z*sin(pi*t)*aa2*cc + 8*x*y*z*sin(pi*t)*(2*x - 2)*cc*(2*y2 - 3*y + 1)) - 4*x3*y2*z2*sin(pi*t)*sin(pi*t)*(2*z - 1)*aa3*cc2*(2*y2 - 3*y + 1)*(2*y2 -3*y + 1) + 2*pi*x2*y*z*cos(pi*t)*aa2*cc*(2*y2 - 3*y + 1) + 8*x3*y2*z2*sin(pi*t)*sin(pi*t)*aa2*cc2*(2*x2 - 3*x + 1)*(2*y2 - 3*y + 1)*(2*y2 -3*y + 1) - 8*x3*y2*z2*sin(pi*t)*sin(pi*t)*(x - z)*aa3*bb2*cc2*(6*y2 - 6*y + 1);

      fy = 8*nu*sin(pi*t)*(x - z)*(x2*y3*y - 2*x2*y3 + 6*x2*y2*z2 - 6*x2*y2*z + x2*y2 - 6*x2*y*z2 + 6*x2*y*z + x2*z2 - x2*z - 2*x*y3*y*z + 4*x*y3*z - 6*x*y2*z2 + 4*x*y2*z + 6*x*y*z2 - 6*x*y*z - x*z2 + x*z + y3*y*z2 - y3*y - 2*y3*z2 + 2*y3 + y2*z2 - y2) + pi*cos(pi*y)*sin(pi*t)*sin(pi*x)*sin(pi*z) + 32*x2*y3*z2*sin(pi*t)*sin(pi*t)*(x - z)*(x - z)*aa2*bb2*cc2*(2*y2 - 3*y + 1) - 8*x2*y3*z2*sin(pi*t)*sin(pi*t)*aa2*bb2*cc2*(2*y2 - 3*y + 1)*(x - 2*z - 2*x*z + 3*z2) + 8*x2*y3*z2*sin(pi*t)*sin(pi*t)*aa2*bb2*cc2*(2*y2 - 3*y + 1)*(2*x - z + 2*x*z - 3*x2) - 4*pi*x*y2*z*cos(pi*t)*(x - z)*aa*bb2*cc;

      fz = nu*(4*y*z2*sin(pi*t)*cc2*(2*y2 - 3*y + 1) + 4*x*y*z2*sin(pi*t)*(x - 1)*(2*y2 - 3*y + 1) + 4*x*y*sin(pi*t)*aa*cc2*(2*y2 - 3*y + 1) + 4*x*z2*sin(pi*t)*(4*y - 3)*(x - 1)*cc2 + 8*x*y*z2*sin(pi*t)*(x - 1)*cc2 + 8*x*y*z*sin(pi*t)*(2*z - 2)*aa*(2*y2 - 3*y + 1)) + pi*cos(pi*z)*sin(pi*t)*sin(pi*x)*sin(pi*y) - 4*x2*y2*z3*sin(pi*t)*sin(pi*t)*(2*x - 1)*aa2*cc3*(2*y2 - 3*y + 1)*(2*y2 - 3*y + 1) - 2*pi*x*y*z2*cos(pi*t)*aa*cc2*(2*y2 - 3*y + 1) + 8*x2*y2*z3*sin(pi*t)*sin(pi*t)*aa2*cc2*(2*y2 - 3*y + 1)*(2*y2 - 3*y + 1)*(2*z2 - 3*z + 1) + 8*x2*y2*z3*sin(pi*t)*sin(pi*t)*(x - z)*aa2*bb2*cc3*(6*y2 - 6*y + 1);
    }

    void get_mf(const double &x, const double &y, const double &z, const double &t,
        double &f1, double &f2, double &f3, double &f4) const
    {
      f1 = 0.0;
      f2 = 0.0;
      f3 = 0.0;
      f4 = 0.0;
    }

    // Get stabilization parameter
    void get_tau(double &tau_m_qua, double &tau_c_qua, const double &dt,
        const double * const &dxi_dx, const double &u, 
        const double &v, const double &w) const;

    // Print the local assembly info on screen
    void print_info() const;
};

#endif
