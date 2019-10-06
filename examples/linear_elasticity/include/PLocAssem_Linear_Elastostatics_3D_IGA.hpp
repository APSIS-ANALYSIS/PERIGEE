#ifndef PLOCASSEM_LINEAR_ELASTOSTATICS_3D_IGA_HPP
#define PLOCASSEM_LINEAR_ELASTOSTATICS_3D_IGA_HPP
// ==================================================================
// PLocAssem_Linear_Elastostatics_3D_IGA.hpp
//
// This is the local assembly routine for the linear static
// elasticity using IGA (i.e., we need to pass extraction operator.)
// 
// This is a NoCache implementation, which implies that the quadrature
// information is not cached in the memory.
//
// Reference:
// 2. T.J.R. Hughes, Linear FEM book, Chapter 2.
//
// Author: Ju Liu
// Date: May 18 2017
// ==================================================================
#include "QuadPts_Gauss.hpp"
#include "QuadPts_bc0.hpp"
#include "QuadPts_bc1.hpp"
#include "IPLocAssem.hpp"

class PLocAssem_Linear_Elastostatics_3D_IGA : public IPLocAssem
{
  public:
    PLocAssem_Linear_Elastostatics_3D_IGA(
        const double &in_mat_E, const double &in_mat_nu,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_sdeg, const int &in_tdeg, const int &in_udeg,
        const int &in_nqpx, const int &in_nqpy, const int &in_nqpz );

    virtual ~PLocAssem_Linear_Elastostatics_3D_IGA();

    virtual int get_dof() const {return dof_per_node;}

    virtual void Zero_Tangent_Residual();

    virtual void Zero_Residual();

    virtual void Assem_Estimate();

    // Assembly Routines
    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
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
        const class AInt_Weight * const &weight ){};

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
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
        FEAElement * const &element,
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
        const class AInt_Weight * const &weight ){};


    
    // Assemble the elements' TOP surface boundary integral
    virtual void Assem_Residual_TopFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element, 
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z );

    // Assemble the elements' BOTTOM surface boundary integral
    virtual void Assem_Residual_BotFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element, 
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z );

    // Assemble the elements' LEFT surface boundary integral
    virtual void Assem_Residual_LefFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element, 
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z );

    // Assemble the elements' RIGHT surface boundary integral
    virtual void Assem_Residual_RigFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element, 
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z );


    // Assemble the elements' FRONT surface boundary integral
    virtual void Assem_Residual_FroFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element, 
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z );


    // Assemble the elements' BACK surface boundary integral
    virtual void Assem_Residual_BacFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element, 
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z );

  private:
    const double E;
    const double nu;
    const double lambda;
    const double mu;
    const double kappa;

    const int nLocBas, dof_per_node, vec_size;
    const int sdeg, tdeg, udeg;

    const int nqp;
    const int nqpx, nqpy, nqpz;

    double * R;
    double * dR_dx;
    double * dR_dy;
    double * dR_dz;

    double ** Sub_Tan;

    // -------------------------------------------------------------- 
    // PRIVATE FUNCTIONS
    // -------------------------------------------------------------- 
    void print_info() const;

    // --------------------------------------------------------------
    // External Forces
    // --------------------------------------------------------------
    // External volumetric body force vector
    void get_f(const double &x, const double &y, const double &z, 
        double &fx, double &fy, double &fz ) const
    {
      const double pi = MATH_T::PI;
      const double lambda = 0.5769231; 
      const double mu = 0.3846154;

      fx = -pi*pi*(2*lambda*cos(pi*x)*cos(2*pi*z)*sin(pi*y) + 2*lambda*cos(2*pi*x)*cos(pi*y)*sin(pi*z) + 2*mu*cos(pi*x)*cos(2*pi*z)*sin(pi*y) + 2*mu*cos(2*pi*x)*cos(pi*y)*sin(pi*z) - lambda*sin(pi*x)*sin(2*pi*y)*sin(pi*z) - 7*mu*sin(pi*x)*sin(2*pi*y)*sin(pi*z)); 

      fy = -pi*pi*(2*lambda*cos(pi*x)*cos(2*pi*y)*sin(pi*z) + 2*lambda*cos(pi*y)*cos(2*pi*z)*sin(pi*x) + 2*mu*cos(pi*x)*cos(2*pi*y)*sin(pi*z) + 2*mu*cos(pi*y)*cos(2*pi*z)*sin(pi*x) - lambda*sin(2*pi*x)*sin(pi*y)*sin(pi*z) - 7*mu*sin(2*pi*x)*sin(pi*y)*sin(pi*z));

      fz = -2*pi*pi*cos(pi*z)*(lambda*cos(pi*x)*cos(pi*y)*sin(pi*x) + lambda*cos(pi*x)*cos(pi*y)*sin(pi*y) + mu*cos(pi*x)*cos(pi*y)*sin(pi*x) + mu*cos(pi*x)*cos(pi*y)*sin(pi*y) - 4*lambda*sin(pi*x)*sin(pi*y)*sin(pi*z) - 10*mu*sin(pi*x)*sin(pi*y)*sin(pi*z));
    }


    // External boundary traction
    // get_top_H : return the top surface's traction vector (gx, gy, gz).
    void get_top_H(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny, 
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }

    // get_bot_H : return the bottom surface's traction vector (gx, gy, gz).
    void get_bot_H(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      gx = 0.0; gy = 0.0; gz = 0.0;
    }

    // get_lef_H : return left surface's traction 
    void get_lef_H(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      const double pi = MATH_T::PI;
      const double lambda = 0.5769231; 
      const double mu = 0.3846154;

      gx = -2*mu*pi*sin(pi*x)*sin(pi*z);

      gy = -pi*sin(2*pi*x)*sin(pi*z)*(lambda + 2*mu);

      gz = -mu*pi*sin(pi*x)*sin(2*pi*z);

    }

    // get_rig_H : return right surface's traction 
    void get_rig_H(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      const double pi = MATH_T::PI;
      const double lambda = 0.5769231; 
      const double mu = 0.3846154;

      gx = 2*mu*pi*sin(pi*x)*sin(pi*z);

      gy = -pi*sin(2*pi*x)*sin(pi*z)*(lambda + 2*mu);

      gz = -mu*pi*sin(pi*x)*sin(2*pi*z);

    }

    // get_fro_H : return front surface's traction 
    void get_fro_H(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      const double pi = MATH_T::PI;
      const double lambda = 0.5769231; 
      const double mu = 0.3846154;

      gx = -pi*sin(2*pi*y)*sin(pi*z)*(lambda + 2*mu);

      gy = 2*mu*pi*sin(pi*y)*sin(pi*z);

      gz = -mu*pi*sin(pi*y)*sin(2*pi*z);

    }

    // get_bac_H : return back surface's traction 
    void get_bac_H(const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      const double pi = MATH_T::PI;
      const double lambda = 0.5769231; 
      const double mu = 0.3846154;

      gx = -pi*sin(2*pi*y)*sin(pi*z)*(lambda + 2*mu); 

      gy = -2*mu*pi*sin(pi*y)*sin(pi*z);

      gz = -mu*pi*sin(pi*y)*sin(2*pi*z);
    }

    void Zero_Sub_Tan()
    {
      for(int ii=0; ii<9; ++ii)
      {
        for(int jj=0; jj<nLocBas * nLocBas; ++jj) Sub_Tan[ii][jj] = 0.0;
      }
    }
};

#endif
