#ifndef PLOCASSEM_NLHEAT_3D_GENALPHA_HPP
#define PLOCASSEM_NLHEAT_3D_GENALPHA_HPP
// ==================================================================
// PLocAssem_NLHeat_3D_GenAlpha.hpp
// This is the local assembly routine for nonlinear heat equation in
// three dimension, using generalized alpha method as time marching
// scheme. See J.A. Cottrell, et. al Isogeometric Analysis, pp 198
// for details.
//
// Date: Dec 3rd 2013
// ==================================================================
#include "Sys_Tools.hpp"
#include "Math_Tools.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "AInt_Weight.hpp"
#include "IPLocAssem.hpp"
#include "FEAElement.hpp"

class PLocAssem_NLHeat_3D_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_NLHeat_3D_GenAlpha(
        const class TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp
        );
    virtual ~PLocAssem_NLHeat_3D_GenAlpha();

    virtual int get_dof() const {return dof_per_node;}

    virtual void Zero_Tangent_Residual();
    virtual void Zero_Residual();

    virtual void Assem_Estimate();
    
    virtual void Assem_Residual(
        double time, double dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const class FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const class AInt_Weight * const &weight );


    virtual void Assem_Tangent_Residual(
        double time, double dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const class FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const class AInt_Weight * const &weight );

    
    virtual void Assem_Mass(
        const class FEAElement * const &element,
        const class AInt_Weight * const &weight );


    virtual void Assem_Mass_Residual(
        const double * const &vec_a,
        const class FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const class AInt_Weight * const &weight );

  private:
    // ! Necessary data for assembly:
    // generalized-alpha method
    double alpha_f, alpha_m, gamma;

    // vec_size = nLocBas * dof_per_node
    int vec_size, nLocBas, dof_per_node;

    // number of quadrature points
    int nqp;

    // ! Dynamic array allocations
    double * R;
    double * dR_dx;
    double * dR_dy;
    double * dR_dz;

    // ! Material properties and external functions:
    // ! define the (nonlinear) conductivity tensor
    void get_k( const double &u, const double &x, const double &y, const double &z,
        double &k11, double &k12, double &k13, double &k21, double &k22,
        double &k23, double &k31, double &k32, double &k33 ) const
    {
      k11 = 1.0 + u;    k12 = 0.0;    k13 = 0.0;
      k21 = 0.0;    k22 = 1.0 + u;    k23 = 0.0;
      k31 = 0.0;    k32 = 0.0;    k33 = 1.0 + u;
    }

    // ! define the derivative the conductivity tensor w.r.t. u
    void get_dk_du( const double &u, const double &x, const double &y, const double &z,
        double &dk11, double &dk12, double &dk13, double &dk21, double &dk22,
        double &dk23, double &dk31, double &dk32, double &dk33 ) const
    {
      dk11 = 1.0;  dk12 = 0.0;  dk13 = 0.0;
      dk21 = 0.0;  dk22 = 1.0;  dk23 = 0.0;
      dk31 = 0.0;  dk32 = 0.0;  dk33 = 1.0;
    }


    // ! define the external heat source
    double get_f( const double &x, const double &y, const double &z, 
        const double &t ) const
    { 
      double pi = MATH_T::PI; 
      const double pi2 = pi * pi;
      double a = sin(pi*x); double a2 = a*a;
      double b = sin(pi*y); double b2 = b*b;
      double c = sin(pi*z); double c2 = c*c;

      double t2 = t*t; double t4 = t2 * t2;
   
      return 2*t*sin(pi*x)*sin(pi*y)*sin(pi*z) - t4*pi2*a2*b2 - t4*pi2*a2*c2 - t4*pi2*b2*c2 + 3*t2*pi2*sin(pi*x)*sin(pi*y)*sin(pi*z) + 6*t4*pi2*a2*b2*c2; 

    } 
};

#endif
