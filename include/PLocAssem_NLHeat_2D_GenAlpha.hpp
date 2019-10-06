#ifndef PLOCASSEM_NLHEAT_2D_GENALPHA_HPP
#define PLOCASSEM_NLHEAT_2D_GENALPHA_HPP
// ==================================================================
// PLocAssem_NLHeat_2D_GenAlpha.hpp
// This is the local assembly routine for nonlinear heat equation in
// two dimension, using generalized alpha method as time marching
// scheme. See J.A. Cottrell, et. al Isogeometric Analysis, pp 198
// for details.
//
// Date: April 17 2014
// ==================================================================
#include "Sys_Tools.hpp"
#include "Math_Tools.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "AInt_Weight.hpp"
#include "IPLocAssem.hpp"
#include "FEAElement.hpp"

class PLocAssem_NLHeat_2D_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_NLHeat_2D_GenAlpha(
        const class TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp
        );
    virtual ~PLocAssem_NLHeat_2D_GenAlpha();

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

    // ! Material properties and external functions:
    // ! define the (nonlinear) conductivity tensor
    void get_k( const double &u, const double &x, const double &y,
        double &k11, double &k12, double &k21, double &k22 ) const
    {
      k11 = 1.0;    k12 = 0.0;
      k21 = 0.0;    k22 = 1.0;
    }

    // ! define the derivative the conductivity tensor w.r.t. u
    void get_dk_du( const double &u, const double &x, const double &y, 
        double &dk11, double &dk12, double &dk21, double &dk22 ) const
    {
      dk11 = 0.0;  dk12 = 0.0;
      dk21 = 0.0;  dk22 = 0.0;
    }


    // ! define the external heat source
    double get_f( const double &x, const double &y, const double &t ) const
    { 
      const double pi = MATH_T::PI;
      const double a = sin(pi*x);
      const double b = sin(pi*y);
      const double pi2 = pi * pi;

      return ( pi * cos(pi*t) + 2.0*pi2*sin(pi*t) )*a*b;
    } 
};

#endif
