#ifndef TIMEMETHOD_GENALPHA_HPP
#define TIMEMETHOD_GENALPHA_HPP
// ==================================================================
// TimeMethod_GenAlpha.hpp
// This is simply a class that gives the time integration scheme
// parameters. 
//
// This class defines the alpha_m, alpha_f, beta, and gamma for the 
// generalized-alpha method.
//
// Chung & Hulbert 1993 showed that 2nd-order accuracy is attained if
//                   gamma = 0.5 - alpha_f + alpha_m;
//                   beta = 0.25 (1 - alpha_f + alpha_m)^2.
//
// Unconditional stability requires
//                   alpha_m >= alpha_f >= 0.5
//
// K. Jansen et al 1999 showed that the above conditions are sufficient for
// 1st-order linear systems (beta condition only pertains to 2nd-order case)
// For 1st-order system, displacement plays no role.
//
// To have strict control over high-frequency damping, alpha_m & 
// alpha_f are parameterized by rho_infty. For 2nd-order system,
//                   alpha_m = 2-rho_inf / (1+rho_inf) 
//                   alpha_f = 1 / (1+rho_inf);
// For 1st-order system,
//                   alpha_m = 0.5 * (3-rho_inf) / (1+rho_inf)
//                   alpha_f = 1 / (1+rho_inf).
//
//
// Reference: Chapter 7, Isogeometric Analysis: Toward Integration of
//            CAD & FEA, J.A. Cottrell, et al. 2009.
// Date: Dec 3rd 2013
// ==================================================================
#include "Sys_Tools.hpp"

class TimeMethod_GenAlpha
{
  public:
    // --------------------------------------------------------------
    // Constructor for the parameters in the Generalized-alpha method
    // for FIRST-order system
    // alpha_m = 0.5 (3-rho) / (1+rho)
    // alpha_f = 1.0 / (1+rho)
    // gamma   = 0.5 + alpha_m - alpha_f
    // beta    =  0.25 * (1.0 - alpha_f + alpha_m)^2 (unused for 1st
    //                                                -order system)
    // This constructor shall be used for 1st-order systems.
    // --------------------------------------------------------------
    TimeMethod_GenAlpha( const double &input_spectral );
    
    // --------------------------------------------------------------
    // Constructor for Generalized-alpha method for FIRST-order system
    // beta    =  0.25 * (1.0 - alpha_f + alpha_m)^2 (unused for 1st
    //                                                -order system)
    // This constructor shall be used for 1st-order systems. The purpose
    // of this constructor is to directly enforce 
    //              alpha_f = alpha_m = gamma = 1.0
    // which gives Backward-Euler method for 1st-order systems.
    // --------------------------------------------------------------
    TimeMethod_GenAlpha( const double &input_alpha_m, 
        const double &input_alpha_f, const double &input_gamma );
    
    // --------------------------------------------------------------
    // Constructor for Generalized-alpha method for either 1st or 2nd
    // order system, with given rho_infty.
    // This constructor gives optimal damping for either 1st or 2nd
    // order systems.
    // --------------------------------------------------------------
    TimeMethod_GenAlpha( const double &input_spectral, const bool &is_2ndorder );   

    ~TimeMethod_GenAlpha() = default;

    // Functions give access to the method parameter
    double get_alpha_m() const {return alpha_m;}
    
    double get_alpha_f() const {return alpha_f;}
    
    double get_gamma() const {return gamma;}
    
    double get_beta() const {return beta;}
    
    double get_rho_infty() const {return rho_infty;}

    bool get_flag() const {return is2nd;}

    void print_info() const;

  private:
    const double rho_infty;
    
    double alpha_m, alpha_f, gamma, beta;
    
    // indicate whether the method is for a 1st-or 2nd-order system
    const bool is2nd;

    // indicate whether the parameters are set from rho_infty
    const bool is_rho_set;
};

#endif
