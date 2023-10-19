#ifndef IGENBC_HPP
#define IGENBC_HPP
// ==================================================================
// IGenBC.hpp
//
// Interface for GenBC classes
//
// This is a pure virtual class for General BC for outflow boundary
// conditions.
// 
// Author: Ju Liu
// Date: Aug 21 2019
// Reference: J. Liu, et al. The nested block preconditioning technique 
// for the incompressible Navier-Stokes equations with emphasis on 
// hemodynamic simulations.
// ==================================================================
#include "Sys_Tools.hpp"

class IGenBC
{
  public:
    IGenBC() = default;

    virtual ~IGenBC() = default;

    virtual void print_info() const = 0;

    // --------------------------------------------------------------
    // return the number of faces with elemental boundary condtiions.
    // This value is assumed to match the num_ebc in ALocal_EBC. 
    // --------------------------------------------------------------
    virtual int get_num_ebc() const = 0;

    // --------------------------------------------------------------
    // Get the dP/dQ for surface ii
    // for implicit BC's, this value is defined via a difference quotient
    // that is, (get_P(Q+epsilon) - get_P(Q)) / epsilon 
    // for simple models, e.g. resistance type bc, this value is just
    // the resistance value on this bc.
    // --------------------------------------------------------------
    virtual double get_m( const int &ii, const double &dot_Q,
       const double &Q ) const = 0;

    // --------------------------------------------------------------
    // Get the dP/d(dot_Q) for surface ii
    // for implicit BC's, this value is defined via a difference quotient
    // that is, (get_P(dot_Q+epsilon) - get_P(dot_Q)) / epsilon
    // for simple models, e.g. inductance type bc, this value is just
    // the inductance value on this bc.
    // --------------------------------------------------------------
    virtual double get_n( const int &ii, const double &dot_Q,
       const double &Q ) const = 0;

    // --------------------------------------------------------------
    // Get the P value for surface ii, the traction on the surface is
    // modeled as h = P I
    // for resistance bc, for example, this value is
    // Resistance x Q + P_offset
    // --------------------------------------------------------------
    virtual double get_P( const int &ii, const double &dot_Q,
       const double &Q, const double &time ) const = 0;

    // --------------------------------------------------------------
    // Return the pressure at the time step n, which is used as the
    // initial value for ODE integration.
    // For Resistance bc, it is Resistance x Q_previous + P_offset
    // For RCR, it is get_P(ii, Q_previous), which is also stored
    // as Pi_0 + Q_0 x Rp
    // --------------------------------------------------------------
    virtual double get_P0( const int &ii ) const = 0;

    // --------------------------------------------------------------
    // Record solution values as initial conditions for the next time step
    // para ii : the outlet face id, ranging from 0 to num_ebc - 1
    // in_Q_0 : the initial value for the ODE integration of flow rate
    // in_P_0 : the initial value for the ODE integration of averaged Pressure
    // For problems like RCR, it is often convenient to integrate the
    // ODE for Pi, the pressure over the capacitor. So, there can be 
    // more data to be initialized. Check the details of each class
    // implementation.
    // --------------------------------------------------------------
    virtual void reset_initial_sol( const int &ii, const double &in_Q_0,
        const double &in_P_0, const double &curr_time, const bool &is_restart )
    {
      SYS_T::print_fatal("Error: IGenBC::reset_initial_sol is not implemented.\n");
    }    

    // --------------------------------------------------------------
    // Write 0D solutions into a file for restarting the simulation from
    // a previous step. This function is needed by Coronary BC, as not all
    // 0D solutions can be restored from a previous 3D solution.
    // For RCR, resistance, and inductance, this function does nothing.
    // --------------------------------------------------------------
    virtual void write_0D_sol( const int &curr_index, const double &curr_time ) const
    {}
};

#endif
