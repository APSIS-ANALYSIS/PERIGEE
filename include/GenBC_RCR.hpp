#ifndef GENBC_RCR_HPP
#define GENBC_RCR_HPP
// ==================================================================
// GenBC_RCR.hpp
// This files defines the General boundary condition in cardiovascular
// flow simulations.
//
// Pi := P - Rp Q is the pressure over the capacitor
//
// Author: Ju Liu
// Date: Aug. 19 2019
// ==================================================================
#include <vector>
#include "IGenBC.hpp"

class GenBC_RCR : public IGenBC
{
  public:
    GenBC_RCR( const std::string &lpn_filename, const int &in_N, 
        const double &dt3d );

    virtual ~GenBC_RCR() = default; 

    virtual void print_info() const;

    virtual int get_num_ebc() const {return num_ebc;}

    virtual double get_m( const int &ii, const double &in_dot_Q,
       const double &in_Q ) const;

    virtual double get_n( const int &ii, const double &in_dot_Q,
       const double &in_Q ) const
    {
      return 0.0;
    }

    // Obtain P in order to the define the outlet traction for the ii-th 
    // outlet surface
    virtual double get_P( const int &ii, const double &in_dot_Q,
       const double &in_Q, const double &time = 0.0 ) const;

    virtual double get_P0( const int &ii ) const
    {
      return Q0[ii] * Rp[ii] + Pi0[ii] + Pd[ii];
    }

    virtual void reset_initial_sol( const int &ii, const double &in_Q_0,
        const double &in_P_0, const double &curr_time, const bool &is_restart )
    {
      reset_initial_sol( ii, in_Q_0, in_P_0 );
    }

  private:
    const int N;
    
    const double h; // delta t = Nh

    // Parameters used to define difference quotient for get_m.
    const double absTol, relTol;

    int num_ebc;

    // Vectors storing the Rp, C, and Rd values on each outlet faces
    // the length of the following 4 vectors are num_ebc
    std::vector<double> Rp, C, Rd, Pd; // R-C-R parameters

    // Vectors storing the Q0 and Pi0 on each outlet faces
    std::vector<double> Q0, Pi0;

    double F(const int &ii, const double &pi, const double &q) const
    {
      return -1.0 * pi / (Rd[ii] * C[ii]) + q / C[ii]; 
    }
    
    virtual void reset_initial_sol( const int &ii, const double &in_Q_0,
        const double &in_P_0 );
};

#endif
