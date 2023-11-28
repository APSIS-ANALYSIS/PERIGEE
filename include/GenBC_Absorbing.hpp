#ifndef GENBC_ABSORBING_HPP
#define GENBC_ABSORBING_HPP
// ============================================================================
// GenBC_Absorbing.hpp
//
// Absorbing pressure boundary condition: P = (sqrt(A) - sqrt(A0)) * beta / PI,
// where A is the outflow area, A0 is the area at t = 0, beta is the parameter.
//
// Assume the outlet surface is a circle, beta = h_s * E / (R^2 * (1 - nu)),
// where R is the radius of outlet and h_s, E, nu is the thickness, Young's modulus
// and the Poisson's ratio of the solid structure, refering to https://doi.org/10.1137/060678439
//
// In this file, the modified version is applied refering to https://doi.org/10.1002/cnm.2756
//
// Author: Xuanming Huang
// Date Created: Nov. 17 2023
// ============================================================================
#include "IGenBC.hpp"
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"

class GenBC_Absorbing : public IGenBC
{
  public:
    // In the lpn_file, the following data should be input for each outlet:
    // ebc_id / initial thickness / inital A0 / steady A_bar / steady P_ref / initial P0
    GenBC_Absorbing( const std::string &lpn_file, const double &solid_E, const double &solid_nu, const double &in_fl_density );

    virtual ~GenBC_Absorbing() = default;

    virtual void print_info() const;

    virtual int get_num_ebc() const
    {
      return num_ebc;
    }

    virtual double get_m( const int &ii, const double &dot_Q, const double &Q ) const
    {
      return 0.0;
    }  // Banned

    virtual double get_n( const int &ii, const double &dot_Q, const double &Q ) const
    {
      return 0.0;
    }  // Banned

    virtual double get_P( const int &ii, const double &dot_Q, const double &Q,
        const double &time) const
    {
      return P[ii];
    }

    virtual double get_P0( const int &ii ) const
    {
      return P0[ii];
    }

    // Update the current_outlet_area in the nonlinear solver
    virtual void reset_initial_sol( const int &ii, const double &in_Q_0,
        const double &in_P_0, const double &curr_time, const bool &is_restart )
    {
      if(!is_restart) // n > 0
        P0[ii] = P[ii];
      else
        ; // Do nothing because P0 has been initialized by initial_P0 from lpn file

      P[ii] = set_P(ii, in_Q_0); // Use Q_n to calculate P_n+1
    }

  private:
    int num_ebc;

    // Parameter ( t*E*pi / ( (1-mu^2) * sqrt(A0) ) ) of each outlet. length num_ebc
    std::vector<double> para_1;

    // Parameter ( sqrt(fl_density/8) / A_bar ) of each outlet. length num_ebc
    std::vector<double> para_2;

    // Steady-state outlet pressure P_ref. length num_ebc
    std::vector<double> P_ref;

    // P_n. length num_ebc
    std::vector<double> P0;

    // P_n+1. length num_ebc
    std::vector<double> P;

    // calculate P_n+1 with Q_n
    virtual double set_P( const int &ii, const double &Q0 ) const;
};

#endif