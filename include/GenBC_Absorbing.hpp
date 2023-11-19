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
    // ebc_id / initial flow radius / initial thickness
    GenBC_Absorbing( const std::string &lpn_file, const double &solid_E, const double &solid_nu );

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

    virtual double get_P( const int &ii, const double &dot_Q, const double &in_Area,
        const double &time) const;

    virtual double get_P0( const int &ii ) const
    {
      return 0.0;
    }   // Banned

    // Update the current_outlet_area in the nonlinear solver
    virtual void reset_initial_sol( const int &ii, const double &in_Area,
        const double &in_P_0, const double &curr_time, const bool &is_restart )
    {
      P0[ii] = get_P(ii, 0.0, in_Area, 0.0);
    }

  private:
    int num_ebc;

    // Parameter beta of each outlet. length num_ebc
    std::vector<double> para_beta;

    // Initial outlet area A0 of each outlet. length num_ebc
    std::vector<double> initial_outlet_area;

    // P0. length num_ebc
    std::vector<double> P0;
};

#endif