#ifndef GENBC_PRESSURE_HPP
#define GENBC_PRESSURE_HPP
// ============================================================================
// GenBC_Pressure.hpp
//
// Pressure boundary condition specified by an analytic function.
//
// Author: Ju Liu
// Date Created: Oct. 30 2021
// ============================================================================
#include "IGenBC.hpp"
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"

class GenBC_Pressure : public IGenBC
{
  public:
    GenBC_Pressure( const std::string &pres_file_name, const double &init_time );

    virtual ~GenBC_Pressure() = default;

    virtual void print_info() const;

    virtual int get_num_ebc() const {return num_ebc;}

    virtual double get_m( const int &ii, const double &dot_Q, const double &Q ) const
    {
      return 0.0;
    }

    virtual double get_n( const int &ii, const double &dot_Q, const double &Q ) const
    {
      return 0.0;
    }

    virtual double get_P( const int &ii, const double &dot_Q, const double &Q,
       const double &time ) const;

    virtual double get_P0( const int &ii ) const
    {
      return P0[ii];
    }

    virtual void reset_initial_sol( const int &ii, const double &in_Q_0,
        const double &in_P_0, const double &curr_time, const bool &is_restart )
    {
      P0[ii] = get_P(ii, 0.0, 0.0, curr_time);
    }

  private:
    int num_ebc;

    // length num_ebc x num_of_mode[ii], 0 <= ii < num_ebc
    std::vector< std::vector<double> > coef_a, coef_b;

    // length num_ebc
    std::vector<int> num_of_mode;

    // length num_ebc
    std::vector<double> w, period;

    // Vectors storing the P0 on each outlet surface, or physically the pressure
    // in the previous time step. length num_ebc
    std::vector<double> P0;
};

#endif
