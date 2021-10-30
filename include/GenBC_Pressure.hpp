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

class GenBC_Pressure : public IGenBC
{
  public:
    GenBC_Pressure( const char * const &pres_file_name );

    virtual ~GenBC_Pressure();

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

    virtual double get_P0( const int &ii ) const;

    virtual void reset_initial_sol( const int &ii, const double &in_Q_0,
        const double &in_P_0, const double &curr_time, const bool &is_restart );

    virtual void write_0D_sol( const int &curr_index, const double &curr_time ) const {};

  private:
    int num_ebc;

    std::vector< std::vector<double> > coef_a, coef_b;

    std::vector<int> num_of_mode;

    std::vector<double> w, period;
};

#endif
