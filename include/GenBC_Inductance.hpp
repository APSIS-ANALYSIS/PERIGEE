#ifndef GENBC_INDUCTANCE_HPP
#define GENBC_INDUCTANCE_HPP
// ==================================================================
// GenBC_Inductance.hpp
//
// Inductance boundary condition.
//
// get_m will return 0.0;
// get_n will return induct[ii];
// get_P will return induct[ii] x dot_Q + pres_offset[ii];
//
// Author: Ju Liu
// Date: Feb. 20 2019
// ==================================================================
#include <vector>
#include "IGenBC.hpp"

class GenBC_Inductance : public IGenBC
{
  public:
    GenBC_Inductance( const std::string &lpn_filename );

    virtual ~GenBC_Inductance() = default;

    virtual void print_info() const;

    virtual int get_num_ebc() const {return num_ebc;}

    // We do not perform boundary check. Users are responsible to
    // make sure 0 <= ii < num_ebc;
    virtual double get_m( const int &ii, const double &dot_Q,
       const double &Q ) const
    {
      return 0.0;
    }
    
    // We do not perform boundary check. Users are responsible to
    // make sure 0 <= ii < num_ebc;
    virtual double get_n( const int &ii, const double &dot_Q,
       const double &Q ) const
    {
      return induct[ii];
    }

    // We do not perform boundary check. Users are responsible to
    // make sure 0 <= ii < num_ebc;
    virtual double get_P( const int &ii, const double &dot_Q,
       const double &Q, const double &time = 0.0 ) const
    {
      return induct[ii] * dot_Q + pres_offset[ii];
    }

    virtual double get_P0( const int &ii ) const
    {
      return P0[ii];
    }

    virtual void reset_initial_sol( const int &ii, const double &in_Q_0,
        const double &in_P_0, const double &curr_time, const bool &is_restart )
    {
      Q0[ii] = in_Q_0;
      P0[ii] = in_P_0;
    }

  private:
    int num_ebc; // number of elemental boundary faces
    
    // vector storing the inductance and pressure offset values
    // on the faces, with length num_ebc
    std::vector<double> induct, pres_offset;

    // Vectors storing the Q0 and P0 on each outlet surface
    std::vector<double> Q0, P0;
};

#endif
