#ifndef GENBC_RESISTANCE_HPP
#define GENBC_RESISTANCE_HPP
// ==================================================================
// GenBC_Resistance.hpp
//
// Resistance boundary condition.
//
// get_m will return resis[ii],
// get_P will return resis[ii] x Q + pres_offset[ii];
//
// Author: Ju Liu
// Date: Aug 21 2019
// ==================================================================
#include "Sys_Tools.hpp"
#include "IGenBC.hpp"

class GenBC_Resistance : public IGenBC
{
  public:
    GenBC_Resistance( const char * const &lpn_filename );

    virtual ~GenBC_Resistance();

    virtual void print_info() const;

    virtual int get_num_ebc() const {return num_ebc;}

    // We do not perform boundary check. Users are responsible to
    // make sure 0 <= ii < num_ebc;
    virtual double get_m( const int &ii, const double &Q ) const
    {
      return resis[ii];
    }

    // We do not perform boundary check. Users are responsible to
    // make sure 0 <= ii < num_ebc;
    virtual double get_P( const int &ii, const double &Q ) const
    {
      return resis[ii] * Q + pres_offset[ii];
    }

    virtual double get_P0( const int &ii ) const
    {
      return P0[ii];
    }

    virtual void reset_initial_sol( const int &ii, const double &in_Q_0,
        const double &in_P_0 )
    {
      Q0[ii] = in_Q_0;
      P0[ii] = in_P_0;
    }

  private:
    int num_ebc; // number of elemental boundary faces
    
    // vector storing the resistance and pressure offset values
    // on the faces, with length num_ebc
    std::vector<double> resis, pres_offset;

    // Vectors storing the Q0 and P0 on each outlet surface
    std::vector<double> Q0, P0;
};

#endif
