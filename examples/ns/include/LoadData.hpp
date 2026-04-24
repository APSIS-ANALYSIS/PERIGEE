#ifndef NS_LOADDATA_HPP
#define NS_LOADDATA_HPP

#include "Vector_3.hpp"

namespace LoadData
{
  // --------------------------------------------------------------------------
  // body_force
  //   Body force density per unit mass, i.e. acceleration-like term b(x,t).
  //   It is multiplied by rho0 during assembly where needed.
  // --------------------------------------------------------------------------
  inline Vector_3 body_force( const Vector_3 &pt, const double &tt )
  {
    return Vector_3(0.0, 0.0, 0.0);
  }

  // --------------------------------------------------------------------------
  // ebc_traction
  //   Surface traction h(x,t,n) applied on an EBC face identified by ebc_id.
  //   The return value is the traction vector in physical coordinates.
  // --------------------------------------------------------------------------
  inline Vector_3 ebc_traction( const int &ebc_id, const Vector_3 &pt,
      const double &tt, const Vector_3 &n_out )
  {
    switch(ebc_id)
    {
      case 0:
      {
        const double p0 = 0.0;
        return Vector_3( p0*n_out.x(), p0*n_out.y(), p0*n_out.z() );
      }
      default:
        return Vector_3(0.0, 0.0, 0.0);
    }
  }
}

#endif
