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
}

#endif
