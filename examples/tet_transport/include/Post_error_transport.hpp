#ifndef POST_ERROR_TRANSPORT_HPP
#define POST_ERROR_TRANSPORT_HPP

#include "FEAElement.hpp"
#include "IQuadPts.hpp"
#include "Math_Tools.hpp"

namespace POST_ERROR_T
{
  double exact_sol( const double &x, const double &y, const double &z, const double &time );

  double exact_sol_x( const double &x, const double &y, const double &z, const double &time );
  
  double exact_sol_y( const double &x, const double &y, const double &z, const double &time );
  
  double exact_sol_z( const double &x, const double &y, const double &z, const double &time );

  double get_manu_sol_error(
      const double &time,
      const double * const &solu,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad );
}

#endif
