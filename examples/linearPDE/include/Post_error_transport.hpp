#ifndef POST_ERROR_TRANSPORT_HPP
#define POST_ERROR_TRANSPORT_HPP
// ==================================================================
// Post_error_transport.hpp
// 
// Error analysis namespace for transport equations
//
// Date: Oct. 26 2023
// ==================================================================
#include "FEAElement.hpp"
#include "IQuadPts.hpp"
#include "Math_Tools.hpp"

namespace POST_ERROR_T
{
  double exact_sol( const double &x, const double &y, const double &z, const double &time );

  double exact_sol_dx( const double &x, const double &y, const double &z, const double &time );

  double exact_sol_dy( const double &x, const double &y, const double &z, const double &time );

  double exact_sol_dz( const double &x, const double &y, const double &z, const double &time );

  double get_manu_sol_errorL2(
      const double &time,
      const double * const &solu,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad );

   double get_manu_sol_errorH1(
      const double &time,
      const double * const &solu,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad );
}

#endif
