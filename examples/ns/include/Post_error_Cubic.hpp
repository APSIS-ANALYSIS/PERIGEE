#ifndef POST_ERROR_CUBIC_HPP
#define POST_ERROR_CUBIC_HPP

#include "FEAElement.hpp"
#include "IQuadPts.hpp"
#include "Math_Tools.hpp"

namespace POST_ERROR_C
{
  inline Vector_3 exact_velocity ( const Vector_3 &pt, const double &Qt, const double &RR );
  inline Vector_3 exact_grad_u ( const Vector_3 &pt, const double &Qt, const double &RR );
  inline Vector_3 exact_grad_v ( const Vector_3 &pt, const double &Qt, const double &RR );
  inline Vector_3 exact_grad_w ( const Vector_3 &pt, const double &Qt, const double &RR );
  inline double exact_pressure ( const Vector_3 &pt, const double &dp_dx, const double &Length );

  double get_manu_sol_u_L2(
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &Qt,
      const double &RR );

  double get_manu_sol_u_errorL2(
      const double * const &sol_u,
      const double * const &sol_v,
      const double * const &sol_w,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &Qt,
      const double &RR );

  double get_manu_sol_u_H1(
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &Qt,
      const double &RR );

  double get_manu_sol_u_errorH1(
      const double * const &sol_u,
      const double * const &sol_v,
      const double * const &sol_w,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &Qt,
      const double &RR );

}

#endif