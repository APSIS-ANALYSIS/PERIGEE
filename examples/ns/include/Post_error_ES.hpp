#ifndef POST_ERROR_ES_HPP
#define POST_ERROR_ES_HPP

#include "FEAElement.hpp"
#include "IQuadPts.hpp"
#include "Math_Tools.hpp"

namespace POST_ERROR_C
{
  inline Vector_3 exact_velocity ( const double &x, const double &y, const double &z, const double &t );
  inline Vector_3 exact_grad_u ( const double &x, const double &y, const double &z, const double &t );
  inline Vector_3 exact_grad_v ( const double &x, const double &y, const double &z, const double &t );
  inline Vector_3 exact_grad_w ( const double &x, const double &y, const double &z, const double &t );
  inline double exact_pressure ( const Vector_3 &pt, const double &t );

double get_manu_sol_p_L2(
      const double &time,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad );

double get_manu_sol_p_errorL2(
      const double &time,
      const double * const &sol_p,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad );

  double get_manu_sol_u_L2(
      const double &time,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad );
      // const double &Qt,
      // const double &RR );

  double get_manu_sol_u_errorL2(
      const double &time,
      const double * const &sol_u,
      const double * const &sol_v,
      const double * const &sol_w,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad );
      // const double &Qt,
      // const double &RR );

  double get_manu_sol_u_H1(
      const double &time,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad );
      // const double &Qt,
      // const double &RR );

  double get_manu_sol_u_errorH1(
      const double &time,
      const double * const &sol_u,
      const double * const &sol_v,
      const double * const &sol_w,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad );
      // const double &Qt,
      // const double &RR );

}

#endif
