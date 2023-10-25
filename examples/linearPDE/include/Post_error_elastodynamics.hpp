#ifndef POST_ERROR_ELASTODYNAMICS_HPP
#define POST_ERROR_ELASTODYNAMICS_HPP

#include "FEAElement.hpp"
#include "IQuadPts.hpp"
#include "Math_Tools.hpp"

namespace POST_ERROR_E
{
  Vector_3 exact_disp( const double &x, const double &y, const double &z, const double &time );

  Tensor2_3D exact_grad_disp( const double &x, const double &y, const double &z, const double &time );

  double get_manu_sol_error(
      const double &time,
      const double * const &solux,
      const double * const &soluy,
      const double * const &soluz,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad );

   double get_manu_sol_errorH1(
      const double &time,
      const double * const &solux,
      const double * const &soluy,
      const double * const &soluz,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad );
}

#endif
