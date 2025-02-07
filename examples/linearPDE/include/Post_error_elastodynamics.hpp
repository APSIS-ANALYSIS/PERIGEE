#ifndef POST_ERROR_ELASTODYNAMICS_HPP
#define POST_ERROR_ELASTODYNAMICS_HPP
// ==================================================================
// Post_error_elastodynamics.hpp
// 
// Error analysis namespace for elastodynamics equations
//
// Date: Oct. 26 2023
// ==================================================================
#include "Tensor2_3D.hpp"
#include "FEAElement.hpp"
#include "IQuadPts.hpp"
#include "Math_Tools.hpp"
#include "Tensor2_3D.hpp"

namespace POST_ERROR_E
{
  const double module_E = 1.0e+7;
  const double nu = 0.4;
  const double lambda = nu * module_E / ((1.0+nu) * (1.0-2.0*nu));
  const double mu = 0.5 * module_E / (1.0+nu);
  const double aa = 1.0e-4;
  const double w = 3.141592653589793;
  
  Vector_3 exact_disp( const double &x, const double &y, const double &z, const double &time );

  Tensor2_3D exact_grad_disp( const double &x, const double &y, const double &z, const double &time );

  double get_manu_sol_errorL2(
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
