#ifndef POST_ERROR_TET4_NS_HPP
#define POST_ERROR_TET4_NS_HPP
// ==================================================================
// Post_error_tet4_ns.hpp
// ------------------------------------------------------------------
// This is a collection of postprocessing -- error calculating at 
// element level.
//
// The exact solution has to be provided in the exact_xxx functions.
//
// Date Created: June 28 2017
// Author: Ju Liu
// ==================================================================
#include "FEAElement.hpp"
#include "AInt_Weight.hpp"
#include "Math_Tools.hpp"
#include "Matrix_3x3.hpp"

namespace POST_T_TET4_NS
{
  double exact_pres(const double &x, const double &y, const double &z,
      const double &t);

  void exact_grad_pres(const double &x, const double &y, const double &z,
      const double &t, double &val_x, double &val_y, double &val_z);

  void exact_velo(const double &x, const double &y, const double &z,
      const double &t, double &val_x, double &val_y, double &val_z);

  double get_pres_l2_error(
      const double * const &sol,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const AInt_Weight * const &weight,
      double * const &R,
      const double &curr );

  double get_pres_h1_error(
      const double * const &sol,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const AInt_Weight * const &weight,
      double * const &R,
      double * const &Rx,
      double * const &Ry,
      double * const &Rz,
      const double &curr );


  double get_velo_l2_error(
      const double * const &solu,
      const double * const &solv,
      const double * const &solw,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const AInt_Weight * const &weight,
      double * const &R,
      const double &curr );

}

#endif
