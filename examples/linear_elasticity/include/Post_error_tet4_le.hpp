#ifndef POST_ERROR_TET4_LE_HPP
#define POST_ERROR_TET4_LE_HPP
// ==================================================================
// Post_error_tet4_le.hpp
// ------------------------------------------------------------------
// This is the header file for a collection of error calculators at
// the element level.
//
// The exact solution has to be provided by the user.
//
// Date: May 11 2017
// Author: Ju Liu
// ==================================================================
#include "FEAElement.hpp"
#include "AInt_Weight.hpp"
#include "Math_Tools.hpp"
#include "Matrix_3x3.hpp"

namespace POST_T_TET4_LE
{
  // ----------------------------------------------------------------
  // ! exact 3d vector : returns the exact vector solution.
  // \para x, y, z: the spatial coordinates
  // \para val_x/y/z : output
  // ----------------------------------------------------------------
  void exact_disp( const double &x, const double &y, const double &z,
      double &val_x, double &val_y, double &val_z );

  // ----------------------------------------------------------------
  // ! exact 3d matrix : returns 3 by 3 matrix values, e.g. strain or
  //                     stress in 3D.
  // ----------------------------------------------------------------
  void exact_3dmatrix( const double &x, const double &y, 
      const double &z, Matrix_3x3 &val );


  // ----------------------------------------------------------------
  // ! get_manu_disp_error : return the l2 error by comparing the 
  //   numerical solution with the given disp exact solution.
  // ----------------------------------------------------------------
  double get_manu_disp_error(
      const double * const &solu,
      const double * const &solv,
      const double * const &solw,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const AInt_Weight * const &weight,
      double * const &R );

  // ----------------------------------------------------------------
  // ! get_manu_stress_error : return the l2 error of the exact 
  //   matrix value and the numerical Cauchy stress value.
  //
  //   Note: 
  //   1. make sure that the matrix given in the exact_3dmatrix
  //   is the exact Cauchy stress.
  //   
  //   2. make sure that the material parameters are given correctly.
  // ----------------------------------------------------------------
  double get_manu_stress_error(
      const double * const &sol_x,
      const double * const &sol_y,
      const double * const &sol_z,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const AInt_Weight * const &weight,
      double * const &R,
      double * const &dR_dx,
      double * const &dR_dy,
      double * const &dR_dz );

}

#endif
