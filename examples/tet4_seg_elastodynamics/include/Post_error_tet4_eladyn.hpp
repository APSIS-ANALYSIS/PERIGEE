#ifndef POST_ERROR_TET4_ELADYN_HPP
#define POST_ERROR_TET4_ELADYN_HPP
// ==================================================================
// Post_error_tet4_eladyn.hpp
// ------------------------------------------------------------------
// This is the header file for a collection of postprocessing
// routines: calculating errors of solutions at element level
// The types of errors include: L2/H1/etc.
// The method include comparing with manufactured solutions, comparing
// with overkill solutions, etc.
//
// The exact solution has to be provided if the user wants to do with
// manufactured solution.
//
// Date: Mar 12 2017
// ==================================================================
#include "FEAElement.hpp"
#include "AInt_Weight.hpp"
#include "Math_Tools.hpp"
#include "Matrix_3x3.hpp"

namespace POST_T_TET4_ELADYN
{
  // ----------------------------------------------------------------
  // ! exact_scalar: returns the scalar exact solution
  // \para x,y,z: the spatial coordinates
  // \para t: time
  // ----------------------------------------------------------------
  double exact_scalar(const double &x, const double &y, const double &z,
      const double &t);
 

  // ----------------------------------------------------------------
  // ! exact_3dvector: returns the exact vector solution or scalar
  //                   function's gradient
  // \para x,y,z: the spatial coordinates
  // \para t: time
  // \para val_x, val_y, val_z: output
  // ----------------------------------------------------------------
  void exact_disp(const double &x, const double &y, const double &z,
      const double &t, double &val_x, double &val_y, double &val_z );

  void exact_velo(const double &x, const double &y, const double &z,
      const double &t, double &val_x, double &val_y, double &val_z );

  
  // ----------------------------------------------------------------
  // ! exact_3dmatrix: returns the exact matrix solution, e.g. strain
  //                   or stress tensors in 3D.
  // ----------------------------------------------------------------
  void exact_3dmatrix( const double &x, const double &y,
      const double &z, const double &t, Matrix_3x3 &val );


  // ----------------------------------------------------------------
  // return the l2 error in local element by comparing with manufactured
  // solution -- scalar case
  // \para sol: the solution vector with this element
  // \para element: the element storing the quadrature of basis
  // \para ectrlPts_x/y/z: element control points' xyz coordinates
  // \para weight: quadrature weights
  // \para R: allocation for basis functions, delete it after usage
  // \para curr: the current time
  // ----------------------------------------------------------------
  double get_manu_scalar_l2_error(
      const double * const &sol,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const AInt_Weight * const &weight,
      double * const &R,
      const double &curr );
  
  
  // ----------------------------------------------------------------
  // return the h1 error in local element by comparing with manufactured
  // solution -- scalar case
  // \para sol: the solution vector with this element
  // \para element: the element storing the quadrature of basis
  // \para ectrlPts_x/y/z: element control points' xyz coordinates
  // \para weight: quadrature weights
  // \para R: allocation for basis functions, delete it after usage
  // \para dR_dx: allocation for basis functions derivatives
  //              delete it after usage
  // \para curr: the current time
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
      double * const &R,
      const double &curr );


  double get_manu_velo_error(
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


  // ----------------------------------------------------------------
  // return the l2 error in local element by comparing with manufactured
  // 3x3 matrix.
  // \para sol_x : the solution vector in this element in x-direction
  // \para sol_y : the solution vector in this element in y-direction
  // \para sol_z : the solution vector in this element in z-direction
  // \para element: the element storing the quadrature of basis
  // \para ectrlPts_x/y/z: element control points' xyz coordinates
  // \para weight: quadrature weights
  // \para R: allocation for basis functions, delete it after usage
  // \para dR_dx: allocation for basis functions derivatives
  //              delete it after usage
  // \para curr: the current time
  // ----------------------------------------------------------------
  double get_manu_matrix_l2_error(
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
      double * const &dR_dz,
      const double &curr );


  double get_manu_sigma_l2_error(
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
      double * const &dR_dz,
      const double &curr );

}

#endif
