#ifndef FE_TOOLS_HPP
#define FE_TOOLS_HPP
// ============================================================================
// FE_Tools.hpp
// This file defines multiple funcitons that assist the implementation of finite
// element routines.
//
// Date Created: Sep. 14 2023
// ============================================================================
#include <stdio.h>
#include <vector> 
#include "Matrix_double_3by3_Array.hpp"

namespace MATH_T
{
  // ----------------------------------------------------------------
  // Generate outward normal vector from a tangential vector.
  // t : the tangential vector
  // p0 : the starting point of the tangential vector
  // p1 : the interior point 
  // n : the normal vector
  // Algorithm: p1->p0 gives the vector m,
  //            n = m - (m,t) t / (t,t).
  // ----------------------------------------------------------------
  void get_n_from_t( const double &tx, const double &ty, const double &tz,
      const double &p0_x, const double &p0_y, const double &p0_z,
      const double &p1_x, const double &p1_y, const double &p1_z,
      double &nx, double &ny, double &nz );


  // ----------------------------------------------------------------
  // Calculate the circumscribing sphere's centre point and radius
  // of four given points
  // ----------------------------------------------------------------
  void get_tet_sphere_info( const double &x0, const double &x1,
      const double &x2, const double &x3, const double &y0, 
      const double &y1, const double &y2, const double &y3,
      const double &z0, const double &z1, const double &z2, 
      const double &z3, double &x, double &y, double &z, double &r );

  Vector_3 get_tet_sphere_info( const Vector_3 &pt0, const Vector_3 &pt1, 
      const Vector_3 &pt2, const Vector_3 &pt3, double &radius );

  // --------------------------------------------------------------------------
  // Projection operator
  // --------------------------------------------------------------------------
  // L2-projection of a function to a piecewise constant (DGP0 space)
  // f : the function f's value evaluated at nqp quadrature points
  // gwts : gwts = detJac(i) * w(i) the Jacobian for the element and the
  //        quadrature weights.
  // nqp : number of quadrature points
  // return a scalar Prof(f) := int_omega f dx / int_omega 1 dx
  //                          = sum(f * gwts) / sum(gwts)
  // --------------------------------------------------------------------------
  double L2Proj_DGP0( const double * const &f, 
      const double * const &gwts, const int &nqp );

  // --------------------------------------------------------------------------
  // L2-projection of a function to a piecewise linear (DGP1 space) in 2D.
  // f : the function value evaluated at nqp quadrature points
  // gwts : gwts = detJac(i) * w(i) the Jacobian for the element and the weights
  // qp_x : the quadrature points x-coordinates
  // qp_y : the quadrature points y-coordinates
  // nqp : the number of quadrature points
  // output: coeff_0, coeff_x, coeff_y.
  // The projected polynomial is
  //         coeff_0 + coeff_x x + coeff_y y.
  // --------------------------------------------------------------------------
  void L2Proj_DGP1_2D( const double * const &f,
      const double * const &gwts,
      const double * const &qp_x,
      const double * const &qp_y,
      const int &nqp,
      double &coeff_0, double &coeff_x, double &coeff_y );

  // --------------------------------------------------------------------------
  // L2-projection of a function to a piecewise linear (DGP1 space) in 3D.
  // f : the function value evaluated at nqp quadrature points
  // gwts : gwts = detJac(i) * w(i) the Jacobian for the element and the weights
  // qp_x : the quadrature points x-coordinates
  // qp_y : the quadrature points y-coordinates
  // qp_z : the quadrature points z-coordinates
  // nqp : the number of quadrature points
  // output: coeff_0, coeff_x, coeff_y, coeff_z.
  // The projected polynomial is
  //         coeff_0 + coeff_x x + coeff_y y + coeff_z z.
  // --------------------------------------------------------------------------
  void L2Proj_DGP1_3D( const double * const &f,
      const double * const &gwts,
      const double * const &qp_x,
      const double * const &qp_y,
      const double * const &qp_z,
      const int &nqp,
      double &coeff_0, double &coeff_x, double &coeff_y, double &coeff_z );
}

#endif
