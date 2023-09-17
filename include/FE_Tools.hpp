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
#include "Vector_3.hpp"

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


  // ============================================================================
  // This is a 3-by-3 matrix class that can calculate LU factorization of the
  // dense matrix. The components are stored in a 1-D array.
  //
  // The array that stores the matrix is mat[9]. Logically, the matrix is 
  //                    
  //                     mat[0], mat[1], mat[2]
  //                     mat[3], mat[4], mat[5]
  //                     mat[6], mat[7], mat[8]
  // 
  // The p[3] array is a pointer for pivoting.
  // ============================================================================

  class Matrix_double_3by3_Array
  {
    public:
      // Defalt constructor: an identity matrix
      Matrix_double_3by3_Array();

      // Explicitly define the matrix components
      Matrix_double_3by3_Array( const double &a11, const double &a12, 
          const double &a13, const double &a21, const double &a22, 
          const double &a23, const double &a31, const double &a32, 
          const double &a33 );

      // Copy constructor
      Matrix_double_3by3_Array( const Matrix_double_3by3_Array &other );

      // Destructor
      ~Matrix_double_3by3_Array();

      // Assignment operator
      Matrix_double_3by3_Array& operator= (const Matrix_double_3by3_Array &input);

      // Parenthesis operator. It allows both access matrix entries as well as
      // assigning values to the entry.
      double& operator()(const int &index) {return mat[index];}

      const double& operator()(const int &index) const {return mat[index];}

      // Generate an identity matrix. Erase all previous values and reset p and
      // invm to default values.
      void gen_id();

      // Generate a matrix with random entries
      // All previous values are erased and p & invm are reset to default.
      void gen_rand(const double &min = -1.0, const double &max = 1.0);

      // Generate a Hilbert matrix
      // all previous values are earsed and p & invm are reset to default.
      void gen_hilb();

      // Perform LU factorization of the matrix object and store the L & U
      // matrix using the mat object.
      // Note: mat is changed after this call.
      void LU_fac();

      // Perform LU solve for the 3 mat x = b equations.
      // LU_fac() has to be called first.
      Vector_3 LU_solve( const Vector_3 &b ) const;

      std::array<double, 3> LU_solve( const std::array<double, 3> &b ) const;

      // Perofrm LU solve for the 3 mat x = b equations
      // LU_fac() has to be called first.
      void LU_solve(const double &b1, const double &b2, const double &b3,
          double &x1, double &x2, double &x3) const;

      // Transpose operation for the matrix 
      void transpose();

      // Inverse of the matrix (based on cofactor). The p and invm are not
      // updated.
      void inverse();

      // determinant
      double det() const;

      // Vector multiplication y = Ax
      // make sure the x y vector has length 3.
      void VecMult( const double * const &x, double * const &y ) const; 

      // Matrix multiplication
      void MatMult( const Matrix_double_3by3_Array &mleft, 
          const Matrix_double_3by3_Array &mright );

      // print mat in matrix format
      void print() const;

      // print mat in matrix format and p & invm0, invm1, invm2
      void print_full() const;

    private:
      // Matrix entries
      double mat[9];

      // Pivoting flag. It is generated in LU-factorization and is used for
      // LU_solve.
      int p[3];

      // Inverse of the diagonal entries. It is generated in LU_fac and is 
      // used for LU_solve.
      double invm0, invm1, invm2;
  };

}

#endif