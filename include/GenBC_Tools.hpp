#ifndef GENBC_TOOLS_HPP
#define GENBC_TOOLS_HPP
// ==================================================================
// GenBC_Tools.hpp
//
// This is a suite of function tools that will be called in the reduced
// modeling of the cardiovascular system
//
// ==================================================================
#include <vector>
#include "Sys_Tools.hpp"

namespace GENBC_T
{
  // ----------------------------------------------------------------
  // ! get_genbc_file_type : read the genbc file and determine what
  //   type of gen bc the file specifies. It will return
  //   0 : unknown type
  //   1 : Resistance
  //   2 : RCR
  //   3 : Inductance
  //   4 : Coronary (including RCR)
  //   5 : Pressure in Fourier modes
  // ----------------------------------------------------------------
  int get_genbc_file_type( const std::string &lpn_filename );

  // ----------------------------------------------------------------
  // Evaluate a Piecewise Cubic Hermite Interpolating Polynomial (PCHIP) 
  // between points x1 and x2 with values f1, f2 and  
  // derivatives d1 and d2 for points xe in [x1, x2] with the number of points
  // ne. PCHIP is given by 
  // f(x) = f1*h00(t) + delta*d1*h10(t) + f2*h01(t) + delta*d2*h11(t),
  // where 
  //   h00(t)=2*t^3-3*t^2+1,
  //   h10(t)=t^3-2*t^2+t,
  //   h01(t)=-2*t^3+3*t^2,
  //   h11(t)=t^3-t^2  
  //   delta=x2-x1,t=(x-x1)/delta. 
  // Outputs: fe with length ne. It stores the values of the Hermite
  // interplation at those xe points.
  // ----------------------------------------------------------------
  void get_cubic_hermite( const double &x1, const double &x2, 
      const double &f1, const double &f2, const double &d1, const double &d2, 
      const int &ne, const std::vector<double> &xe, std::vector<double> &fe );
  
  // ----------------------------------------------------------------
  // Evaluate the derivatives of a Piecewise Cubic Hermite Interpolating 
  // Polynomial (PCHIP) between points x1 and x2 with values f1, f2 and  
  // derivatives d1 and d2 for points xe in [x1, x2] with the number of 
  // points ne. 
  // Outputs: de with length ne. It stores the derivatives of the Hermite
  // interplation at those xe points.
  // ----------------------------------------------------------------
  void get_cubic_hermite_der( const double &x1, const double &x2, 
      const double &f1, const double &f2, const double &d1, const double &d2, 
      const int &ne, const std::vector<double> &xe, std::vector<double> &de );
       
  // ----------------------------------------------------------------
  // This function sets derivatives for user provided data points forming PCHIP.
  // This function is modified from John Burkardt's C++ version of the original 
  // Fortran program by Fred Fritsch under the GNU LGPL license.
  //
  // Input: np is the number of points
  //        xp is the corresponding (e.g. time) point values, with length np
  //        fp is the corresponding (e.g. pressure) point values, with length np
  // Output: dp stores the derivative at xp, with length np
  //
  // Reference:
  //
  // Fred Fritsch, Ralph Carlson,
  //    Monotone Piecewise Cubic Interpolation,
  //    SIAM Journal on Numerical Analysis,
  //    Volume 17, Number 2, April 1980, pages 238-246.
  //
  //    Fred Fritsch, Judy Butland,
  //    A Method for Constructing Local Monotone Piecewise
  //    Cubic Interpolants,
  //    SIAM Journal on Scientific and Statistical Computing,
  //    Volume 5, Number 2, 1984, pages 300-304.
  // ----------------------------------------------------------------
  void set_pchip( const int &np, const std::vector<double> &xp, 
      const std::vector<double> &fp, std::vector<double> &dp );
      
  // ----------------------------------------------------------------
  // This function performs a sign test.
  // This function is modified from John Burkardt's C++ version of the original 
  // Fortran program by Fred Fritsch under the GNU LGPL license.
  //
  // return -1.0, if arg1 and arg2 are of opposite sign.
  // return  0.0, if either argument is zero.
  // return +1.0, if arg1 and arg2 are of the same sign.  
  //
  // The function is to do this without multiplying ARG1 * ARG2, to avoid possible 
  // over/underflow problems.
  // ----------------------------------------------------------------
  double sign_test( const double &arg1, const double &arg2 );
}

#endif
