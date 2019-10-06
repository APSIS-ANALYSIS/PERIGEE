#ifndef POST_ERROR_HPP
#define POST_ERROR_HPP
// ==================================================================
// Post_error.hpp
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
// Date: Dec 12 2013
// ==================================================================
#include "FEAElement.hpp"
#include "AInt_Weight.hpp"
#include "Math_Tools.hpp"

namespace POST_T
{
  // ----------------------------------------------------------------
  // ! exact_scalar: returns the scalar exact solution
  // \para x,y,z: the spatial coordinates
  // \para t: time
  // ----------------------------------------------------------------
  double exact_scalar(const double &x, const double &y, const double &z,
      const double &t);
 

  // ----------------------------------------------------------------
  // ! exact_scalar: returns the scalar exact solution
  // \para x,y: the spatial coordinates
  // \para t: time
  // ----------------------------------------------------------------
  double exact_scalar(const double &x, const double &y, const double &t);

  
  // ----------------------------------------------------------------
  // ! exact_3dvector: returns the exact vector solution or scalar
  //                   function's gradient
  // \para x,y,z: the spatial coordinates
  // \para t: time
  // \para val_x, val_y, val_z: output
  // ----------------------------------------------------------------
  void exact_3dvector(const double &x, const double &y, const double &z,
      const double &t, double &val_x, double &val_y, double &val_z );

  
  // ----------------------------------------------------------------
  // ! exact_2dvector: returns the exact vector solution or scalar
  //                   function's gradient in 2D
  // \para x,y: the spatial coordinates
  // \para t: time
  // \para val_x, val_y: output
  // ----------------------------------------------------------------
  void exact_2dvector(const double &x, const double &y,
      const double &t, double &val_x, double &val_y );
  
  
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
  
  
  double get_manu_scalar_l2_error(
      const double * const &sol,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
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
  double get_manu_scalar_h1_error(
      const double * const &sol,
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
  
  double get_manu_scalar_h1_error(
      const double * const &sol,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const AInt_Weight * const &weight,
      double * const &R,
      double * const &dR_dx,
      double * const &dR_dy,
      const double &curr );

  
  // ----------------------------------------------------------------
	// ! get_INSK_mass: return the /int_e \rho dx in the element, viz, the
	//   mass in this element
	// \para sol: the vector of local density values
	// \para element: the element storing basis
	// \para weight: the quadrature weights
	// \para R: the allocation for holding basis quadrature info, delete it 
	//          after usage
  // ----------------------------------------------------------------
	double get_INSK_mass(
      const double * const &sol,
      const FEAElement * const &element,
      const AInt_Weight * const &weight,
      double * const &R );


  // ----------------------------------------------------------------
  // ! get_SINSK_Lyapunov: return the Lyapunov functional of the isothermal
  //   Navier-Stokes-Korteweg problem in one single element.
  //   The formula for this functiona is:
  //        W(\rho) + |\nabla \rho|^2 / (2We) + \rho |u|^2 / 2
  //        W(\rho) = 8*theta * rho * log(\rho / (1-\rho))/27 - rho^2
  // \para theta  : the nondimensional temperature
  // \para We     : the Webner number, for definition, see JCP 248:47-86
  // \para sol_rho: the vector of local density solutions
  // \para sol_u  : the vector of local velocity in x-dir solutions
  // \para sol_v  : the vector of local velocity in y-dir solutions
  // \para sol_w  : the vector of local velocity in z-dir solutions
  // \para element: the element that stores element basis
  // \para weight : the quadrature weights
  // \para R      : the allocation for holding basis quadrature info, users
  //                need to delete it after usage.
  // ----------------------------------------------------------------
  double get_SINSK_Lyapunov(
      const double &theta,
      const double &We,
      const double * const &sol_rho,
      const double * const &sol_u,
      const double * const &sol_v,
      const double * const &sol_w,
      const FEAElement * const &element,
      const AInt_Weight * const &weight,
      double * const &R,
      double * const &R_x,
      double * const &R_y,
      double * const &R_z ); 

   
  // -----------------------------------------------------------------
  // ! get_SINSK_free_energy: return the free energy function W(\rho) in the
  //   Lyapunov functional
  //      W(rho) = 8 * theta * rho * log(rho / (1-rho)) - rho^2
  // -----------------------------------------------------------------
  double get_SINSK_free_energy(
      const double &theta,
      const double * const &sol_rho,
      const FEAElement * const &element,
      const AInt_Weight * const &weight,
      double * const &R );


  // -----------------------------------------------------------------
  // ! get_TNSK_uz_theta : return the volume integral of uz * theta.
  //   This value is used for the evaluation of Nusselt number in thermal
  //   analysis.
  //   users should allocate R for holding quadrature values and delete after
  //   this function call.
  // -----------------------------------------------------------------
  double get_TNSK_uz_theta(
      const double * const &uz_theta,
      const double * const &inv_theta,
      const FEAElement * const &element,
      const AInt_Weight * const &weight,
      double * const &R );

  // -----------------------------------------------------------------
  // ! get_volume : return the physical domain size
  // -----------------------------------------------------------------
  double get_volume( const FEAElement * const &element,
      const AInt_Weight * const &weight );


  // -----------------------------------------------------------------
  // ! get_TNSK_rhou : return the square root of the volume integral of 
  //                   rho u \cdot rho u.
  //   This value is used for the evaluation of Reynolds number in TNSK
  //   analysis.
  //   users should allocate R for holding quadrature values and delete after
  //   this function call.
  // ----------------------------------------------------------------
  double get_TNSK_rhou(
      const double * const &rho_vec,
      const double * const &ux_theta,
      const double * const &uy_theta,
      const double * const &uz_theta,
      const double * const &inv_theta,
      const FEAElement * const &element,
      const AInt_Weight * const &weight,
      double * const &R );

  double get_TNSK_rhou(
      const double * const &rho_vec,
      const double * const &ux_theta,
      const double * const &uy_theta,
      const double * const &inv_theta,
      const FEAElement * const &element,
      const AInt_Weight * const &weight,
      double * const &R );
}

#endif
