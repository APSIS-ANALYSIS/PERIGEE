#ifndef SEG_SOL_TOOLS_HPP
#define SEG_SOL_TOOLS_HPP
// ==================================================================
// Seg_Sol_Tools.hpp
// ------------------------------------------------------------------
// SEG_SOL_T namespace contains a suite of frequently used tools for
// the manipulation of solution vectors in the segregated algorithms.
// 
// Date Created: May 21 2017
// ==================================================================
#include "PDNSolution.hpp"
#include "Math_Tools.hpp"

namespace SEG_SOL_T
{
  // ----------------------------------------------------------------
  // ! PlusAiPV: The solution vector sol is assumed to have 7 dofs
  //               [disp_x, disp_y, disp_z, pres, vel_x, vel_y, vel_z].
  //             The PV vector is assumed to have 4 dofs:
  //               [dp, dv_x, dv_y, dv_z].
  //             The updated solution takes the results:
  //                        disp_x += a dv_x
  //                        disp_y += a dv_y
  //                        disp_z += a dv_z
  //                        pres   += b dp
  //                        vel_x  += c dv_x
  //                        vel_y  += c dv_y
  //                        vel_z  += c dv_z.
  //             E.G., taking a = 0, b = -1, c = -1, one gets the
  //             update in the initialization of dot_sol vector using
  //             consistent mass matrix for the Generalized-alpha.
  // ----------------------------------------------------------------
  void PlusAiPV(const double &aa, const double &bb, const double &cc,
      const PDNSolution * const &step,
      PDNSolution * const &sol ); 


  // ----------------------------------------------------------------
  // ! PlusAiUPV: The solution vector sol is assumed to have 7 dofs
  //               [disp_x, disp_y, disp_z, pres, vel_x, vel_y, vel_z].
  //             The UPV vector is assumed to have 7 dofs:
  //               [dd_x, dd_y, dd_z, dp, dv_x, dv_y, dv_z].
  //             The updated solution takes the results:
  //                        disp_x += a dd_x
  //                        disp_y += a dd_y
  //                        disp_z += a dd_z
  //                        pres   += b dp
  //                        vel_x  += c dv_x
  //                        vel_y  += c dv_y
  //                        vel_z  += c dv_z.
  // ----------------------------------------------------------------
  void PlusAiUPV( const double &aa, const double &bb, const double &cc,
      const PDNSolution * const &step, PDNSolution * const &sol );


  // ----------------------------------------------------------------
  // ! PlusAiVPV: The solution vector sol is assumed to have 7 dofs
  //               [disp_x, disp_y, disp_z, pres, vel_x, vel_y, vel_z].
  //             The UPV vector is assumed to have 7 dofs:
  //               [dd_x, dd_y, dd_z, dp, dv_x, dv_y, dv_z].
  //             The updated solution takes the results:
  //                        disp_x += a dv_x
  //                        disp_y += a dv_y
  //                        disp_z += a dv_z
  //                        pres   += b dp
  //                        vel_x  += c dv_x
  //                        vel_y  += c dv_y
  //                        vel_z  += c dv_z.
  // ----------------------------------------------------------------
  void PlusAiVPV( const double &aa, const double &bb, const double &cc,
      const PDNSolution * const &step, PDNSolution * const &sol );


  // ----------------------------------------------------------------
  // ! CheckUV: This routine checks the dot{U}_n+alpha_m - V_n+alpha_f
  //            = 0 for all the consistent Newton-Raphson iterations.
  //            This means that, using the update procedure proposed
  //            in the segregated algorithm, the kinematic equation
  //            is guaranteed to be satisfied.
  //            If the equation is not satisfied up to a tolerance,
  //            an error message will be sent.
  // ----------------------------------------------------------------
  void CheckUV(const PDNSolution * const &dotU, 
      const PDNSolution * const &V );

}

#endif
