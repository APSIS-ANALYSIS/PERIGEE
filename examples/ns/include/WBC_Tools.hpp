#ifndef WBC_TOOLS_HPP
#define WBC_TOOLS_HPP
// ============================================================================
// WBC_Tools.hpp
// ----------------------------------------------------------------------------
// This is a temporary file for the implementation of weakly enforced Dirichlet
// boundary condition.
// WBC_T namespace contains the temporary functions which are necessary for the 
// formulation or the algorithm. Perhaps many of them will be sent to specific
// classes as member functions, and this file will be removed in future.
// ============================================================================

#include "Sys_Tools.hpp"
#include "Vector_3.hpp"

namespace WBC_T
{
  // ----------------------------------------------------------------
  // ! get_tau_B : Calculate the coefficient tau_B := [u*]^2 / ||u_tan||,
  //               by solving the non-linear equation of [u+]:
  // g([u+]) = [u+] + 0.1108 * (exp(0.4*[u+]) - 1 - 0.4*[u+] - (0.4*[u+])^2 / 2 - (0.4*[u+])^3 / 6) - [y+]
  //         = 0,
  //               where [u+] := ||u_tan|| / [u*],
  //                     [y+] := y * [u*] / mu = y * ||u_tan|| / (mu * [u+]),
  //               according to Spalding's paper in 1961.
  // Input: \para u_tan : the tangential velocity vector relative to the wall.
  //        \para yy    : the distance from the wall i.e. the 'y' in the formulation.
  //        \para fl_mu : the fluid viscosity i.e the 'mu' in the formulation.
  // ----------------------------------------------------------------
  double get_tau_B(const Vector_3 &u_tan, const double &yy, const double &fl_mu)
  {
    // Use Newton-Raphson method to solve g([u+]) = 0.
    // When [u+] > 0 and [y+] > 0, g([u+]) is monotonically increasing, and there is a unique root.
    const double u_t = u_tan.norm2();

    double u_p0 = 0.0;  // [u+]_i

    double u_p = 1.0;   // [u+]_(i+1)

    do
    {
      u_p0 = u_p;

      const double g_0 = u_p0
      + 0.110803158362334 * (std::exp(0.4*u_p0) - 1.0 - 0.4*u_p0 - 0.08*u_p0*u_p0 - 0.032*u_p0*u_p0*u_p0* (1.0/3.0))
      - yy * u_t * (1.0 / (fl_mu * u_p0));    // g([u+]_i)

      const double g_der_0 = 1
      + 0.110803158362334 * (0.4 * std::exp(0.4*u_p0) - 0.4 - 0.16*u_p0 - 0.032*u_p0*u_p0)
      + yy * u_t * (1.0 / (fl_mu * u_p0 * u_p0)); // dg/d[u+] at [u+]_i

      u_p = u_p0 - g_0 * (1.0 / g_der_0);

    } while (std::abs(u_p - u_p0) > 1.0e-13);
    
    return u_t * (1.0 / (u_p * u_p));   //  tau_B = [u*]^2 / ||u_tan|| = ||u_tan|| / [u+]^2
  }

}

#endif