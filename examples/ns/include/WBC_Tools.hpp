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
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Hex.hpp"
#include "QuadPts_Gauss_Quad.hpp"

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

// ----------------------------------------------------------------
// class QuadPts_on_face : express a triangular quadrature rule on a tet's face
//                         or express a quadrilateral quadrature rule on a hex's face
//                         with volume coordinate
//  
// Input: \para qp_surface : the quadrature rule of a surface element
//        \para face_id    : on which face of the volume element
// ----------------------------------------------------------------
class QuadPts_on_face : public IQuadPts
{ 
  public:
    // The input IQuadPts should be QuadPts_Gauss_Quad or QuadPts_Gauss_Triangle
    QuadPts_on_face(IQuadPts * qp_surface, const int &face_id)
    { 
      num_pts = qp_surface->get_num_quadPts();
      dim = qp_surface->get_dim() + 1;

      if(dim == 4) // tet
      {
        qp.assign(4 * num_pts, 0.0);
        switch (face_id)
        {
          case 0: // r = 0
          {
            for(unsigned int ii{0}; ii < num_pts; ++ii)
            {
              qp[4*ii + 1] = qp_surface->get_qp(3*ii + 0);
              qp[4*ii + 2] = qp_surface->get_qp(3*ii + 1);
              qp[4*ii + 3] = qp_surface->get_qp(3*ii + 2);
            }
          } break;
          case 1: // s = 0
          {
            for(unsigned int ii{0}; ii < num_pts; ++ii)
            {
              qp[4*ii + 0] = qp_surface->get_qp(3*ii + 0);
              qp[4*ii + 2] = qp_surface->get_qp(3*ii + 1);
              qp[4*ii + 3] = qp_surface->get_qp(3*ii + 2);
            }
          } break;
          case 2: // t = 0
          {
            for(unsigned int ii{0}; ii < num_pts; ++ii)
            {
              qp[4*ii + 0] = qp_surface->get_qp(3*ii + 0);
              qp[4*ii + 1] = qp_surface->get_qp(3*ii + 1);
              qp[4*ii + 3] = qp_surface->get_qp(3*ii + 2);
            }
          } break;
          case 3: // u = 0
          {
            for(unsigned int ii{0}; ii < num_pts; ++ii)
            {
              qp[4*ii + 0] = qp_surface->get_qp(3*ii + 0);
              qp[4*ii + 1] = qp_surface->get_qp(3*ii + 1);
              qp[4*ii + 2] = qp_surface->get_qp(3*ii + 2);
            }
          } break;
          default:
            SYS_T::print_fatal("Error: QuadPts_on_face, wrong face id for tet.");
            break;
        }
      }
      else if(dim == 3) // hex
      {
        ; // Unimplemented
      }
      else
        SYS_T::print_fatal("Error: QuadPts_on_face, unsupported IQuadPts.\n");

      qw.resize(num_pts);
      for(unsigned int ii{0}; ii < num_pts; ++ii)
        qw[ii] = qp_surface->get_qw(ii);
    }

    ~QuadPts_on_face(){}

    virtual int get_dim() const {return dim;};

    virtual double get_qp(unsigned int ii, unsigned int comp) const
    {return qp[dim*ii + comp];}

    virtual double get_qw(unsigned int ii) const
    {return qw[ii];}

  private:
    int num_pts;
    std::vector<double> qp {};
    std::vector<double> qw {};

    int dim;
};

#endif