#ifndef PLOCASSEM_VMS_NS_GENALPHA_WEAKBC_HPP
#define PLOCASSEM_VMS_NS_GENALPHA_WEAKBC_HPP
// ==================================================================
// PLocAssem_VMS_NS_GenAlpha_WeakBC.hpp
// 
// Parallel Local Assembly routine for VMS, weakly enforced no-slip
// boundary condition and Gen-alpha based NS solver.
//
// Author: Xuanming Huang
// Date: Feb. 4 2024
// ==================================================================
#include "PLocAssem_VMS_NS_GenAlpha.hpp"
#include "Math_Tools.hpp"

class PLocAssem_VMS_NS_GenAlpha_WeakBC : public PLocAssem_VMS_NS_GenAlpha
{
  public:
    PLocAssem_VMS_NS_GenAlpha_WeakBC(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas, const double &in_rho, 
        const double &in_vis_mu, const double &in_beta,
        const FEType &elemtype, const double &angular,
        const Vector_3 &point_xyz, const Vector_3 &angular_direc,
        const double &in_ct = 4.0, const double &in_ctauc = 1.0,
        const double &in_C_bI = 4.0 );

    virtual ~PLocAssem_VMS_NS_GenAlpha_WeakBC();

    virtual void print_info() const;

    virtual void Assem_Residual_Weak(
        const double &time, const double &dt,
        const double * const &sol,
        const double * const &local_mvelo,
        const double * const &local_mdisp,
        FEAElement * const &elementvs,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quads,
        const int &face_id);

    virtual void Assem_Tangent_Residual_Weak(
        const double &time, const double &dt,
        const double * const &sol,
        const double * const &local_mvelo,
        const double * const &local_mdisp,
        FEAElement * const &elementvs,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quads,
        const int &face_id);

  protected:
    const double C_bI;

    virtual Vector_3 get_g_weak(const Vector_3 &pt, const double &tt)
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    // ----------------------------------------------------------------
    // ! get_h_b : Calculate the coefficient h_b for weak BC
    // Input: \para dxi_dx : the inverse Jacobian
    //        \para n_out  : the outward normal 
    // ----------------------------------------------------------------
    virtual double get_h_b(const std::array<double, 9> &dxi_dx, const Vector_3 &n_out)
    {
      const Tensor2_3D inv_Jac (dxi_dx[0], dxi_dx[1], dxi_dx[2],
                                dxi_dx[3], dxi_dx[4], dxi_dx[5],
                                dxi_dx[6], dxi_dx[7], dxi_dx[8]);
      
      const Vector_3 temp_vec = inv_Jac.VecMult(n_out);

      const double nT_G_n = temp_vec.dot_product(temp_vec);

      return 2.0 / std::sqrt(nT_G_n);
    }

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
    //
    // Ref: Y.Basilevs et al. Isogeometric variational multiscale modeling of wall-bounded turbulent flows 
    //      with weakly enforced boundary conditions on unstretched meshes, CMAME 2010
    // ----------------------------------------------------------------
    virtual double get_tau_B(const Vector_3 &u_tan, const double &yy, const double &fl_mu)
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
};

#endif