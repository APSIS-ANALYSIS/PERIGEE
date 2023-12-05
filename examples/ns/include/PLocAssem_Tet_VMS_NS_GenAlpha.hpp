#ifndef PLOCASSEM_TET_VMS_NS_GENALPHA_HPP
#define PLOCASSEM_TET_VMS_NS_GENALPHA_HPP
// ==================================================================
// PLocAssem_Tet_VMS_NS_GenAlpha.hpp
// 
// Parallel Local Assembly routine for VMS and Gen-alpha based NS
// solver.
//
// Author: Ju Liu
// Date: Feb. 10 2020
// ==================================================================
#include "IPLocAssem.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "SymmTensor2_3D.hpp"
#include "Math_Tools.hpp"

class PLocAssem_Tet_VMS_NS_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_Tet_VMS_NS_GenAlpha(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas, const double &in_rho, 
        const double &in_vis_mu, const double &in_beta,
        const double &in_ct = 4.0, const double &in_ctauc = 1.0, 
        const int &elemtype = 501 );

    virtual ~PLocAssem_Tet_VMS_NS_GenAlpha();

    virtual int get_dof() const {return 4;}

    virtual int get_dof_mat() const {return 4;}

    virtual double get_model_para_1() const {return alpha_f;}

    virtual double get_model_para_2() const {return gamma;}

    virtual void Zero_Tangent_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
      for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 0.0;
    }

    virtual void Zero_sur_Tangent_Residual()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Residual[ii] = 0.0;
      for(int ii=0; ii<sur_size*sur_size; ++ii) sur_Tangent[ii] = 0.0;
    }

    virtual void Zero_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
    }

    virtual void Zero_sur_Residual()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Residual[ii] = 0.0;
    }

    virtual void Assem_Estimate()
    {
      for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 1.0;
    }

    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Mass_Residual(
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual double get_flowrate( const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void get_pressure_area( const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad,
        double &pres, double &area );

    virtual void Assem_Residual_EBC_Resistance(
        const int &ebc_id, const double &val,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_BackFlowStab(
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual_BackFlowStab(
        const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_Weak1(
        const double &time, const double &dt,
        const double * const &sol,
        FEAElement * const &elementv,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quads,
        const int &face_id,
        const double &C_bI);

    virtual void Assem_Tangent_Residual_Weak1(
        const double &time, const double &dt,
        const double * const &sol,
        FEAElement * const &elementv,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quads,
        const int &face_id,
        const double &C_bI);

  private:
    // Private data
    const double rho0, vis_mu, alpha_f, alpha_m, gamma, beta;

    const int nqp; // number of quadrature points

    double CI; // Constants for stabilization parameters
    
    const double CT; // Constants for stabilization parameters

    const double Ctauc; // Constant scaling factor for tau_C

    int nLocBas, snLocBas, vec_size, sur_size;

    // Private functions
    void print_info() const;

    SymmTensor2_3D get_metric( const std::array<double, 9> &dxi_dx ) const;

    // Return tau_m and tau_c in RB-VMS
    std::array<double, 2> get_tau( const double &dt, 
        const std::array<double, 9> &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    // Return tau_bar := (v' G v')^-0.5 x rho0, 
    //        which scales like Time x Density
    // Users can refer to Int. J. Numer. Meth. Fluids 2001; 35: 93–116 for more details
    double get_DC( const std::array<double, 9> &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    Vector_3 get_f(const Vector_3 &pt, const double &tt) const
    {
      const double x = pt.x(), y = pt.y(), R = 0.1, v_ave = 3.0, fl_mu = 4.0e-2, fl_density = 1.0;
      const double r = std::sqrt(x*x + y*y);
      const double dp_dz = -1000.0;

      double fz = 5.0 * v_ave * fl_mu * (3 * r) / (R*R*R * fl_density) + dp_dz / fl_density;
      return Vector_3( 0.0, 0.0, fz );
    }

    Vector_3 get_H1(const Vector_3 &pt, const double &tt, 
        const Vector_3 &n_out ) const
    {
      const double p0 = 0.0;
      return Vector_3( p0*n_out.x(), p0*n_out.y(), p0*n_out.z() );
    }

    typedef Vector_3 ( PLocAssem_Tet_VMS_NS_GenAlpha::*locassem_tet_vms_ns_funs )( 
        const Vector_3 &pt, const double &tt, const Vector_3 &n_out ) const;

    locassem_tet_vms_ns_funs * flist;

    Vector_3 get_ebc_fun( const int &ebc_id, const Vector_3 &pt, 
        const double &tt, const Vector_3 &n_out ) const
    {
      return ((*this).*(flist[ebc_id]))(pt, tt, n_out);
    }

    Vector_3 get_Poiseuille_traction(const Vector_3 &pt, const double &tt, const Vector_3 &n_out) const
    {
      const double delta_P = 153.6; // pressure
      const double Length = 1.2; // tube length
      const double fl_mu = 4.0e-2; // viscosity

      // ux = 0, uy = 0, uz = delta_P * (R^2 - x^2 - y^2) / (4 * fl_mu * L);
      const double x = pt.x(), y = pt.y();

      Tensor2_3D velo_grad ( 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 
                             -delta_P * x / (2 * fl_mu * Length),
                             -delta_P * y / (2 * fl_mu * Length),
                             0.0);

      Tensor2_3D velo_grad_T = velo_grad;
      velo_grad_T.transpose();

      Tensor2_3D S = velo_grad + velo_grad_T;

      S *= fl_mu;
      
      return S * n_out;
    }

    Vector_3 get_cubic_velo_traction(const Vector_3 &pt, const double &tt, const Vector_3 &n_out) const
    {
      const double Q = 0.03 * MATH_T::PI; // 0.03 * pi, v_ave = 3
      const double R = 0.1;

      const double A = R * R * MATH_T::PI; // area
      const double fl_mu = 4.0e-2;

      // ux = 0, uy = 0, uz = v_max * (1 - r^3 / R^3);
      const double x = pt.x(), y = pt.y();
      const double v_max = (5 * Q/ (3 * A));
      const double ra = std::sqrt(x * x + y * y);

      Tensor2_3D velo_grad ( 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 
                             -3 * v_max * x * ra / (R * R * R),
                             -3 * v_max * y * ra / (R * R * R),
                             0.0);

      Tensor2_3D velo_grad_T = velo_grad;
      velo_grad_T.transpose();

      Tensor2_3D S = velo_grad + velo_grad_T;

      S *= fl_mu;
      
      return S * n_out;
    }

    Vector_3 get_g_weak(const Vector_3 &pt, const double &tt)
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    // ----------------------------------------------------------------
    // ! get_h_b : Calculate the coefficient h_b for weak BC
    // Input: \para dx_dxi : the inverse Jacobian
    //        \para n_out  : the outward normal 
    // ----------------------------------------------------------------
    double get_h_b(const std::array<double, 9> &dx_dxi, const Vector_3 &n_out)
    {
      const Tensor2_3D inv_Jac (dx_dxi[0], dx_dxi[1], dx_dxi[2],
                                dx_dxi[3], dx_dxi[4], dx_dxi[5],
                                dx_dxi[6], dx_dxi[7], dx_dxi[8]);
      
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
};
#endif
