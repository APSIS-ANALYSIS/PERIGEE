#ifndef PLOCASSEM_TET_CMM_GENALPHA_HPP
#define PLOCASSEM_TET_CMM_GENALPHA_HPP
// ==================================================================
// PLocAssem_Tet_CMM_GenAlpha.hpp
// 
// Parallel Local Assembly routine for VMS and Gen-alpha based
// solver for the CMM type FSI problem.
// ==================================================================
#include <complex_bessel.h>
#include "IPLocAssem.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_Tet_CMM_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_Tet_CMM_GenAlpha(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas, const double &in_rho, 
        const double &in_vis_mu, const double &in_beta,
        const double &in_wall_rho, const double &in_nu,
        const double &in_kappa, const double &in_ctauc = 1.0,
        const int &elemtype = 501 );

    virtual ~PLocAssem_Tet_CMM_GenAlpha();

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

    // **** PRESTRESS TODO: additional arg ele_prestress
    virtual void Assem_Residual_EBC_Wall(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol_wall_disp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &ele_thickness,
        const double * const &ele_youngsmod,
        const IQuadPts * const &quad );

    // **** PRESTRESS TODO: additional arg ele_prestress
    virtual void Assem_Tangent_Residual_EBC_Wall(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol_wall_disp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &ele_thickness,
        const double * const &ele_youngsmod,
        const IQuadPts * const &quad );

    // **** PRESTRESS TODO
    // **** Could replace relevant code in Assem_(Tangent_)Residual_EBC_Wall
    // virtual void get_Wall_CauchyStress(
    //     const double * const &sol_wall_disp,
    //     FEAElement * const &element,
    //     const double * const &eleCtrlPts_x,
    //     const double * const &eleCtrlPts_y,
    //     const double * const &eleCtrlPts_z,
    //     const double * const &ele_thickness,
    //     const double * const &ele_youngsmod,
    //     const IQuadPts * const &quad,
    //     double * const &stress );

  private:
    // Private data
    const double rho0, vis_mu, alpha_f, alpha_m, gamma, beta;

    // wall properties: density, Poisson ratio, shear correction factor (kappa)
    const double rho_w, nu_w, kappa_w; 

    const int nqp; // number of quadrature points

    double CI, CT; // Constants for stabilization parameters

    const double Ctauc; // Constant scaling factor for tau_C

    int nLocBas, snLocBas, vec_size, sur_size;

    std::vector<double> R, dR_dx, dR_dy, dR_dz;

    std::vector<double> d2R_dxx, d2R_dyy, d2R_dzz;

    double dxi_dx[9];

    // Private functions
    void print_info() const;

    void get_metric( const double * const &dxi_dx,
        double &G11, double &G12, double &G13,
        double &G22, double &G23, double &G33 ) const;

    // Return tau_m and tau_c in RB-VMS
    void get_tau( double &tau_m_qua, double &tau_c_qua,
        const double &dt, const double * const &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    // Return tau_bar := (v' G v')^-0.5 x rho0, 
    //        which scales like Time x Density
    void get_DC( double &dc_tau, const double * const &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    // Return body force acting on the fluid domain
    void get_f( const double &x, const double &y, const double &z,
        const double &t, double &fx, double &fy, double &fz ) const
    {
      fx = 0.0; fy = 0.0; fz = 0.0;
    }

    // Return body force acting on the wall domain
    void get_fw( const double &x, const double &y, const double &z,
        const double &t, double &fw_x, double &fw_y, double &fw_z ) const
    {
      fw_x = 0.0; fw_y = 0.0; fw_z = 0.0;
    }

    void get_H1( const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      const double p0 = 0.0;
      gx = p0*nx; gy = p0*ny; gz = p0*nz;
    }

    typedef void ( PLocAssem_Tet_CMM_GenAlpha::*locassem_tet_cmm_funs )( const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const;

    locassem_tet_cmm_funs * flist;

    void get_ebc_fun( const int &ebc_id,
        const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz, double &gx, double &gy, double &gz ) const
    {
      // ==== WOMERSLEY CHANGES BEGIN ====
      const double R_pipe = 0.3;                                                 // pipe radius
      const double omega  = MATH_T::PI * 2.0 / 1.1;                              // freqency
      const std::complex<double> i1(0.0, 1.0);
      const std::complex<double> i1_1d5(-0.707106781186547, 0.707106781186547);
      const auto Omega    = std::sqrt(rho0 * omega / vis_mu) * R_pipe;           // womersley number 
      const auto Lambda   = i1_1d5 * Omega;
      const double r      = std::sqrt(x*x+y*y);                                  // radial coord
      const auto   xi     = Lambda * r / R_pipe;

      const double k0 = -21.0469;                                                // mean pressure gradient
      const std::complex<double> B1(-4.926286624202966e3, -4.092542965905093e3); // pressure Fourier coeff
      const std::complex<double> c1(8.863128942479001e2,   2.978553160539686e1); // wave speed
      const std::complex<double> G1(0.829733473284180,      -0.374935589823809); // elasticity factor

      const auto bes0_xi     = sp_bessel::besselJ(0, xi);
      const auto bes1_xi     = sp_bessel::besselJ(1, xi);
      const auto bes2_xi     = sp_bessel::besselJ(2, xi);
      const auto bes0_Lambda = sp_bessel::besselJ(0, Lambda);

      // axial velo gradient
      const double w_x = k0 * x / (2.0*vis_mu) 
          + std::real( B1 * G1 * i1_1d5 * Omega * x * bes1_xi / (rho0 * c1 * R_pipe * r * bes0_Lambda) * exp(i1*omega*(t-z/c1)) );
      const double w_y = k0 * y / (2.0*vis_mu)
          + std::real( B1 * G1 * i1_1d5 * Omega * y * bes1_xi / (rho0 * c1 * R_pipe * r * bes0_Lambda) * exp(i1*omega*(t-z/c1)) );
      const double w_z = std::real( -i1 * omega * B1 / (rho0 * c1 * c1) * (1.0 - G1 * bes0_xi / bes0_Lambda) * exp(i1*omega*(t-z/c1)) );
      const double p = k0 * z + std::real( B1 * exp(i1*omega*(t-z/c1)) );

      // radial velo
      const double vr = std::real( i1 * omega * R_pipe * B1 / ( 2.0 * rho0 * c1 * c1 )
          * ( r / R_pipe - 2.0 * G1 * bes1_xi / (Lambda * bes0_Lambda) ) * exp(i1*omega*(t-z/c1)) );

      // radial velo gradient
      const auto coef = 1.0 - 2.0 * G1 / bes0_Lambda * ( bes1_xi / xi - bes2_xi );
      const double vr_r = std::real( i1 * omega * B1 / (2.0 * rho0 * c1 * c1) * coef * exp(i1*omega*(t-z/c1)) );
      const double vr_z = std::real( B1 * omega * omega * R_pipe / (2.0 * rho0 * c1 * c1 * c1)
          * ( r / R_pipe - 2.0 * G1 * bes1_xi / (Lambda * bes0_Lambda) ) * exp(i1*omega*(t-z/c1)) );

      // polar to cartesian transformation
      const double theta = std::atan2(y, x);
      const double sin_theta = std::sin(theta);
      const double cos_theta = std::cos(theta);

      const double u_x = cos_theta * cos_theta * vr_r + vr * sin_theta * sin_theta / r;
      const double u_y = sin_theta * cos_theta * ( vr_r - vr / r );
      const double u_z = x / r * vr_z;

      const double v_x = sin_theta * cos_theta * ( vr_r - vr / r );
      const double v_y = sin_theta * sin_theta * vr_r + vr * cos_theta * cos_theta / r;
      const double v_z = y / r * vr_z;

      gx = MATH_T::dot3d(-p + 2.0*vis_mu * u_x,    vis_mu*(u_y + v_x),    vis_mu*(u_z + w_x), nx, ny, nz); 
      gy = MATH_T::dot3d(   vis_mu*(v_x + u_y), -p + 2.0*vis_mu * v_y,    vis_mu*(v_z + w_y), nx, ny, nz); 
      gz = MATH_T::dot3d(   vis_mu*(w_x + u_z),    vis_mu*(w_y + v_z), -p + 2.0*vis_mu * w_z, nx, ny, nz); 

      // return ((*this).*(flist[ebc_id]))(x,y,z,t,nx,ny,nz,gx,gy,gz);
      // ==== WOMERSLEY CHANGES END ====
    }
};

#endif
