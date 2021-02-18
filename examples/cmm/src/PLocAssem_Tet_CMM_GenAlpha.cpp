#include "PLocAssem_Tet_CMM_GenAlpha.hpp"

PLocAssem_Tet_CMM_GenAlpha::PLocAssem_Tet_CMM_GenAlpha(
    const TimeMethod_GenAlpha * const &tm_gAlpha,
    const int &in_nlocbas, const int &in_nqp,
    const int &in_snlocbas,
    const double &in_rho, const double &in_vis_mu,
    const double &in_beta, const double &in_wall_rho,
    const double &in_nu, const double &in_kappa,
    const double &in_ctauc, const int &elemtype )
: rho0( in_rho ), vis_mu( in_vis_mu ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), beta(in_beta), rho_w(in_wall_rho),
  nu_w(in_nu), kappa_w(in_kappa), nqp(in_nqp), Ctauc( in_ctauc )
{
  if(elemtype == 501)
  {
    // 501 is linear element
    CI = 36.0; CT = 4.0;
    nLocBas = 4; snLocBas = 3;
  }
  else if(elemtype == 502)
  {
    // 502 is quadratic element
    CI = 60.0; CT = 4.0;
    nLocBas = 10; snLocBas = 6;
  }
  else SYS_T::print_fatal("Error: unknown elem type.\n");

  vec_size = nLocBas * 4; // dof_per_node = 4
  sur_size = snLocBas * 4;

  R.resize(nLocBas);
  dR_dx.resize(nLocBas);
  dR_dy.resize(nLocBas);
  dR_dz.resize(nLocBas);
  d2R_dxx.resize(nLocBas);
  d2R_dyy.resize(nLocBas);
  d2R_dzz.resize(nLocBas);

  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  sur_Tangent = new PetscScalar[sur_size * sur_size];
  sur_Residual = new PetscScalar[sur_size];

  Zero_Tangent_Residual();

  Zero_sur_Tangent_Residual();

  print_info();
}


PLocAssem_Tet_CMM_GenAlpha::~PLocAssem_Tet_CMM_GenAlpha()
{
  delete [] Tangent; Tangent = nullptr; 
  delete [] Residual; Residual = nullptr;
  delete [] sur_Tangent; sur_Tangent = nullptr;
  delete [] sur_Residual; sur_Residual = nullptr;
}


void PLocAssem_Tet_CMM_GenAlpha::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("  Three-dimensional Incompressible Navier-Stokes equations: \n");
  if(nLocBas == 4)
    SYS_T::commPrint("  FEM: 4-node Tetrahedral element \n");
  else if(nLocBas == 10)
    SYS_T::commPrint("  FEM: 10-node Tetrahedral element \n");
  else SYS_T::print_fatal("Error: unknown elem type.\n");
  SYS_T::commPrint("  Spatial: Residual-based VMS \n");
  SYS_T::commPrint("  Temporal: Generalized-alpha Method \n");
  SYS_T::commPrint("  Density rho = %e \n", rho0);
  SYS_T::commPrint("  Dynamic Viscosity mu = %e \n", vis_mu);
  SYS_T::commPrint("  Kienmatic Viscosity nu = %e \n", vis_mu / rho0);
  SYS_T::commPrint("  Stabilization para CI = %e \n", CI);
  SYS_T::commPrint("  Stabilization para CT = %e \n", CT);
  SYS_T::commPrint("  Scaling factor for tau_C = %e \n", Ctauc);
  SYS_T::commPrint("  Backflow Stab. para beta = %e \n", beta);
  SYS_T::commPrint("  Note: \n");
  SYS_T::commPrint("  1. Consistent tangent matrix used. \n");
  SYS_T::commPrint("  2. Nonlinear quadratic term is in advective form. \n");
  SYS_T::commPrint("  3. Pressure is evaluated at n+alpha_f rather than n+1. \n");
  SYS_T::commPrint("----------------------------------------------------------- \n");
}


void PLocAssem_Tet_CMM_GenAlpha::get_metric(
    const double * const &f,
    double &G11, double &G12, double &G13,
    double &G22, double &G23, double &G33 ) const
{
  // PHASTA definition 
  const double coef = 0.6299605249474365;
  const double diag = 2.0;
  const double offd = 1.0;

  const double fk0 = diag * f[0] + offd * (f[3] + f[6]);
  const double fk1 = diag * f[3] + offd * (f[0] + f[6]);
  const double fk2 = diag * f[6] + offd * (f[0] + f[3]);
  const double fk3 = diag * f[1] + offd * (f[4] + f[7]);
  const double fk4 = diag * f[4] + offd * (f[1] + f[7]);
  const double fk5 = diag * f[7] + offd * (f[1] + f[4]);
  const double fk6 = diag * f[2] + offd * (f[5] + f[8]);
  const double fk7 = diag * f[5] + offd * (f[2] + f[8]);
  const double fk8 = diag * f[8] + offd * (f[2] + f[5]);

  G11 = coef * ( fk0 * f[0] + fk1 * f[3] + fk2 * f[6] );
  G12 = coef * ( fk0 * f[1] + fk1 * f[4] + fk2 * f[7] );
  G13 = coef * ( fk0 * f[2] + fk1 * f[5] + fk2 * f[8] );
  G22 = coef * ( fk3 * f[1] + fk4 * f[4] + fk5 * f[7] );
  G23 = coef * ( fk3 * f[2] + fk4 * f[5] + fk5 * f[8] );
  G33 = coef * ( fk6 * f[2] + fk7 * f[5] + fk8 * f[8] );
}


void PLocAssem_Tet_CMM_GenAlpha::get_tau(
    double &tau_m_qua, double &tau_c_qua,
    const double &dt, const double * const &dxi_dx,
    const double &u, const double &v, const double &w ) const
{
  // Use K matrix to correct the metric
  double G11, G12, G13, G22, G23, G33;
  get_metric( dxi_dx, G11, G12, G13, G22, G23, G33 );

  const double GdG = G11 * G11 + 2.0 * G12 * G12 + 2.0 * G13 * G13
    + G22 * G22 + 2.0 * G23 * G23 + G33 * G33;

  const double uGu = G11 * u * u + 2.0 * G12 * u * v + 2.0 * G13 * u * w
    + G22 * v * v + 2.0 * G23 * v * w + G33 * w * w;

  const double g_dot_g = G11 + G22 + G33;

  const double temp_nu = vis_mu / rho0;

  const double denom_m = CT / (dt*dt) + uGu + CI * temp_nu * temp_nu * GdG;

  tau_m_qua = 1.0 / ( rho0 * sqrt(denom_m) );

  const double denom_c = tau_m_qua * g_dot_g;

  tau_c_qua = Ctauc / denom_c;
}


void PLocAssem_Tet_CMM_GenAlpha::get_DC(
    double &dc_tau, const double * const &dxi_dx,
    const double &u, const double &v, const double &w ) const
{
  //double G11, G12, G13, G22, G23, G33;
  //get_metric( dxi_dx, G11, G12, G13, G22, G23, G33 );

  //dc_tau = G11 * u * u + 2.0 * G12 * u * v + 2.0 * G13 * u * w + G22 * v * v
  //  + 2.0 * G23 * v * w + G33 * w * w;

  //if(dc_tau > 1.0e-15) dc_tau = rho0 * std::pow(dc_tau, -0.5);
  //else dc_tau = 0.0;

  dc_tau = 0.0;
}


void PLocAssem_Tet_CMM_GenAlpha::Assem_Residual(
    const double &time, const double &dt,
    const double * const &dot_sol,
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  double f1, f2, f3;

  double tau_m, tau_c, tau_dc;

  const double two_mu = 2.0 * vis_mu;

  const double curr = time + alpha_f * dt;

  Zero_Residual();

  for(int qua=0; qua<nqp; ++qua)
  {
    double u = 0.0, u_t = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v = 0.0, v_t = 0.0, v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w = 0.0, w_t = 0.0, w_x = 0.0, w_y = 0.0, w_z = 0.0;
    double p = 0.0, coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;
    double p_x = 0.0, p_y = 0.0, p_z = 0.0;
    double u_xx = 0.0, u_yy = 0.0, u_zz = 0.0;
    double v_xx = 0.0, v_yy = 0.0, v_zz = 0.0;
    double w_xx = 0.0, w_yy = 0.0, w_zz = 0.0;

    element->get_3D_R_gradR_LaplacianR( qua, &R[0], &dR_dx[0], 
        &dR_dy[0], &dR_dz[0], &d2R_dxx[0], &d2R_dyy[0], &d2R_dzz[0] );

    element->get_invJacobian( qua, dxi_dx );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int ii4 = 4 * ii;

      u_t += dot_sol[ii4+1] * R[ii];
      v_t += dot_sol[ii4+2] * R[ii];
      w_t += dot_sol[ii4+3] * R[ii];

      u += sol[ii4+1] * R[ii];
      v += sol[ii4+2] * R[ii];
      w += sol[ii4+3] * R[ii];
      p += sol[ii4+0] * R[ii];

      u_x += sol[ii4+1] * dR_dx[ii];
      v_x += sol[ii4+2] * dR_dx[ii];
      w_x += sol[ii4+3] * dR_dx[ii];
      p_x += sol[ii4+0] * dR_dx[ii];

      u_y += sol[ii4+1] * dR_dy[ii];
      v_y += sol[ii4+2] * dR_dy[ii];
      w_y += sol[ii4+3] * dR_dy[ii];
      p_y += sol[ii4+0] * dR_dy[ii];

      u_z += sol[ii4+1] * dR_dz[ii];
      v_z += sol[ii4+2] * dR_dz[ii];
      w_z += sol[ii4+3] * dR_dz[ii];
      p_z += sol[ii4+0] * dR_dz[ii];

      u_xx += sol[ii4+1] * d2R_dxx[ii];
      u_yy += sol[ii4+1] * d2R_dyy[ii];
      u_zz += sol[ii4+1] * d2R_dzz[ii];

      v_xx += sol[ii4+2] * d2R_dxx[ii];
      v_yy += sol[ii4+2] * d2R_dyy[ii];
      v_zz += sol[ii4+2] * d2R_dzz[ii];

      w_xx += sol[ii4+3] * d2R_dxx[ii];
      w_yy += sol[ii4+3] * d2R_dyy[ii];
      w_zz += sol[ii4+3] * d2R_dzz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    // Get the tau_m and tau_c
    get_tau(tau_m, tau_c, dt, dxi_dx, u, v, w);

    const double tau_m_2 = tau_m * tau_m;

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    // Get the body force
    get_f(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    const double u_lap = u_xx + u_yy + u_zz;
    const double v_lap = v_xx + v_yy + v_zz;
    const double w_lap = w_xx + w_yy + w_zz;

    const double rx = rho0 * ( u_t + u_x * u + u_y * v + u_z * w - f1 ) + p_x - vis_mu * u_lap;
    const double ry = rho0 * ( v_t + v_x * u + v_y * v + v_z * w - f2 ) + p_y - vis_mu * v_lap;
    const double rz = rho0 * ( w_t + w_x * u + w_y * v + w_z * w - f3 ) + p_z - vis_mu * w_lap;

    const double div_vel = u_x + v_y + w_z;

    const double u_prime = -1.0 * tau_m * rx;
    const double v_prime = -1.0 * tau_m * ry;
    const double w_prime = -1.0 * tau_m * rz;

    const double r_dot_gradu = u_x * rx + u_y * ry + u_z * rz;
    const double r_dot_gradv = v_x * rx + v_y * ry + v_z * rz;
    const double r_dot_gradw = w_x * rx + w_y * ry + w_z * rz;

    // Get the Discontinuity Capturing tau
    get_DC( tau_dc, dxi_dx, u_prime, v_prime, w_prime );

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];
      const double velo_dot_gradR = NA_x * u + NA_y * v + NA_z * w;
      const double r_dot_gradR = NA_x * rx + NA_y * ry + NA_z * rz;
      const double velo_prime_dot_gradR = NA_x * u_prime + NA_y * v_prime + NA_z * w_prime;

      Residual[4*A] += gwts * ( NA * div_vel + tau_m * r_dot_gradR );

      Residual[4*A+1] += gwts * ( NA * rho0 * u_t
          + NA * rho0 * (u * u_x + v * u_y + w * u_z)
          - NA_x * p
          + NA_x * two_mu * u_x
          + NA_y * vis_mu * (u_y + v_x)
          + NA_z * vis_mu * (u_z + w_x)
          + velo_dot_gradR * tau_m * rho0 * rx
          - NA * tau_m * rho0 * r_dot_gradu
          + NA_x * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * rho0 * rx
          + velo_prime_dot_gradR * tau_dc 
          * (u_prime * u_x + v_prime * u_y + w_prime * u_z)
          - NA * rho0 * f1 );

      Residual[4*A+2] += gwts * ( NA * rho0 * v_t
          + NA * rho0 * (u * v_x + v * v_y + w * v_z)
          - NA_y * p
          + NA_x * vis_mu * (u_y + v_x)
          + NA_y * two_mu * v_y
          + NA_z * vis_mu * (v_z + w_y)
          + velo_dot_gradR * tau_m * rho0 * ry
          - NA * tau_m * rho0 * r_dot_gradv
          + NA_y * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * rho0 * ry
          + velo_prime_dot_gradR * tau_dc
          * (u_prime * v_x + v_prime * v_y + w_prime * v_z)
          - NA * rho0 * f2 );

      Residual[4*A+3] += gwts * (NA * rho0 * w_t
          + NA * rho0 * (u * w_x + v * w_y + w * w_z)
          - NA_z * p
          + NA_x * vis_mu * (u_z + w_x)
          + NA_y * vis_mu * (w_y + v_z)
          + NA_z * two_mu * w_z
          + velo_dot_gradR * tau_m * rho0 * rz
          - NA * tau_m * rho0 * r_dot_gradw
          + NA_z * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * rho0 * rz
          + velo_prime_dot_gradR * tau_dc
          * (u_prime * w_x + v_prime * w_y + w_prime * w_z)
          - NA * rho0 * f3 );
    }
  }
}


void PLocAssem_Tet_CMM_GenAlpha::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &dot_sol,
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  double f1, f2, f3;

  double tau_m, tau_c, tau_dc;

  const double two_mu = 2.0 * vis_mu;

  const double rho0_2 = rho0 * rho0;

  const double curr = time + alpha_f * dt;

  const double dd_dv = alpha_f * gamma * dt;

  Zero_Tangent_Residual();

  for(int qua=0; qua<nqp; ++qua)
  {
    double u = 0.0, u_t = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v = 0.0, v_t = 0.0, v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w = 0.0, w_t = 0.0, w_x = 0.0, w_y = 0.0, w_z = 0.0;
    double p = 0.0, coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;
    double p_x = 0.0, p_y = 0.0, p_z = 0.0;
    double u_xx = 0.0, u_yy = 0.0, u_zz = 0.0;
    double v_xx = 0.0, v_yy = 0.0, v_zz = 0.0;
    double w_xx = 0.0, w_yy = 0.0, w_zz = 0.0;

    element->get_3D_R_gradR_LaplacianR( qua, &R[0], &dR_dx[0], 
        &dR_dy[0], &dR_dz[0], &d2R_dxx[0], &d2R_dyy[0], &d2R_dzz[0] );

    element->get_invJacobian( qua, dxi_dx );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int ii4 = 4 * ii;

      u_t += dot_sol[ii4+1] * R[ii];
      v_t += dot_sol[ii4+2] * R[ii];
      w_t += dot_sol[ii4+3] * R[ii];

      u += sol[ii4+1] * R[ii];
      v += sol[ii4+2] * R[ii];
      w += sol[ii4+3] * R[ii];
      p += sol[ii4+0] * R[ii];

      u_x += sol[ii4+1] * dR_dx[ii];
      v_x += sol[ii4+2] * dR_dx[ii];
      w_x += sol[ii4+3] * dR_dx[ii];
      p_x += sol[ii4+0] * dR_dx[ii];

      u_y += sol[ii4+1] * dR_dy[ii];
      v_y += sol[ii4+2] * dR_dy[ii];
      w_y += sol[ii4+3] * dR_dy[ii];
      p_y += sol[ii4+0] * dR_dy[ii];

      u_z += sol[ii4+1] * dR_dz[ii];
      v_z += sol[ii4+2] * dR_dz[ii];
      w_z += sol[ii4+3] * dR_dz[ii];
      p_z += sol[ii4+0] * dR_dz[ii];

      u_xx += sol[ii4+1] * d2R_dxx[ii];
      u_yy += sol[ii4+1] * d2R_dyy[ii];
      u_zz += sol[ii4+1] * d2R_dzz[ii];

      v_xx += sol[ii4+2] * d2R_dxx[ii];
      v_yy += sol[ii4+2] * d2R_dyy[ii];
      v_zz += sol[ii4+2] * d2R_dzz[ii];

      w_xx += sol[ii4+3] * d2R_dxx[ii];
      w_yy += sol[ii4+3] * d2R_dyy[ii];
      w_zz += sol[ii4+3] * d2R_dzz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    get_tau(tau_m, tau_c, dt, dxi_dx, u, v, w);

    const double tau_m_2 = tau_m * tau_m;

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua); 

    get_f(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    const double u_lap = u_xx + u_yy + u_zz;
    const double v_lap = v_xx + v_yy + v_zz;
    const double w_lap = w_xx + w_yy + w_zz;

    const double rx = rho0 * ( u_t + u_x * u + u_y * v + u_z * w - f1 ) + p_x - vis_mu * u_lap;
    const double ry = rho0 * ( v_t + v_x * u + v_y * v + v_z * w - f2 ) + p_y - vis_mu * v_lap ;
    const double rz = rho0 * ( w_t + w_x * u + w_y * v + w_z * w - f3 ) + p_z - vis_mu * w_lap;

    const double div_vel = u_x + v_y + w_z;

    const double u_prime = -1.0 * tau_m * rx;
    const double v_prime = -1.0 * tau_m * ry;
    const double w_prime = -1.0 * tau_m * rz;

    const double r_dot_gradu = u_x * rx + u_y * ry + u_z * rz;
    const double r_dot_gradv = v_x * rx + v_y * ry + v_z * rz;
    const double r_dot_gradw = w_x * rx + w_y * ry + w_z * rz;

    get_DC( tau_dc, dxi_dx, u_prime, v_prime, w_prime );

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      const double velo_dot_gradR = NA_x * u + NA_y * v + NA_z * w;
      const double r_dot_gradR = NA_x * rx + NA_y * ry + NA_z * rz;
      const double velo_prime_dot_gradR = NA_x * u_prime + NA_y * v_prime + NA_z * w_prime;

      Residual[4*A] += gwts * ( NA * div_vel + tau_m * r_dot_gradR );

      Residual[4*A+1] += gwts * ( NA * rho0 * u_t
          + NA * rho0 * (u * u_x + v * u_y + w * u_z)
          - NA_x * p
          + NA_x * two_mu * u_x
          + NA_y * vis_mu * (u_y + v_x)
          + NA_z * vis_mu * (u_z + w_x)
          + velo_dot_gradR * tau_m * rho0 * rx
          - NA * tau_m * rho0 * r_dot_gradu
          + NA_x * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * rho0 * rx
          + velo_prime_dot_gradR * tau_dc 
          * (u_prime * u_x + v_prime * u_y + w_prime * u_z)
          - NA * rho0 * f1 );

      Residual[4*A+2] += gwts * ( NA * rho0 * v_t
          + NA * rho0 * (u * v_x + v * v_y + w * v_z)
          - NA_y * p
          + NA_x * vis_mu * (u_y + v_x)
          + NA_y * two_mu * v_y
          + NA_z * vis_mu * (v_z + w_y)
          + velo_dot_gradR * tau_m * rho0 * ry
          - NA * tau_m * rho0 * r_dot_gradv
          + NA_y * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * rho0 * ry
          + velo_prime_dot_gradR * tau_dc
          * (u_prime * v_x + v_prime * v_y + w_prime * v_z)
          - NA * rho0 * f2 );

      Residual[4*A+3] += gwts * (NA * rho0 * w_t
          + NA * rho0 * (u * w_x + v * w_y + w * w_z)
          - NA_z * p
          + NA_x * vis_mu * (u_z + w_x)
          + NA_y * vis_mu * (w_y + v_z)
          + NA_z * two_mu * w_z
          + velo_dot_gradR * tau_m * rho0 * rz
          - NA * tau_m * rho0 * r_dot_gradw
          + NA_z * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * rho0 * rz
          + velo_prime_dot_gradR * tau_dc
          * (u_prime * w_x + v_prime * w_y + w_prime * w_z)
          - NA * rho0 * f3 );

      for(int B=0; B<nLocBas; ++B)
      {
        const double NB = R[B], NB_x = dR_dx[B], NB_y = dR_dy[B], NB_z = dR_dz[B];
        const double NB_xx = d2R_dxx[B], NB_yy = d2R_dyy[B], NB_zz = d2R_dzz[B];
        const double NB_lap = NB_xx + NB_yy + NB_zz;
        const double velo_dot_gradNB = u * NB_x + v * NB_y + w * NB_z;
        const double velo_prime_dot_gradNB = u_prime * NB_x + v_prime * NB_y + w_prime * NB_z;

        const double NANB  = NA*NB, NANBx = NA*NB_x, NANBy = NA*NB_y, NANBz = NA*NB_z;
        const double NAxNB = NA_x*NB, NAxNBx = NA_x*NB_x, NAxNBy = NA_x*NB_y, NAxNBz = NA_x*NB_z;
        const double NAyNB = NA_y*NB, NAyNBx = NA_y*NB_x, NAyNBy = NA_y*NB_y, NAyNBz = NA_y*NB_z;
        const double NAzNB = NA_z*NB, NAzNBx = NA_z*NB_x, NAzNBy = NA_z*NB_y, NAzNBz = NA_z*NB_z;

        const double drx_du_B = rho0 * ( u_x * NB + velo_dot_gradNB ) - vis_mu * NB_lap; 
        const double drx_dv_B = rho0 * u_y * NB;
        const double drx_dw_B = rho0 * u_z * NB;

        const double dry_du_B = rho0 * v_x * NB;
        const double dry_dv_B = rho0 * ( v_y * NB + velo_dot_gradNB ) - vis_mu * NB_lap;
        const double dry_dw_B = rho0 * v_z * NB;

        const double drz_du_B = rho0 * w_x * NB;
        const double drz_dv_B = rho0 * w_y * NB;
        const double drz_dw_B = rho0 * ( w_z * NB + velo_dot_gradNB ) - vis_mu * NB_lap;

        // Continuity equation with respect to p, u, v, w
        Tangent[16*nLocBas*A+4*B] += gwts * dd_dv * tau_m * (NAxNBx + NAyNBy + NAzNBz);

        Tangent[16*nLocBas*A+4*B+1] += gwts * ( alpha_m * tau_m * rho0 * NAxNB
            + dd_dv * ( NANBx + tau_m * NA_x * drx_du_B
              + tau_m * NA_y * dry_du_B + tau_m * NA_z * drz_du_B ) );

        Tangent[16*nLocBas*A+4*B+2] += gwts * ( alpha_m * tau_m * rho0 * NAyNB
            + dd_dv * ( NANBy + tau_m * NA_x * drx_dv_B
              + tau_m * NA_y * dry_dv_B + tau_m * NA_z * drz_dv_B ) );

        Tangent[16*nLocBas*A+4*B+3] += gwts * ( alpha_m * tau_m * rho0 * NAzNB
            + dd_dv * ( NANBz + tau_m * NA_x * drx_dw_B
              + tau_m * NA_y * dry_dw_B + tau_m * NA_z * drz_dw_B ) );

        // Momentum-x with respect to p, u, v, w
        Tangent[4*nLocBas*(4*A+1)+4*B] += gwts * dd_dv * ((-1.0) * NAxNB
            + velo_dot_gradR * tau_m * rho0 * NB_x
            - NA * tau_m * rho0 * (u_x * NB_x + u_y * NB_y + u_z * NB_z)
            - 2.0 * tau_m_2 * rho0 * rx * NAxNBx
            - tau_m_2 * rho0 * NA_y * (rx * NB_y + ry * NB_x)
            - tau_m_2 * rho0 * NA_z * (rx * NB_z + rz * NB_x) );

        Tangent[4*nLocBas*(4*A+1)+4*B+1] += gwts * ( 
            alpha_m * ( rho0 * NANB + velo_dot_gradR * rho0_2 * tau_m * NB
              - rho0_2 * tau_m * u_x * NANB
              - rho0_2 * tau_m_2 * rx * NAxNB
              - rho0_2 * tau_m_2 * (rx * NAxNB + ry * NAyNB + rz * NAzNB) )
            + dd_dv * ( NA * rho0 * velo_dot_gradNB + NANB * rho0 * u_x
              + vis_mu * (2.0*NAxNBx + NAyNBy + NAzNBz)
              + velo_dot_gradR * rho0 * tau_m * drx_du_B
              + rho0 * tau_m * rx * NAxNB
              - rho0 * tau_m * (rx * NANBx + ry * NANBy + rz * NANBz)
              - rho0 * tau_m * NA * (u_x * drx_du_B 
                + u_y * dry_du_B + u_z * drz_du_B ) 
              + tau_c * NAxNBx
              - 2.0 * rho0 * tau_m_2 * rx  * NA_x * drx_du_B
              - rho0 * tau_m_2 * ry * NA_y * drx_du_B
              - rho0 * tau_m_2 * rz * NA_z * drx_du_B
              - rho0 * tau_m_2 * rx * NA_y * dry_du_B
              - rho0 * tau_m_2 * rx * NA_z * drz_du_B
              + velo_prime_dot_gradR * tau_dc * velo_prime_dot_gradNB ) );

        Tangent[4*nLocBas*(4*A+1)+4*B+2] += gwts * ( 
            alpha_m * (-1.0) * rho0_2 * (tau_m * u_y * NANB + tau_m_2 * rx * NAyNB)
            + dd_dv * ( NANB * rho0 * u_y + vis_mu * NAyNBx 
              + rho0 * tau_m * rx * NAyNB
              + velo_dot_gradR * rho0 * tau_m * drx_dv_B
              - rho0 * tau_m * NA * (u_x*drx_dv_B + u_y*dry_dv_B + u_z*drz_dv_B)
              + tau_c * NAxNBy
              - 2.0 * rho0 * tau_m_2 * rx * NA_x * drx_dv_B
              - rho0 * tau_m_2 * NA_y * (rx * dry_dv_B + ry * drx_dv_B)
              - rho0 * tau_m_2 * NA_z * (rx * drz_dv_B + rz * drx_dv_B) ) );

        Tangent[4*nLocBas*(4*A+1)+4*B+3] += gwts * (
            alpha_m * (-1.0) * rho0_2 * (tau_m * u_z * NANB + tau_m_2 * rx * NAzNB)
            + dd_dv * ( NANB * rho0 * u_z + vis_mu * NAzNBx 
              + rho0 * tau_m * rx * NAzNB
              + velo_dot_gradR * rho0 * tau_m * drx_dw_B
              - rho0 * tau_m * NA * (u_x*drx_dw_B + u_y*dry_dw_B + u_z*drz_dw_B)
              + tau_c * NAxNBz
              - 2.0 * rho0 * tau_m_2 * rx * NA_x * drx_dw_B
              - rho0 * tau_m_2 * NA_y * (rx * dry_dw_B + ry * drx_dw_B)
              - rho0 * tau_m_2 * NA_z * (rx * drz_dw_B + rz * drx_dw_B) ) );

        // Momentum-y with respect to p u v w
        Tangent[4*nLocBas*(4*A+2)+4*B] += gwts * dd_dv * ( (-1.0) * NAyNB
            + velo_dot_gradR * tau_m * rho0 * NB_y
            - NA * tau_m * rho0 * (v_x * NB_x + v_y * NB_y + v_z * NB_z)
            - tau_m_2 * rho0 * NA_x * (rx * NB_y + ry * NB_x)
            - 2.0 * tau_m_2 * rho0 * ry * NAyNBy
            - tau_m_2 * rho0 * NA_z * (ry * NB_z + rz * NB_y) );

        Tangent[4*nLocBas*(4*A+2)+4*B+1] += gwts * (
            alpha_m * (-1.0) * rho0_2 * (tau_m * v_x * NANB + tau_m_2 * ry * NAxNB)
            + dd_dv * ( NANB * rho0 * v_x + vis_mu * NAxNBy
              + rho0 * tau_m * ry * NAxNB
              + velo_dot_gradR * rho0 * tau_m * dry_du_B
              - rho0 * tau_m * NA * (v_x*drx_du_B + v_y*dry_du_B + v_z*drz_du_B)
              + tau_c * NAyNBx
              - rho0 * tau_m_2 * NA_x * (ry * drx_du_B + rx * dry_du_B)
              - 2.0 * rho0 * tau_m_2 * ry * NA_y * dry_du_B
              - rho0 * tau_m_2 * NA_z * (ry * drz_du_B + rz * dry_du_B) ) );

        Tangent[4*nLocBas*(4*A+2)+4*B+2] += gwts * (
            alpha_m * ( rho0 * NANB + velo_dot_gradR * rho0_2 * tau_m * NB
              - rho0_2 * tau_m * v_y * NANB
              - rho0_2 * tau_m_2 * ry * NAyNB
              - rho0_2 * tau_m_2 * (rx * NAxNB + ry * NAyNB + rz * NAzNB) )
            + dd_dv * ( NA * rho0 * velo_dot_gradNB + NANB * rho0 * v_y
              + vis_mu * (NAxNBx + 2.0 * NAyNBy + NAzNBz)
              + velo_dot_gradR * rho0 * tau_m * dry_dv_B
              + rho0 * tau_m * ry * NAyNB
              - rho0 * tau_m * ( rx * NANBx + ry * NANBy + rz * NANBz )
              - rho0 * tau_m * NA * (v_x * drx_dv_B + v_y * dry_dv_B + v_z * drz_dv_B)
              + tau_c * NAyNBy
              - rho0 * tau_m_2 * NA_x * (rx * dry_dv_B + ry * drx_dv_B)
              - 2.0 * rho0 * tau_m_2 * ry * NA_y * dry_dv_B
              - rho0 * tau_m_2 * NA_z * (ry * drz_dv_B + rz * dry_dv_B)
              + velo_prime_dot_gradR * tau_dc * velo_prime_dot_gradNB ) );

        Tangent[4*nLocBas*(4*A+2)+4*B+3] += gwts * (
            alpha_m * (-1.0) * rho0_2 * ( tau_m * v_z * NANB + tau_m_2 * ry * NAzNB ) 
            + dd_dv * ( NANB * rho0 * v_z + vis_mu * NAzNBy
              + rho0 * tau_m * ry * NAzNB
              + velo_dot_gradR * rho0 * tau_m * dry_dw_B
              - rho0 * tau_m * NA * (v_x*drx_dw_B + v_y*dry_dw_B + v_z*drz_dw_B)
              + tau_c * NAyNBz
              - rho0 * tau_m_2 * NA_x * (rx * dry_dw_B + ry * drx_dw_B)
              - rho0 * tau_m_2 * 2.0 * ry * NA_y * dry_dw_B
              - rho0 * tau_m_2 * NA_z * (ry * drz_dw_B + rz * dry_dw_B) ) );

        // Momentum-z with respect to p u v w
        Tangent[4*nLocBas*(4*A+3)+4*B] += gwts * dd_dv * ( (-1.0) * NAzNB
            + velo_dot_gradR * tau_m * rho0 * NB_z
            - NA * tau_m * rho0 * (w_x * NB_x + w_y * NB_y + w_z * NB_z)
            - tau_m_2 * rho0 * NA_x * (rx * NB_z + rz * NB_x)
            - tau_m_2 * rho0 * NA_y * (ry * NB_z + rz * NB_y)
            - 2.0 * tau_m_2 * rho0 * rz * NAzNBz );

        Tangent[4*nLocBas*(4*A+3)+4*B+1] += gwts * (
            alpha_m * (-1.0) * rho0_2 * (tau_m * w_x * NANB + tau_m_2 * rz * NAxNB)
            + dd_dv * ( NANB * rho0 * w_x + vis_mu * NAxNBz
              + rho0 * tau_m * rz * NAxNB
              + velo_dot_gradR * rho0 * tau_m * drz_du_B
              - rho0 * tau_m * NA * (w_x*drx_du_B + w_y*dry_du_B + w_z*drz_du_B)
              + tau_c * NAzNBx
              - rho0 * tau_m_2 * NA_x * (rx * drz_du_B + rz * drx_du_B)
              - rho0 * tau_m_2 * NA_y * (ry * drz_du_B + rz * dry_du_B)
              - 2.0 * rho0 * tau_m_2 * rz * NA_z * drz_du_B ) );

        Tangent[4*nLocBas*(4*A+3)+4*B+2] += gwts * (
            alpha_m * (-1.0) * rho0_2 * (tau_m * w_y * NANB + tau_m_2 * rz * NAyNB)
            + dd_dv * ( NANB * rho0 * w_y + vis_mu * NAyNBz
              + rho0 * tau_m * rz * NAyNB
              + velo_dot_gradR * rho0 * tau_m * drz_dv_B
              - rho0 * tau_m * NA * (w_x*drx_dv_B + w_y*dry_dv_B + w_z*drz_dv_B)
              + tau_c * NAzNBy
              - rho0 * tau_m_2 * NA_x * (rx * drz_dv_B + rz * drx_dv_B)
              - rho0 * tau_m_2 * NA_y * (ry * drz_dv_B + rz * dry_dv_B)
              - 2.0 * rho0 * tau_m_2 * rz * NA_z * drz_dv_B ) );

        Tangent[4*nLocBas*(4*A+3)+4*B+3] += gwts * (
            alpha_m * ( rho0 * NANB + velo_dot_gradR * rho0_2 * tau_m * NB
              - rho0_2 * tau_m * w_z * NANB
              - rho0_2 * tau_m_2 * rz * NAzNB
              - rho0_2 * tau_m_2 * (rx*NAxNB + ry*NAyNB + rz * NAzNB) )
            + dd_dv * ( rho0 * NA * velo_dot_gradNB + NANB * rho0 * w_z
              + vis_mu * (NAxNBx + NAyNBy + 2.0 * NAzNBz)
              + velo_dot_gradR * rho0 * tau_m * drz_dw_B
              + rho0 * tau_m * rz * NAzNB
              - rho0 * tau_m * (rx*NANBx + ry * NANBy + rz * NANBz)
              - rho0 * tau_m * NA * (w_x*drx_dw_B + w_y*dry_dw_B + w_z*drz_dw_B)
              + tau_c * NAzNBz
              - rho0 * tau_m_2 * NA_x * (rx * drz_dw_B + rz * drx_dw_B)
              - rho0 * tau_m_2 * NA_y * (ry * drz_dw_B + rz * dry_dw_B)
              - 2.0 * rho0 * tau_m_2 * NA_z * rz * drz_dw_B 
              + velo_prime_dot_gradR * tau_dc * velo_prime_dot_gradNB ) );
      } // B-loop
    } // A-loop
  } // qua-loop
  // ----------------------------------------------------------------
  // The local `stiffness' matrix 
  //            K[p][q] = Sub_Tan[4*ii + jj][A*nLocBas+B],
  // where p = 4*A+ii, q = 4*B+jj, and K has 4*nLocBas rows/columns.
  // Tangent is a 1D vector storing K by rows:
  // Tangent[4*nLocBas*p + q] = K[p][q] = Sub_Tan[4*ii+jj][A*nLocBas+B]
  // ----------------------------------------------------------------
}


void PLocAssem_Tet_CMM_GenAlpha::Assem_Mass_Residual(
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  double f1, f2, f3;

  const double two_mu = 2.0 * vis_mu;

  const double curr = 0.0;

  Zero_Tangent_Residual();

  for(int qua=0; qua<nqp; ++qua)
  {
    double u = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v = 0.0, v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w = 0.0, w_x = 0.0, w_y = 0.0, w_z = 0.0;
    double p = 0.0, coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    element->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int ii4 = ii * 4;

      u += sol[ii4+1] * R[ii];
      v += sol[ii4+2] * R[ii];
      w += sol[ii4+3] * R[ii];
      p += sol[ii4+0] * R[ii];

      u_x += sol[ii4+1] * dR_dx[ii];
      v_x += sol[ii4+2] * dR_dx[ii];
      w_x += sol[ii4+3] * dR_dx[ii];

      u_y += sol[ii4+1] * dR_dy[ii];
      v_y += sol[ii4+2] * dR_dy[ii];
      w_y += sol[ii4+3] * dR_dy[ii];

      u_z += sol[ii4+1] * dR_dz[ii];
      v_z += sol[ii4+2] * dR_dz[ii];
      w_z += sol[ii4+3] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    get_f(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      Residual[4*A+1] += gwts * ( NA * rho0 * (u*u_x + v*u_y + w*u_z) 
          - NA_x * p
          + two_mu * NA_x * u_x
          + vis_mu * NA_y * (u_y + v_x)
          + vis_mu * NA_z * (u_z + w_x)
          - NA * rho0 * f1 );

      Residual[4*A+2] += gwts * ( NA * rho0 * (u*v_x + v*v_y + w*v_z) 
          - NA_y * p
          + vis_mu * NA_x * (u_y + v_x)
          + two_mu * NA_y * v_y
          + vis_mu * NA_z * (v_z + w_y)
          - NA * rho0 * f2 );

      Residual[4*A+3] += gwts * ( NA * rho0 * (u*w_x + v*w_y + w*w_z) 
          - NA_z * p
          + vis_mu * NA_x * (u_z + w_x)
          + vis_mu * NA_y * (w_y + v_z)
          + two_mu * NA_z * w_z
          - NA * rho0 * f3 );

      for(int B=0; B<nLocBas; ++B)
      {
        Tangent[4*nLocBas*(4*A) + 4*B] += gwts * rho0 * NA * R[B];
        Tangent[4*nLocBas*(4*A+1) + 4*B+1] += gwts * rho0 * NA * R[B];
        Tangent[4*nLocBas*(4*A+2) + 4*B+2] += gwts * rho0 * NA * R[B];
        Tangent[4*nLocBas*(4*A+3) + 4*B+3] += gwts * rho0 * NA * R[B];
      }
    }
  }
}


void PLocAssem_Tet_CMM_GenAlpha::Assem_Residual_EBC(
    const int &ebc_id,
    const double &time, const double &dt,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  const double curr = time + alpha_f * dt;

  double gx, gy, gz, nx, ny, nz, surface_area;

  Zero_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, &R[0]);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);

    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    get_ebc_fun( ebc_id, coor_x, coor_y, coor_z, curr, nx, ny, nz,
        gx, gy, gz );

    for(int A=0; A<snLocBas; ++A)
    {
      Residual[4*A+1] -= surface_area * quad -> get_qw(qua) * R[A] * gx;
      Residual[4*A+2] -= surface_area * quad -> get_qw(qua) * R[A] * gy;
      Residual[4*A+3] -= surface_area * quad -> get_qw(qua) * R[A] * gz;
    }
  }
}


void PLocAssem_Tet_CMM_GenAlpha::Assem_Residual_EBC_Resistance(
    const int &ebc_id,
    const double &val,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  double nx, ny, nz, surface_area;

  Zero_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, &R[0]);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);

    for(int A=0; A<snLocBas; ++A)
    {
      Residual[4*A+1] += surface_area * quad -> get_qw(qua) * R[A] * nx * val;
      Residual[4*A+2] += surface_area * quad -> get_qw(qua) * R[A] * ny * val;
      Residual[4*A+3] += surface_area * quad -> get_qw(qua) * R[A] * nz * val;
    }
  }
}


void PLocAssem_Tet_CMM_GenAlpha::Assem_Residual_BackFlowStab(
    const double * const &dot_sol,
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  double nx, ny, nz, surface_area, factor;

  Zero_sur_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, &R[0]);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);

    double u = 0.0, v = 0.0, w = 0.0;;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      const int ii4 = ii * 4;
      u += sol[ii4+1] * R[ii];
      v += sol[ii4+2] * R[ii];
      w += sol[ii4+3] * R[ii];
    }

    const double temp = u * nx + v * ny + w * nz;

    if(temp < 0.0) factor = temp * rho0 * beta;
    else factor = 0.0;

    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual[4*A+1] -= surface_area * quad -> get_qw(qua) * R[A] * factor * u;
      sur_Residual[4*A+2] -= surface_area * quad -> get_qw(qua) * R[A] * factor * v;
      sur_Residual[4*A+3] -= surface_area * quad -> get_qw(qua) * R[A] * factor * w;
    }
  }
}


void PLocAssem_Tet_CMM_GenAlpha::Assem_Tangent_Residual_BackFlowStab(
    const double &dt,
    const double * const &dot_sol,
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  const double dd_dv = alpha_f * gamma * dt;

  double nx, ny, nz, surface_area, factor;

  Zero_sur_Tangent_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, &R[0]);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);

    double u = 0.0, v = 0.0, w = 0.0;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      u += sol[ii*4+1] * R[ii];
      v += sol[ii*4+2] * R[ii];
      w += sol[ii*4+3] * R[ii];
    }

    const double temp = u * nx + v * ny + w * nz;

    if(temp < 0.0) factor = temp * rho0 * beta;
    else factor = 0.0;

    const double gwts = surface_area * quad -> get_qw(qua);

    // snLocBas = 3 for linear tri element
    //            6 for quadratic tri element
    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual[4*A+1] -= gwts * R[A] * factor * u;
      sur_Residual[4*A+2] -= gwts * R[A] * factor * v;
      sur_Residual[4*A+3] -= gwts * R[A] * factor * w;

      for(int B=0; B<snLocBas; ++B)
      {
        // index := A *snLocBas+B here ranges from 0 to 8 for linear triangle
        //                        0 to 35 for quadratic triangle
        sur_Tangent[ 4*snLocBas*(4*A+1) + 4*B+1 ] -= gwts * dd_dv * R[A] * factor * R[B];
        sur_Tangent[ 4*snLocBas*(4*A+2) + 4*B+2 ] -= gwts * dd_dv * R[A] * factor * R[B];
        sur_Tangent[ 4*snLocBas*(4*A+3) + 4*B+3 ] -= gwts * dd_dv * R[A] * factor * R[B];
      }
    }
  }
}


double PLocAssem_Tet_CMM_GenAlpha::get_flowrate( const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  double nx, ny, nz, surface_area;

  double flrate = 0.0;

  for(int qua =0; qua< face_nqp; ++qua)
  {
    element->get_R(qua, &R[0]);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);

    double u = 0.0, v = 0.0, w = 0.0;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      u += sol[ii*4+1] * R[ii];
      v += sol[ii*4+2] * R[ii];
      w += sol[ii*4+3] * R[ii];
    }

    flrate += surface_area * quad->get_qw(qua) * ( u * nx + v * ny + w * nz );
  }

  return flrate;
}


void PLocAssem_Tet_CMM_GenAlpha::get_pressure_area( 
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad,
    double &pres, double &area )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  double nx, ny, nz, surface_area;

  // Initialize the two variables to be passed out
  pres = 0.0; area = 0.0;

  for(int qua =0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, &R[0]);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);

    double pp = 0.0;
    for(int ii=0; ii<snLocBas; ++ii) pp += sol[4*ii+0] * R[ii];

    pres += surface_area * quad->get_qw(qua) * pp;
    area += surface_area * quad->get_qw(qua);
  }
}


void PLocAssem_Tet_CMM_GenAlpha::Assem_Residual_EBC_Wall(
    const double &time, const double &dt,
    const double * const &dot_sol,
    const double * const &sol_wall_disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &ele_thickness,
    const double * const &ele_youngsmod,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z ); 

  const int dim = 3;

  const int face_nqp = quad -> get_num_quadPts();
  const double curr = time + alpha_f * dt;

  Zero_sur_Residual();

  for(int qua=0; qua<face_nqp; ++qua)
  {
    // Lamina displacements
    double * sol_wall_disp_l = new double [ snLocBas * dim ];

    // For membrane elements, basis function gradients are computed
    // with respect to lamina coords
    double * dR_dxl = new double [ snLocBas ];
    double * dR_dyl = new double [ snLocBas ];
    element->get_R_gradR( qua, &R[0], &dR_dxl[0], &dR_dyl[0] );

    // Global-to-local rotation matrix Q
    Matrix_3x3 Q = Matrix_3x3();
    element->get_rotationMatrix(qua, Q);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      sol_wall_disp_l[dim*ii]   = sol_wall_disp[dim*ii] * Q(0, 0) 
        + sol_wall_disp[dim*ii+1] * Q(0, 1) + sol_wall_disp[dim*ii+2] * Q(0, 2);

      sol_wall_disp_l[dim*ii+1] = sol_wall_disp[dim*ii] * Q(1, 0) 
        + sol_wall_disp[dim*ii+1] * Q(1, 1) + sol_wall_disp[dim*ii+2] * Q(1, 2);

      sol_wall_disp_l[dim*ii+2] = sol_wall_disp[dim*ii] * Q(2, 0) 
        + sol_wall_disp[dim*ii+1] * Q(2, 1) + sol_wall_disp[dim*ii+2] * Q(2, 2);
    }

    double u_t = 0.0, v_t = 0.0, w_t = 0.0;
    double h_w = 0.0, E_w = 0.0;
    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;
    double u1l_xl = 0.0, u2l_xl = 0.0, u3l_xl = 0.0;
    double u1l_yl = 0.0, u2l_yl = 0.0, u3l_yl = 0.0;

    for(int ii=0; ii<snLocBas; ++ii)
    {
      const int ii4 = 4 * ii;

      u_t += dot_sol[ii4+1] * R[ii];
      v_t += dot_sol[ii4+2] * R[ii];
      w_t += dot_sol[ii4+3] * R[ii];

      h_w += ele_thickness[ii] * R[ii];
      E_w += ele_youngsmod[ii] * R[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];

      u1l_xl += sol_wall_disp_l[dim*ii]   * dR_dxl[ii];
      u1l_yl += sol_wall_disp_l[dim*ii]   * dR_dyl[ii];
      u2l_xl += sol_wall_disp_l[dim*ii+1] * dR_dxl[ii];
      u2l_yl += sol_wall_disp_l[dim*ii+1] * dR_dyl[ii];
      u3l_xl += sol_wall_disp_l[dim*ii+2] * dR_dxl[ii];
      u3l_yl += sol_wall_disp_l[dim*ii+2] * dR_dyl[ii];
    }

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    // Body force acting on the wall
    double f1, f2, f3;
    get_fw(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    const double coef = E_w / (1.0 - nu_w * nu_w);

    // Lamina Cauchy stress
    Matrix_3x3 sigma = Matrix_3x3(
      u1l_xl + nu_w*u2l_yl,               0.5*(1.0-nu_w) * (u1l_yl + u2l_xl), 0.5*kappa_w*(1.0-nu_w) * u3l_xl,
      0.5*(1.0-nu_w) * (u1l_yl + u2l_xl), nu_w*u1l_xl + u2l_yl,               0.5*kappa_w*(1.0-nu_w) * u3l_yl,
      0.5*kappa_w*(1.0-nu_w) * u3l_xl,    0.5*kappa_w*(1.0-nu_w) * u3l_yl,    0.0 );
    sigma *= coef;

    // Global Cauchy stress: Q^T * lamina_Cauchy * Q
    sigma.MatRot(Q);

    // Basis function gradients with respect to global coords
    // dR/dx_{i} = Q_{ji} * dR/dxl_{j}. Note that dR/dzl = 0.0
    for(int ii=0; ii<snLocBas; ++ii)
    {
      dR_dx[ii] = Q(0, 0) * dR_dxl[ii] + Q(1, 0) * dR_dyl[ii];
      dR_dy[ii] = Q(0, 1) * dR_dxl[ii] + Q(1, 1) * dR_dyl[ii];
      dR_dz[ii] = Q(0, 2) * dR_dxl[ii] + Q(1, 2) * dR_dyl[ii];
    }

    for(int A=0; A<snLocBas; ++A)
    {
      const double NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      sur_Residual[4*A+1] += gwts * h_w * ( R[A] * rho_w * ( u_t - f1 )
          + NA_x * sigma(0, 0) + NA_y * sigma(0, 1) + NA_z * sigma(0, 2) ); 
      sur_Residual[4*A+2] += gwts * h_w * ( R[A] * rho_w * ( v_t - f2 )
          + NA_x * sigma(1, 0) + NA_y * sigma(1, 1) + NA_z * sigma(1, 2) ); 
      sur_Residual[4*A+3] += gwts * h_w * ( R[A] * rho_w * ( w_t - f3 )
          + NA_x * sigma(2, 0) + NA_y * sigma(2, 1) + NA_z * sigma(2, 2) ); 
    }

    delete [] sol_wall_disp_l; delete [] dR_dxl; delete [] dR_dyl;
    sol_wall_disp_l = nullptr; dR_dxl = nullptr; dR_dyl = nullptr;

  } // end qua loop
}


void PLocAssem_Tet_CMM_GenAlpha::Assem_Tangent_Residual_EBC_Wall(
    const double &time, const double &dt,
    const double * const &dot_sol,
    const double * const &sol_wall_disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &ele_thickness,
    const double * const &ele_youngsmod,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z ); 

  const int dim = 3;

  const int face_nqp = quad -> get_num_quadPts();

  const double dd_dv = alpha_f * gamma * dt;
  const double dd_du = dd_dv * dd_dv / alpha_m;

  const double curr = time + alpha_f * dt;

  Zero_sur_Tangent_Residual();

  for(int qua=0; qua<face_nqp; ++qua)
  {
    // Lamina displacements
    double * sol_wall_disp_l = new double [ snLocBas * dim ];

    // For membrane elements, basis function gradients are computed
    // with respect to lamina coords
    double * dR_dxl = new double [ snLocBas ];
    double * dR_dyl = new double [ snLocBas ];
    element->get_R_gradR( qua, &R[0], &dR_dxl[0], &dR_dyl[0] );

    // Lamina and global stiffness matrices
    double * Kl = new double [ (snLocBas*dim) * (snLocBas*dim) ] {};
    double * Kg = new double [ (snLocBas*dim) * (snLocBas*dim) ] {};

    // Global-to-local rotation matrix Q
    Matrix_3x3 Q = Matrix_3x3();
    element->get_rotationMatrix(qua, Q);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      sol_wall_disp_l[dim*ii]   = sol_wall_disp[dim*ii] * Q(0, 0) 
        + sol_wall_disp[dim*ii+1] * Q(0, 1) + sol_wall_disp[dim*ii+2] * Q(0, 2);

      sol_wall_disp_l[dim*ii+1] = sol_wall_disp[dim*ii] * Q(1, 0) 
        + sol_wall_disp[dim*ii+1] * Q(1, 1) + sol_wall_disp[dim*ii+2] * Q(1, 2);

      sol_wall_disp_l[dim*ii+2] = sol_wall_disp[dim*ii] * Q(2, 0) 
        + sol_wall_disp[dim*ii+1] * Q(2, 1) + sol_wall_disp[dim*ii+2] * Q(2, 2);
    }

    double u_t = 0.0, v_t = 0.0, w_t = 0.0;
    double h_w = 0.0, E_w = 0.0;
    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;
    double u1l_xl = 0.0, u2l_xl = 0.0, u3l_xl = 0.0;
    double u1l_yl = 0.0, u2l_yl = 0.0, u3l_yl = 0.0;

    for(int ii=0; ii<snLocBas; ++ii)
    {
      const int ii4 = 4 * ii;

      u_t += dot_sol[ii4+1] * R[ii];
      v_t += dot_sol[ii4+2] * R[ii];
      w_t += dot_sol[ii4+3] * R[ii];

      h_w += ele_thickness[ii] * R[ii];
      E_w += ele_youngsmod[ii] * R[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];

      u1l_xl += sol_wall_disp_l[dim*ii]   * dR_dxl[ii];
      u1l_yl += sol_wall_disp_l[dim*ii]   * dR_dyl[ii];
      u2l_xl += sol_wall_disp_l[dim*ii+1] * dR_dxl[ii];
      u2l_yl += sol_wall_disp_l[dim*ii+1] * dR_dyl[ii];
      u3l_xl += sol_wall_disp_l[dim*ii+2] * dR_dxl[ii];
      u3l_yl += sol_wall_disp_l[dim*ii+2] * dR_dyl[ii];
    }

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    // Body force acting on the wall
    double f1, f2, f3;
    get_fw(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    const double coef = E_w / (1.0 - nu_w * nu_w);

    // Lamina Cauchy stress
    Matrix_3x3 sigma = Matrix_3x3(
      u1l_xl + nu_w*u2l_yl,               0.5*(1.0-nu_w) * (u1l_yl + u2l_xl), 0.5*kappa_w*(1.0-nu_w) * u3l_xl,
      0.5*(1.0-nu_w) * (u1l_yl + u2l_xl), nu_w*u1l_xl + u2l_yl,               0.5*kappa_w*(1.0-nu_w) * u3l_yl,
      0.5*kappa_w*(1.0-nu_w) * u3l_xl,    0.5*kappa_w*(1.0-nu_w) * u3l_yl,    0.0 );
    sigma *= coef;

    // Global Cauchy stress: Q^T * lamina_Cauchy * Q
    sigma.MatRot(Q);

    // Basis function gradients with respect to global coords
    // dR/dx_{i} = Q_{ji} * dR/dxl_{j}. Note that dR/dzl = 0.0
    for(int ii=0; ii<snLocBas; ++ii)
    {
      dR_dx[ii] = Q(0, 0) * dR_dxl[ii] + Q(1, 0) * dR_dyl[ii];
      dR_dy[ii] = Q(0, 1) * dR_dxl[ii] + Q(1, 1) * dR_dyl[ii];
      dR_dz[ii] = Q(0, 2) * dR_dxl[ii] + Q(1, 2) * dR_dyl[ii];
    }

    // Stiffness tensor in lamina coords
    // Bl^T * D * Bl = Bl_{ki} * D_{kl} * Bl_{lj}
    for(int A=0; A<snLocBas; ++A)
    {
      const double NA_xl = dR_dxl[A], NA_yl = dR_dyl[A];

      for(int B=0; B<snLocBas; ++B)
      {
        const double NB_xl = dR_dxl[B], NB_yl = dR_dyl[B];

        // Momentum-x with respect to u1, u2 
        Kl[(snLocBas*dim)*(A*dim) + (B*dim)]     += coef * ( NA_xl * NB_xl
            + 0.5*(1.0-nu_w) * NA_yl * NB_yl );
        Kl[(snLocBas*dim)*(A*dim) + (B*dim+1)]   += coef * ( nu_w * NA_xl * NB_yl
            + 0.5*(1.0-nu_w) * NA_yl * NB_xl );

        // Momentum-y with respect to u1, u2 
        Kl[(snLocBas*dim)*(A*dim+1) + (B*dim)]   += coef * ( nu_w * NA_yl * NB_xl
            + 0.5*(1.0-nu_w) * NA_xl * NB_yl );
        Kl[(snLocBas*dim)*(A*dim+1) + (B*dim+1)] += coef * ( NA_yl * NB_yl
            + 0.5*(1.0-nu_w) * NA_xl * NB_xl );

        // Momentum-z with respect to u3 
        Kl[(snLocBas*dim)*(A*dim+2) + (B*dim+2)] += coef * 0.5*kappa_w*(1.0-nu_w) * (
            NA_xl * NB_xl + NA_yl * NB_yl );
      }
    }

    // Stiffness tensor in global coords
    // theta^T * Kl * theta, where theta = [Q, 0, 0; 0, Q, 0; 0, 0, Q]
    // or Q^T * Kl_[AB] * Q = Q_{ki} * Kl_[AB]{kl} * Q_{lj}
    for(int A=0; A<snLocBas; ++A)
    {
      for(int B=0; B<snLocBas; ++B)
      {
        for(int ii=0; ii<dim; ++ii)
        {
          for(int jj=0; jj<dim; ++jj)
          {
            // Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(kk,ii) * Kl[ (A*dim+kk)*(snLocBas*dim) + (B*dim+ll) ] * Q(ll,jj);
            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(0,ii) * Kl[ (A*dim+0)*(snLocBas*dim) + (B*dim+0) ] * Q(0,jj);
            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(0,ii) * Kl[ (A*dim+0)*(snLocBas*dim) + (B*dim+1) ] * Q(1,jj);
            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(0,ii) * Kl[ (A*dim+0)*(snLocBas*dim) + (B*dim+2) ] * Q(2,jj);

            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(1,ii) * Kl[ (A*dim+1)*(snLocBas*dim) + (B*dim+0) ] * Q(0,jj);
            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(1,ii) * Kl[ (A*dim+1)*(snLocBas*dim) + (B*dim+1) ] * Q(1,jj);
            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(1,ii) * Kl[ (A*dim+1)*(snLocBas*dim) + (B*dim+2) ] * Q(2,jj);

            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(2,ii) * Kl[ (A*dim+2)*(snLocBas*dim) + (B*dim+0) ] * Q(0,jj);
            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(2,ii) * Kl[ (A*dim+2)*(snLocBas*dim) + (B*dim+1) ] * Q(1,jj);
            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(2,ii) * Kl[ (A*dim+2)*(snLocBas*dim) + (B*dim+2) ] * Q(2,jj);
          }
        }
      }
    }

    for(int A=0; A<snLocBas; ++A)
    {
      const double NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      sur_Residual[4*A+1] += gwts * h_w * ( R[A] * rho_w * ( u_t - f1 )
          + NA_x * sigma(0, 0) + NA_y * sigma(0, 1) + NA_z * sigma(0, 2) ); 
      sur_Residual[4*A+2] += gwts * h_w * ( R[A] * rho_w * ( v_t - f2 )
          + NA_x * sigma(1, 0) + NA_y * sigma(1, 1) + NA_z * sigma(1, 2) ); 
      sur_Residual[4*A+3] += gwts * h_w * ( R[A] * rho_w * ( w_t - f3 )
          + NA_x * sigma(2, 0) + NA_y * sigma(2, 1) + NA_z * sigma(2, 2) ); 

      for(int B=0; B<snLocBas; ++B)
      {
        // Momentum-x with respect to u, v, w
        sur_Tangent[ 4*snLocBas*(4*A+1) + 4*B+1 ] += gwts * h_w * (
            alpha_m * rho_w * R[A] * R[B]
            + dd_du * Kg[ (snLocBas*dim)*(A*dim) + (B*dim) ] );

        sur_Tangent[ 4*snLocBas*(4*A+1) + 4*B+2 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim) + (B*dim+1) ] );

        sur_Tangent[ 4*snLocBas*(4*A+1) + 4*B+3 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim) + (B*dim+2) ] );

        // Momentum-y with respect to u, v, w
        sur_Tangent[ 4*snLocBas*(4*A+2) + 4*B+1 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim+1) + (B*dim) ] );

        sur_Tangent[ 4*snLocBas*(4*A+2) + 4*B+2 ] += gwts * h_w * (
            alpha_m * rho_w * R[A] * R[B]
            + dd_du * Kg[ (snLocBas*dim)*(A*dim+1) + (B*dim+1) ] );

        sur_Tangent[ 4*snLocBas*(4*A+2) + 4*B+3 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim+1) + (B*dim+2) ] );

        // Momentum-z with respect to u, v, w
        sur_Tangent[ 4*snLocBas*(4*A+3) + 4*B+1 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim+2) + (B*dim) ] );

        sur_Tangent[ 4*snLocBas*(4*A+3) + 4*B+2 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim+2) + (B*dim+1) ] );

        sur_Tangent[ 4*snLocBas*(4*A+3) + 4*B+3 ] += gwts * h_w * (
            alpha_m * rho_w * R[A] * R[B]
            + dd_du * Kg[ (snLocBas*dim)*(A*dim+2) + (B*dim+2) ] );

      } // end B loop
    } // end A loop

    delete [] sol_wall_disp_l; delete [] dR_dxl; delete [] dR_dyl;
    delete [] Kl; delete [] Kg;
    sol_wall_disp_l = nullptr; dR_dxl = nullptr; dR_dyl = nullptr;
    Kl = nullptr; Kg = nullptr;

  } // end qua loop
}

// EOF
