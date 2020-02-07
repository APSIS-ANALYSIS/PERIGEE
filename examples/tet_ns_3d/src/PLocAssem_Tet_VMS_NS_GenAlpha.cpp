#include "PLocAssem_Tet_VMS_NS_GenAlpha.hpp"

PLocAssem_Tet_VMS_NS_GenAlpha::PLocAssem_Tet_VMS_NS_GenAlpha(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas,
        const double &in_rho, const double &in_vis_mu,
        const double &in_beta )
: rho0( in_rho ), vis_mu( in_vis_mu ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), beta(in_beta),
  dof_per_node(4), nqp(in_nqp)
{
  // for linear elements
  CI = 36.0; CT = 4.0;
  nLocBas = 4; snLocBas = 3;

  vec_size = nLocBas * dof_per_node;
  sur_size = snLocBas * dof_per_node;

  R.resize(nLocBas);
  dR_dx.resize(nLocBas);
  dR_dy.resize(nLocBas);
  dR_dz.resize(nLocBas);
  
  Sub_Tan.resize(16); Sub_sur_Tan.resize(16);
  for(int ii=0; ii<16; ++ii)
  {
    Sub_Tan[ii].resize( nLocBas * nLocBas );
    Sub_sur_Tan[ii].resize( snLocBas * snLocBas );
  }

  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  sur_Tangent = new PetscScalar[sur_size * sur_size];
  sur_Residual = new PetscScalar[sur_size];

  Zero_Tangent_Residual();

  Zero_sur_Tangent_Residual();

  print_info();
}



PLocAssem_Tet_VMS_NS_GenAlpha::~PLocAssem_Tet_VMS_NS_GenAlpha()
{
  delete [] Tangent; Tangent = nullptr; 
  delete [] Residual; Residual = nullptr;
  delete [] sur_Tangent; sur_Tangent = nullptr;
  delete [] sur_Residual; sur_Residual = nullptr;
}


void PLocAssem_Tet_VMS_NS_GenAlpha::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("  Three-dimensional Incompressible Navier-Stokes with rho in momentum equations: \n");
  SYS_T::commPrint("  FEM: 4/10-node Tetrahedral \n");
  SYS_T::commPrint("  Spatial: VMS \n");
  SYS_T::commPrint("  Temporal: Generalized-alpha Method \n");
  SYS_T::commPrint("  Density rho = %e \n", rho0);
  SYS_T::commPrint("  Dynamic Viscosity mu = %e \n", vis_mu);
  SYS_T::commPrint("  Kienmatic Viscosity nu = %e \n", vis_mu / rho0);
  SYS_T::commPrint("  Stabilization para CI = %e \n", CI);
  SYS_T::commPrint("  Stabilization para CT = %e \n", CT);
  SYS_T::commPrint("  Backflow Stab. para beta = %e \n", beta);
  SYS_T::commPrint("  Consistent tangent matrix used. \n");
  SYS_T::commPrint("  Nonlinear quadratic term is in advective form. \n");
  SYS_T::commPrint("  Pressure is evaluated at n+alpha_f rather than n+1. \n");
  SYS_T::commPrint("----------------------------------------------------------- \n");
}


void PLocAssem_Tet_VMS_NS_GenAlpha::get_metric(
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


void PLocAssem_Tet_VMS_NS_GenAlpha::get_tau(
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

  tau_c_qua = 1.0 / denom_c;
}


void PLocAssem_Tet_VMS_NS_GenAlpha::get_DC(
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


void PLocAssem_Tet_VMS_NS_GenAlpha::Assem_Residual(
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

  int ii, qua, A, ii4;
  double u, u_t, u_x, u_y, u_z;
  double v, v_t, v_x, v_y, v_z;
  double w, w_t, w_x, w_y, w_z;
  double p, f1, f2, f3;
  double p_x, p_y, p_z;
  double gwts, coor_x, coor_y, coor_z;

  double rx, ry, rz;

  double tau_m, tau_c, tau_dc;

  double u_prime, v_prime, w_prime;

  const double two_mu = 2.0 * vis_mu;
  double NA, NA_x, NA_y, NA_z;
  double velo_dot_gradR, div_vel, r_dot_gradR;
  double tau_m_2;

  double r_dot_gradu, r_dot_gradv, r_dot_gradw;
  double velo_prime_dot_gradR; // v' dot grad NA

  const double curr = time + alpha_f * dt;

  Zero_Residual();

  for(qua=0; qua<nqp; ++qua)
  {
    u = 0.0; u_t = 0.0; u_x = 0.0; u_y = 0.0; u_z = 0.0;
    v = 0.0; v_t = 0.0; v_x = 0.0; v_y = 0.0; v_z = 0.0;
    w = 0.0; w_t = 0.0; w_x = 0.0; w_y = 0.0; w_z = 0.0;
    p = 0.0; coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    p_x = 0.0; p_y = 0.0; p_z = 0.0;

    element->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );
    element->get_invJacobian( qua, dxi_dx );

    for(ii=0; ii<nLocBas; ++ii)
    {
      ii4 = 4 * ii;

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

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    get_tau(tau_m, tau_c, dt, dxi_dx, u, v, w);

    tau_m_2 = tau_m * tau_m;

    gwts = element->get_detJac(qua) * quad->get_qw(qua);

    get_f(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    rx = rho0 * ( u_t + u_x * u + u_y * v + u_z * w - f1 ) + p_x;
    ry = rho0 * ( v_t + v_x * u + v_y * v + v_z * w - f2 ) + p_y;
    rz = rho0 * ( w_t + w_x * u + w_y * v + w_z * w - f3 ) + p_z;

    div_vel = u_x + v_y + w_z;

    u_prime = -1.0 * tau_m * rx;
    v_prime = -1.0 * tau_m * ry;
    w_prime = -1.0 * tau_m * rz;

    get_DC( tau_dc, dxi_dx, u_prime, v_prime, w_prime );

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      velo_dot_gradR = NA_x * u + NA_y * v + NA_z * w;
      r_dot_gradR = NA_x * rx + NA_y * ry + NA_z * rz;
      r_dot_gradu = u_x * rx + u_y * ry + u_z * rz;
      r_dot_gradv = v_x * rx + v_y * ry + v_z * rz;
      r_dot_gradw = w_x * rx + w_y * ry + w_z * rz;
      velo_prime_dot_gradR = NA_x * u_prime + NA_y * v_prime + NA_z * w_prime;

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


// EOF
