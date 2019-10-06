#include "PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha.hpp"

PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha(
    const TimeMethod_GenAlpha * const &tm_gAlpha,
    const int &in_nlocbas, const int &in_nqp,
    const int &in_snlocbas,
    const double &in_rho, const double &in_vis_mu,
    const double &in_beta )
: rho0( in_rho ), vis_mu( in_vis_mu ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), beta(in_beta),
  nLocBas(4), dof_per_node(7), vec_size(16), sur_size(12),
  nqp(in_nqp), snLocBas(in_snlocbas), CI(36.0), CT(4.0)
{
  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  // Set up the surface assembly residual and tangent
  sur_Tangent = new PetscScalar[sur_size * sur_size];
  sur_Residual = new PetscScalar[sur_size];

  Zero_Tangent_Residual();
  
  Zero_sur_Tangent_Residual();

  // Print the infomation on screen
  print_info();
}


PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::~PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha()
{
  delete [] Tangent; Tangent = NULL; delete [] Residual; Residual = NULL;
  delete [] sur_Tangent; sur_Tangent = NULL; 
  delete [] sur_Residual; sur_Residual = NULL;
}


void PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("  Three-dimensional Incompressible Navier-Stokes with rho in momentum equations: \n");
  SYS_T::commPrint("  FEM: 4-node Tetrahedral \n");
  SYS_T::commPrint("  Spatial: ALE-VMS \n");
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
  SYS_T::commPrint("  Input solution vector should have 7 dofs including the mesh motion information in the first three slots. \n");
  SYS_T::commPrint("----------------------------------------------------------- \n");
}


void PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::get_currPts( 
    const double * const &ept_x,
    const double * const &ept_y,
    const double * const &ept_z,
    const double * const &sol )
{
  for(int ii=0; ii<nLocBas; ++ii)
  {
    curPt_x[ii] = ept_x[ii] + sol[7*ii];
    curPt_y[ii] = ept_y[ii] + sol[7*ii+1];
    curPt_z[ii] = ept_z[ii] + sol[7*ii+2];
  }
}


void PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::get_currBCPts( 
    const double * const &ept_x,
    const double * const &ept_y,
    const double * const &ept_z,
    const double * const &sol )
{
  for(int ii=0; ii<snLocBas; ++ii)
  {
    curPt_x[ii] = ept_x[ii] + sol[7*ii];
    curPt_y[ii] = ept_y[ii] + sol[7*ii+1];
    curPt_z[ii] = ept_z[ii] + sol[7*ii+2];
  }
}

void PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::get_metric(
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


void PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::get_tau( 
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


void PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::get_DC( 
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


void PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::Assem_Residual(
    const double &time, const double &dt,
    const double * const &velo,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  int ii, qua, A, ii7;
  double u, u_t, u_x, u_y, u_z;
  double v, v_t, v_x, v_y, v_z;
  double w, w_t, w_x, w_y, w_z;
  double p, f1, f2, f3;
  double p_x, p_y, p_z;
  double gwts, coor_x, coor_y, coor_z;

  double mu, mv, mw; // mesh velocity, i.e. hat(v)
  double cu, cv, cw; // v - hat(v)

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
    mu = 0.0; mv = 0.0; mw = 0.0;

    element->get_R_gradR( qua, R, dR_dx, dR_dy, dR_dz );
    element->get_invJacobian( qua, dxi_dx );

    for(ii=0; ii<nLocBas; ++ii)
    {
      ii7 = 7 * ii;

      mu += velo[ii7+0] * R[ii];
      mv += velo[ii7+1] * R[ii];
      mw += velo[ii7+2] * R[ii];

      u_t += velo[ii7+4] * R[ii];
      v_t += velo[ii7+5] * R[ii];
      w_t += velo[ii7+6] * R[ii];

      u += disp[ii7+4] * R[ii];
      v += disp[ii7+5] * R[ii];
      w += disp[ii7+6] * R[ii];
      p += disp[ii7+3] * R[ii];

      u_x += disp[ii7+4] * dR_dx[ii];
      v_x += disp[ii7+5] * dR_dx[ii];
      w_x += disp[ii7+6] * dR_dx[ii];
      p_x += disp[ii7+3] * dR_dx[ii];

      u_y += disp[ii7+4] * dR_dy[ii];
      v_y += disp[ii7+5] * dR_dy[ii];
      w_y += disp[ii7+6] * dR_dy[ii];
      p_y += disp[ii7+3] * dR_dy[ii];

      u_z += disp[ii7+4] * dR_dz[ii];
      v_z += disp[ii7+5] * dR_dz[ii];
      w_z += disp[ii7+6] * dR_dz[ii];
      p_z += disp[ii7+3] * dR_dz[ii];

      coor_x += curPt_x[ii] * R[ii];
      coor_y += curPt_y[ii] * R[ii];
      coor_z += curPt_z[ii] * R[ii];
    }

    cu = u - mu;
    cv = v - mv;
    cw = w - mw;

    get_tau(tau_m, tau_c, dt, dxi_dx, cu, cv, cw);

    tau_m_2 = tau_m * tau_m;

    gwts = element->get_detJac(qua) * quad->get_qw(qua);

    get_f(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    rx = rho0 * ( u_t + u_x * cu + u_y * cv + u_z * cw - f1 ) + p_x;
    ry = rho0 * ( v_t + v_x * cu + v_y * cv + v_z * cw - f2 ) + p_y;
    rz = rho0 * ( w_t + w_x * cu + w_y * cv + w_z * cw - f3 ) + p_z;

    div_vel = u_x + v_y + w_z;

    u_prime = -1.0 * tau_m * rx;
    v_prime = -1.0 * tau_m * ry;
    w_prime = -1.0 * tau_m * rz;

    get_DC( tau_dc, dxi_dx, u_prime, v_prime, w_prime );

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      velo_dot_gradR = NA_x * cu + NA_y * cv + NA_z * cw;
      r_dot_gradR = NA_x * rx + NA_y * ry + NA_z * rz;
      r_dot_gradu = u_x * rx + u_y * ry + u_z * rz;
      r_dot_gradv = v_x * rx + v_y * ry + v_z * rz;
      r_dot_gradw = w_x * rx + w_y * ry + w_z * rz;
      velo_prime_dot_gradR = NA_x * u_prime + NA_y * v_prime + NA_z * w_prime;

      Residual[4*A] += gwts * ( NA * div_vel + tau_m * r_dot_gradR );

      Residual[4*A+1] += gwts * ( NA * rho0 * u_t
          + NA * rho0 * (cu * u_x + cv * u_y + cw * u_z)
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
          + NA * rho0 * (cu * v_x + cv * v_y + cw * v_z)
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
          + NA * rho0 * (cu * w_x + cv * w_y + cw * w_z)
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


void PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &velo,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  int ii, qua, A, ii7, B, jj, index;
  double u, u_t, u_x, u_y, u_z;
  double v, v_t, v_x, v_y, v_z;
  double w, w_t, w_x, w_y, w_z;
  double p, f1, f2, f3, p_x, p_y, p_z;
  double gwts, coor_x, coor_y, coor_z;

  double mu, mv, mw; // mesh velocity: hat(v)
  double cu, cv, cw; // c = v - hat(v)

  double rx, ry, rz; // Momentum residual

  double tau_m, tau_c, tau_dc;

  double u_prime, v_prime, w_prime;

  const double two_mu = 2.0 * vis_mu;
  const double rho0_2 = rho0 * rho0;
  double NA, NA_x, NA_y, NA_z;
  double velo_dot_gradR, div_vel, r_dot_gradR;
  double tau_m_2;
  double r_dot_gradu, r_dot_gradv, r_dot_gradw;

  double velo_prime_dot_gradR;

  double NB, NB_x, NB_y, NB_z;
  double NANB, NAxNB, NAyNB, NAzNB;
  double NANBx, NAxNBx, NAyNBx, NAzNBx;
  double NANBy, NAxNBy, NAyNBy, NAzNBy;
  double NANBz, NAxNBz, NAyNBz, NAzNBz;

  double drx_du_B, drx_dv_B, drx_dw_B;
  double dry_du_B, dry_dv_B, dry_dw_B;
  double drz_du_B, drz_dv_B, drz_dw_B;

  double velo_dot_gradNB, velo_prime_dot_gradNB;

  const double curr = time + alpha_f * dt;

  const double dd_dv = alpha_f * gamma * dt;

  Zero_Tangent_Residual();

  Zero_Sub_Tan();

  for(qua=0; qua<nqp; ++qua)
  {
    u = 0.0; u_t = 0.0; u_x = 0.0; u_y = 0.0; u_z = 0.0;
    v = 0.0; v_t = 0.0; v_x = 0.0; v_y = 0.0; v_z = 0.0;
    w = 0.0; w_t = 0.0; w_x = 0.0; w_y = 0.0; w_z = 0.0;
    p = 0.0; coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    p_x = 0.0; p_y = 0.0; p_z = 0.0;
    mu = 0.0; mv = 0.0; mw = 0.0;

    element->get_R_gradR( qua, R, dR_dx, dR_dy, dR_dz );
    element->get_invJacobian( qua, dxi_dx );

    for(ii=0; ii<nLocBas; ++ii)
    {
      ii7 = 7 * ii;

      mu += velo[ii7+0] * R[ii];
      mv += velo[ii7+1] * R[ii];
      mw += velo[ii7+2] * R[ii];

      u_t += velo[ii7+4] * R[ii];
      v_t += velo[ii7+5] * R[ii];
      w_t += velo[ii7+6] * R[ii];

      u += disp[ii7+4] * R[ii];
      v += disp[ii7+5] * R[ii];
      w += disp[ii7+6] * R[ii];
      p += disp[ii7+3] * R[ii];

      u_x += disp[ii7+4] * dR_dx[ii];
      v_x += disp[ii7+5] * dR_dx[ii];
      w_x += disp[ii7+6] * dR_dx[ii];
      p_x += disp[ii7+3] * dR_dx[ii];

      u_y += disp[ii7+4] * dR_dy[ii];
      v_y += disp[ii7+5] * dR_dy[ii];
      w_y += disp[ii7+6] * dR_dy[ii];
      p_y += disp[ii7+3] * dR_dy[ii];

      u_z += disp[ii7+4] * dR_dz[ii];
      v_z += disp[ii7+5] * dR_dz[ii];
      w_z += disp[ii7+6] * dR_dz[ii];
      p_z += disp[ii7+3] * dR_dz[ii];

      coor_x += curPt_x[ii] * R[ii];
      coor_y += curPt_y[ii] * R[ii];
      coor_z += curPt_z[ii] * R[ii];
    }

    cu = u - mu;
    cv = v - mv;
    cw = w - mw;

    get_tau(tau_m, tau_c, dt, dxi_dx, cu, cv, cw);

    tau_m_2 = tau_m * tau_m;

    gwts = element->get_detJac(qua) * quad->get_qw(qua); 

    get_f(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    rx = rho0 * ( u_t + u_x * cu + u_y * cv + u_z * cw - f1 ) + p_x;
    ry = rho0 * ( v_t + v_x * cu + v_y * cv + v_z * cw - f2 ) + p_y;
    rz = rho0 * ( w_t + w_x * cu + w_y * cv + w_z * cw - f3 ) + p_z;

    div_vel = u_x + v_y + w_z;

    u_prime = -1.0 * tau_m * rx;
    v_prime = -1.0 * tau_m * ry;
    w_prime = -1.0 * tau_m * rz;

    get_DC( tau_dc, dxi_dx, u_prime, v_prime, w_prime );

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      velo_dot_gradR = NA_x * cu + NA_y * cv + NA_z * cw;
      r_dot_gradR = NA_x * rx + NA_y * ry + NA_z * rz;
      r_dot_gradu = u_x * rx + u_y * ry + u_z * rz;
      r_dot_gradv = v_x * rx + v_y * ry + v_z * rz;
      r_dot_gradw = w_x * rx + w_y * ry + w_z * rz;
      velo_prime_dot_gradR = NA_x * u_prime + NA_y * v_prime + NA_z * w_prime;

      Residual[4*A] += gwts * ( NA * div_vel + tau_m * r_dot_gradR );

      Residual[4*A+1] += gwts * ( NA * rho0 * u_t
          + NA * rho0 * (cu * u_x + cv * u_y + cw * u_z)
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
          + NA * rho0 * (cu * v_x + cv * v_y + cw * v_z)
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
          + NA * rho0 * (cu * w_x + cv * w_y + cw * w_z)
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

      for(B=0; B<nLocBas; ++B)
      {
        index = B + A * nLocBas;
        NB = R[B]; NB_x = dR_dx[B]; NB_y = dR_dy[B]; NB_z = dR_dz[B];
        velo_dot_gradNB = cu * NB_x + cv * NB_y + cw * NB_z;
        velo_prime_dot_gradNB = u_prime * NB_x + v_prime * NB_y + w_prime * NB_z;

        NANB  = NA*NB; NANBx = NA*NB_x; NANBy = NA*NB_y; NANBz = NA*NB_z;
        NAxNB = NA_x*NB; NAxNBx = NA_x*NB_x; NAxNBy = NA_x*NB_y; NAxNBz = NA_x*NB_z;
        NAyNB = NA_y*NB; NAyNBx = NA_y*NB_x; NAyNBy = NA_y*NB_y; NAyNBz = NA_y*NB_z;
        NAzNB = NA_z*NB; NAzNBx = NA_z*NB_x; NAzNBy = NA_z*NB_y; NAzNBz = NA_z*NB_z;

        drx_du_B = rho0 * ( u_x * NB + velo_dot_gradNB ); 
        drx_dv_B = rho0 * u_y * NB;
        drx_dw_B = rho0 * u_z * NB;

        dry_du_B = rho0 * v_x * NB;
        dry_dv_B = rho0 * ( v_y * NB + velo_dot_gradNB );
        dry_dw_B = rho0 * v_z * NB;

        drz_du_B = rho0 * w_x * NB;
        drz_dv_B = rho0 * w_y * NB;
        drz_dw_B = rho0 * (w_z * NB + velo_dot_gradNB);

        // Continuity equation with respect to p, u, v, w
        Sub_Tan[0][index] += gwts * dd_dv * tau_m * (NAxNBx + NAyNBy + NAzNBz);

        Sub_Tan[1][index] += gwts * ( alpha_m * tau_m * rho0 * NAxNB
            + dd_dv * ( NANBx + tau_m * NA_x * drx_du_B
              + tau_m * NA_y * dry_du_B + tau_m * NA_z * drz_du_B ) );

        Sub_Tan[2][index] += gwts * ( alpha_m * tau_m * rho0 * NAyNB
            + dd_dv * ( NANBy + tau_m * NA_x * drx_dv_B
              + tau_m * NA_y * dry_dv_B + tau_m * NA_z * drz_dv_B ) );

        Sub_Tan[3][index] += gwts * ( alpha_m * tau_m * rho0 * NAzNB
            + dd_dv * ( NANBz + tau_m * NA_x * drx_dw_B
              + tau_m * NA_y * dry_dw_B + tau_m * NA_z * drz_dw_B ) );

        // Momentum-x with respect to p, u, v, w
        Sub_Tan[4][index] += gwts * dd_dv * ((-1.0) * NAxNB
            + velo_dot_gradR * tau_m * rho0 * NB_x
            - NA * tau_m * rho0 * (u_x * NB_x + u_y * NB_y + u_z * NB_z)
            - 2.0 * tau_m_2 * rho0 * rx * NAxNBx
            - tau_m_2 * rho0 * NA_y * (rx * NB_y + ry * NB_x)
            - tau_m_2 * rho0 * NA_z * (rx * NB_z + rz * NB_x) );

        Sub_Tan[5][index] += gwts * ( 
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

        Sub_Tan[6][index] += gwts * ( 
            alpha_m * (-1.0) * rho0_2 * (tau_m * u_y * NANB + tau_m_2 * rx * NAyNB)
            + dd_dv * ( NANB * rho0 * u_y + vis_mu * NAyNBx 
              + rho0 * tau_m * rx * NAyNB
              + velo_dot_gradR * rho0 * tau_m * drx_dv_B
              - rho0 * tau_m * NA * (u_x*drx_dv_B + u_y*dry_dv_B + u_z*drz_dv_B)
              + tau_c * NAxNBy
              - 2.0 * rho0 * tau_m_2 * rx * NA_x * drx_dv_B
              - rho0 * tau_m_2 * NA_y * (rx * dry_dv_B + ry * drx_dv_B)
              - rho0 * tau_m_2 * NA_z * (rx * drz_dv_B + rz * drx_dv_B) ) );

        Sub_Tan[7][index] += gwts * (
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
        Sub_Tan[8][index] += gwts * dd_dv * ( (-1.0) * NAyNB
            + velo_dot_gradR * tau_m * rho0 * NB_y
            - NA * tau_m * rho0 * (v_x * NB_x + v_y * NB_y + v_z * NB_z)
            - tau_m_2 * rho0 * NA_x * (rx * NB_y + ry * NB_x)
            - 2.0 * tau_m_2 * rho0 * ry * NAyNBy
            - tau_m_2 * rho0 * NA_z * (ry * NB_z + rz * NB_y) );

        Sub_Tan[9][index] += gwts * (
            alpha_m * (-1.0) * rho0_2 * (tau_m * v_x * NANB + tau_m_2 * ry * NAxNB)
            + dd_dv * ( NANB * rho0 * v_x + vis_mu * NAxNBy
              + rho0 * tau_m * ry * NAxNB
              + velo_dot_gradR * rho0 * tau_m * dry_du_B
              - rho0 * tau_m * NA * (v_x*drx_du_B + v_y*dry_du_B + v_z*drz_du_B)
              + tau_c * NAyNBx
              - rho0 * tau_m_2 * NA_x * (ry * drx_du_B + rx * dry_du_B)
              - 2.0 * rho0 * tau_m_2 * ry * NA_y * dry_du_B
              - rho0 * tau_m_2 * NA_z * (ry * drz_du_B + rz * dry_du_B) ) );

        Sub_Tan[10][index] += gwts * (
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

        Sub_Tan[11][index] += gwts * (
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
        Sub_Tan[12][index] += gwts * dd_dv * ( (-1.0) * NAzNB
            + velo_dot_gradR * tau_m * rho0 * NB_z
            - NA * tau_m * rho0 * (w_x * NB_x + w_y * NB_y + w_z * NB_z)
            - tau_m_2 * rho0 * NA_x * (rx * NB_z + rz * NB_x)
            - tau_m_2 * rho0 * NA_y * (ry * NB_z + rz * NB_y)
            - 2.0 * tau_m_2 * rho0 * rz * NAzNBz );

        Sub_Tan[13][index] += gwts * (
            alpha_m * (-1.0) * rho0_2 * (tau_m * w_x * NANB + tau_m_2 * rz * NAxNB)
            + dd_dv * ( NANB * rho0 * w_x + vis_mu * NAxNBz
              + rho0 * tau_m * rz * NAxNB
              + velo_dot_gradR * rho0 * tau_m * drz_du_B
              - rho0 * tau_m * NA * (w_x*drx_du_B + w_y*dry_du_B + w_z*drz_du_B)
              + tau_c * NAzNBx
              - rho0 * tau_m_2 * NA_x * (rx * drz_du_B + rz * drx_du_B)
              - rho0 * tau_m_2 * NA_y * (ry * drz_du_B + rz * dry_du_B)
              - 2.0 * rho0 * tau_m_2 * rz * NA_z * drz_du_B ) );

        Sub_Tan[14][index] += gwts * (
            alpha_m * (-1.0) * rho0_2 * (tau_m * w_y * NANB + tau_m_2 * rz * NAyNB)
            + dd_dv * ( NANB * rho0 * w_y + vis_mu * NAyNBz
              + rho0 * tau_m * rz * NAyNB
              + velo_dot_gradR * rho0 * tau_m * drz_dv_B
              - rho0 * tau_m * NA * (w_x*drx_dv_B + w_y*dry_dv_B + w_z*drz_dv_B)
              + tau_c * NAzNBy
              - rho0 * tau_m_2 * NA_x * (rx * drz_dv_B + rz * drx_dv_B)
              - rho0 * tau_m_2 * NA_y * (ry * drz_dv_B + rz * dry_dv_B)
              - 2.0 * rho0 * tau_m_2 * rz * NA_z * drz_dv_B ) );

        Sub_Tan[15][index] += gwts * (
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

  for(ii=0; ii<4; ++ii)
  {
    for(jj=0; jj<4; ++jj)
    {
      for(A=0; A<nLocBas; ++A)
      {
        for(B=0; B<nLocBas; ++B)
        {
          Tangent[ 4*nLocBas*(4*A+ii) + 4*B + jj  ] =
            Sub_Tan[ii*4+jj][A*nLocBas + B];
        }
      }
    }
  }
}


void PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::Assem_Mass_Residual(
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  int ii, jj, qua, A, B, index;
  double u, u_x, u_y, u_z;
  double v, v_x, v_y, v_z;
  double w, w_x, w_y, w_z;
  double p, f1, f2, f3;
  double gwts, coor_x, coor_y, coor_z;

  double NA, NA_x, NA_y, NA_z;

  const double two_mu = 2.0 * vis_mu;

  double curr = 0.0;

  Zero_Tangent_Residual();

  Zero_Sub_Tan();

  for(qua=0; qua<nqp; ++qua)
  {
    u = 0.0; u_x = 0.0; u_y = 0.0; u_z = 0.0;
    v = 0.0; v_x = 0.0; v_y = 0.0; v_z = 0.0;
    w = 0.0; w_x = 0.0; w_y = 0.0; w_z = 0.0;
    p = 0.0; coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    element->get_R_gradR( qua, R, dR_dx, dR_dy, dR_dz );

    for(ii=0; ii<nLocBas; ++ii)
    {
      u += disp[7*ii+4] * R[ii];
      v += disp[7*ii+5] * R[ii];
      w += disp[7*ii+6] * R[ii];
      p += disp[7*ii+3] * R[ii];

      u_x += disp[7*ii+4] * dR_dx[ii];
      v_x += disp[7*ii+5] * dR_dx[ii];
      w_x += disp[7*ii+6] * dR_dx[ii];

      u_y += disp[7*ii+4] * dR_dy[ii];
      v_y += disp[7*ii+5] * dR_dy[ii];
      w_y += disp[7*ii+6] * dR_dy[ii];

      u_z += disp[7*ii+4] * dR_dz[ii];
      v_z += disp[7*ii+5] * dR_dz[ii];
      w_z += disp[7*ii+6] * dR_dz[ii];

      coor_x += curPt_x[ii] * R[ii];
      coor_y += curPt_y[ii] * R[ii];
      coor_z += curPt_z[ii] * R[ii];
    }

    gwts = element->get_detJac(qua) * quad->get_qw(qua);

    get_f(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      Residual[4*A+1] += gwts * ( NA * rho0 * (u*u_x + v*v_y + w*u_z) 
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

      for(B=0; B<nLocBas; ++B)
      {
        index = A * nLocBas + B;
        Sub_Tan[0][index]  += gwts * rho0 * NA * R[B];
        Sub_Tan[5][index]  += gwts * rho0 * NA * R[B];
        Sub_Tan[10][index] += gwts * rho0 * NA * R[B];
        Sub_Tan[15][index] += gwts * rho0 * NA * R[B];
      }
    }
  }

  for(ii=0; ii<4; ++ii)
  {
    for(jj=0; jj<4; ++jj)
    {
      for(A=0; A<nLocBas; ++A)
      {
        for(B=0; B<nLocBas; ++B)
        {
          Tangent[ 4*nLocBas*(4*A+ii) + 4*B + jj ] =
            Sub_Tan[ii*4+jj][A*nLocBas + B];
        }
      }
    }
  }
}


void PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::Assem_Residual_EBC(
    const int &ebc_id,
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  get_currBCPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, vec_b);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const int face_nqp = quad -> get_num_quadPts();

  int ii, qua, A;
  double gwts, coor_x, coor_y, coor_z, gx, gy, gz, nx, ny, nz, surface_area;
  const double curr = time + alpha_f * dt;

  Zero_Residual();

  for(qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, R);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    for(ii=0; ii<snLocBas; ++ii)
    {
      coor_x += curPt_x[ii] * R[ii];
      coor_y += curPt_y[ii] * R[ii];
      coor_z += curPt_z[ii] * R[ii];
    }

    get_ebc_fun( ebc_id, coor_x, coor_y, coor_z, curr, nx, ny, nz,
        gx, gy, gz );

    gwts = surface_area * quad -> get_qw(qua);

    for(A=0; A<snLocBas; ++A)
    {
      Residual[4*A+1] -= gwts * R[A] * gx;
      Residual[4*A+2] -= gwts * R[A] * gy;
      Residual[4*A+3] -= gwts * R[A] * gz;
    }
  }
}


double PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::get_flowrate( 
    const double * const &vec,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  get_currBCPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, vec);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const int face_nqp = quad -> get_num_quadPts();

  int ii, qua;
  double gwts, nx, ny, nz, surface_area, u, v, w;

  double flrate = 0.0;

  for(qua =0; qua< face_nqp; ++qua)
  {
    element->get_R(qua, R);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);

    u = 0.0; v = 0.0; w = 0.0;
    for(ii=0; ii<snLocBas; ++ii)
    {
      u += vec[7*ii+4] * R[ii];
      v += vec[7*ii+5] * R[ii];
      w += vec[7*ii+6] * R[ii];
    }
    gwts = surface_area * quad->get_qw(qua);
    flrate += gwts * ( u * nx + v * ny + w * nz );
  }

  return flrate;
}


void PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::get_pressure_area( 
    const double * const &vec,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad,
    double &pres, double &area )
{
  get_currBCPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, vec);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const int face_nqp = quad -> get_num_quadPts();

  int ii, qua;
  double gwts, nx, ny, nz, surface_area, pp;

  // Initialize the two variables to be passed out
  pres = 0.0;
  area = 0.0;

  for(qua =0; qua< face_nqp; ++qua)
  {
    element->get_R(qua, R);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);

    pp = 0.0;
    for(ii=0; ii<snLocBas; ++ii) pp += vec[7*ii+3] * R[ii];
    
    gwts = surface_area * quad->get_qw(qua);
    
    pres += gwts * pp; 
    area += gwts;
  }
}


void PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::Assem_Residual_EBC_Resistance(
    const int &ebc_id,
    const double &val,
    const double * const &vec_b,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  get_currBCPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, vec_b);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const int face_nqp = quad -> get_num_quadPts();

  double gwts, nx, ny, nz, surface_area;

  Zero_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, R);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);
    gwts = surface_area * quad -> get_qw(qua);

    for(int A=0; A<snLocBas; ++A)
    {
      Residual[4*A+1] += gwts * R[A] * nx * val;
      Residual[4*A+2] += gwts * R[A] * ny * val;
      Residual[4*A+3] += gwts * R[A] * nz * val;
    }
  }
}


void PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::Assem_Residual_BackFlowStab(
    const double * const &velo,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  get_currBCPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const int face_nqp = quad -> get_num_quadPts();

  int ii, qua, A;
  double gwts, nx, ny, nz, surface_area;
  double u, v, w, mu, mv, mw, cu, cv, cw;
  double factor;
  
  Zero_sur_Residual();

  for(qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, R);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);

    u = 0.0; v = 0.0; w = 0.0; mu = 0.0; mv = 0.0; mw = 0.0;
    for(ii=0; ii<snLocBas; ++ii)
    {
      mu += velo[ii*7+0] * R[ii];
      mv += velo[ii*7+1] * R[ii];
      mw += velo[ii*7+2] * R[ii];

      u += disp[ii*7+4] * R[ii];
      v += disp[ii*7+5] * R[ii];
      w += disp[ii*7+6] * R[ii];
    }

    cu = u - mu;
    cv = v - mv;
    cw = w - mw;

    const double temp = cu * nx + cv * ny + cw * nz;

    if(temp < 0.0) factor = temp * rho0 * beta;
    else factor = 0.0;

    gwts = surface_area * quad -> get_qw(qua);

    for(A=0; A<snLocBas; ++A)
    {
      sur_Residual[4*A+1] -= gwts * R[A] * factor * u;
      sur_Residual[4*A+2] -= gwts * R[A] * factor * v;
      sur_Residual[4*A+3] -= gwts * R[A] * factor * w;
    }
  }
}


void PLocAssem_Tet4_ALE_VMS_NS_mom_3D_GenAlpha::Assem_Tangent_Residual_BackFlowStab(
    const double &dt,
    const double * const &velo,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  get_currBCPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const int face_nqp = quad -> get_num_quadPts();

  int ii, qua, A, B;
  double gwts, nx, ny, nz, surface_area;
  double u, v, w, mu, mv, mw, cu, cv, cw;
  double factor;
  
  const double dd_dv = alpha_f * gamma * dt;

  Zero_sur_Tangent_Residual();

  Zero_Sub_sur_Tan();

  for(qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, R);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);

    u = 0.0; v = 0.0; w = 0.0; mu = 0.0; mv = 0.0; mw = 0.0;
    for(ii=0; ii<snLocBas; ++ii)
    {
      mu += velo[ii*7+0] * R[ii];
      mv += velo[ii*7+1] * R[ii];
      mw += velo[ii*7+2] * R[ii];

      u += disp[ii*7+4] * R[ii];
      v += disp[ii*7+5] * R[ii];
      w += disp[ii*7+6] * R[ii];
    }

    cu = u - mu;
    cv = v - mv;
    cw = w - mw;

    const double temp = cu * nx + cv * ny + cw * nz;

    if(temp < 0.0) factor = temp * rho0 * beta;
    else factor = 0.0;

    gwts = surface_area * quad -> get_qw(qua);

    // snLocBas = 3 for linear tet/tri element
    for(A=0; A<snLocBas; ++A)
    {
      sur_Residual[4*A+1] -= gwts * R[A] * factor * u;
      sur_Residual[4*A+2] -= gwts * R[A] * factor * v;
      sur_Residual[4*A+3] -= gwts * R[A] * factor * w;

      for(B=0; B<snLocBas; ++B)
      {
        const int index = B + A * snLocBas; // index here ranges 0 to 8

        Sub_sur_Tan[5][index]  -= gwts * dd_dv * R[A] * factor * R[B];
        Sub_sur_Tan[10][index] -= gwts * dd_dv * R[A] * factor * R[B];
        Sub_sur_Tan[15][index] -= gwts * dd_dv * R[A] * factor * R[B];
      }
    }
  }
  
  for(A=0; A<snLocBas; ++A)
  {
    for(B=0; B<snLocBas; ++B)
    {
      // ii = jj = 1
      sur_Tangent[ 4*snLocBas*(4*A+1) + 4*B+1 ] = Sub_sur_Tan[5][A*snLocBas + B];

      // ii = jj = 2
      sur_Tangent[ 4*snLocBas*(4*A+2) + 4*B+2 ] = Sub_sur_Tan[10][A*snLocBas + B];
      
      // ii = jj = 3
      sur_Tangent[ 4*snLocBas*(4*A+3) + 4*B+3 ] = Sub_sur_Tan[15][A*snLocBas + B];
    }
  }
}

// EOF
