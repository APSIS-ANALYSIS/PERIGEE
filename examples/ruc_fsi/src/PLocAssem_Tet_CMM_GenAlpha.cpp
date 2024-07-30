#include "PLocAssem_Tet_CMM_GenAlpha.hpp"

PLocAssem_Tet_CMM_GenAlpha::PLocAssem_Tet_CMM_GenAlpha(
    const TimeMethod_GenAlpha * const &tm_gAlpha,
    const int &in_nqp, const int &in_face_nqp,
    const double &in_rho, const double &in_vis_mu,
    const double &in_beta, const double &in_wall_rho,
    const double &in_nu, const double &in_kappa,
    const double &in_ctauc, const int &elemtype )
: rho0( in_rho ), vis_mu( in_vis_mu ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), beta(in_beta), rho_w(in_wall_rho),
  nu_w(in_nu), kappa_w(in_kappa), nqp(in_nqp), face_nqp(in_face_nqp), 
  Ctauc( in_ctauc )
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
  SYS_T::commPrint("  Three-dimensional FSI with membrane wall: \n");
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
  SYS_T::commPrint("  Wall density = %e \n", rho_w);
  SYS_T::commPrint("  Wall Poisson ratio = %e \n", nu_w);
  SYS_T::commPrint("  Wall transverse shearing moduli = %e \n", kappa_w);
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


SymmTensor2_3D PLocAssem_Tet_CMM_GenAlpha::get_metric(
    const std::array<double,9> &f ) const
{
  // PHASTA definition 
  const double coef = 0.6299605249474365;

  const double fk0 = 2.0 * f[0] + (f[3] + f[6]);
  const double fk1 = 2.0 * f[3] + (f[0] + f[6]);
  const double fk2 = 2.0 * f[6] + (f[0] + f[3]);
  const double fk3 = 2.0 * f[1] + (f[4] + f[7]);
  const double fk4 = 2.0 * f[4] + (f[1] + f[7]);
  const double fk5 = 2.0 * f[7] + (f[1] + f[4]);
  const double fk6 = 2.0 * f[2] + (f[5] + f[8]);
  const double fk7 = 2.0 * f[5] + (f[2] + f[8]);
  const double fk8 = 2.0 * f[8] + (f[2] + f[5]);

  return SymmTensor2_3D(coef * ( fk0 * f[0] + fk1 * f[3] + fk2 * f[6] ),
  coef * ( fk3 * f[1] + fk4 * f[4] + fk5 * f[7] ),
  coef * ( fk6 * f[2] + fk7 * f[5] + fk8 * f[8] ),
  coef * ( fk3 * f[2] + fk4 * f[5] + fk5 * f[8] ),
  coef * ( fk0 * f[2] + fk1 * f[5] + fk2 * f[8] ),
  coef * ( fk0 * f[1] + fk1 * f[4] + fk2 * f[7] ));
}


std::array<double, 2> PLocAssem_Tet_CMM_GenAlpha::get_tau(
    const double &dt, const std::array<double,9> &dxi_dx,
    const double &u, const double &v, const double &w ) const
{
  // Use K matrix to correct the metric
  const SymmTensor2_3D G = get_metric( dxi_dx );

  const Vector_3 velo_vec( u, v, w );

  const double temp_nu = vis_mu / rho0;

  const double denom_m = std::sqrt(CT / (dt*dt) + 
  G.VecMatVec( velo_vec, velo_vec ) + 
  CI * temp_nu * temp_nu * G.MatContraction( G ));

  return {{1.0 / ( rho0 * denom_m ), Ctauc * rho0 * denom_m / G.tr()}};
}


double PLocAssem_Tet_CMM_GenAlpha::get_DC(
    const std::array<double, 9> &dxi_dx,
    const double &u, const double &v, const double &w ) const
{
  // const SymmTensor2_3D G = get_metric( dxi_dx );
  // const Vector_3 velo_vec( u, v, w );
  // double dc_tau = G.VecMatVec( velo_vec, velo_vec );

  // if(dc_tau > 1.0e-15) dc_tau = rho0 * std::pow(dc_tau, -0.5);
  // else dc_tau = 0.0;

  const double dc_tau = 0.0;

  return dc_tau;
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

  const double two_mu = 2.0 * vis_mu;

  const double curr = time + alpha_f * dt;

  Zero_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);
  std::vector<double> d2R_dxx(nLocBas, 0.0), d2R_dyy(nLocBas, 0.0), d2R_dzz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    double u = 0.0, u_t = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v = 0.0, v_t = 0.0, v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w = 0.0, w_t = 0.0, w_x = 0.0, w_y = 0.0, w_z = 0.0;
    double p = 0.0, p_x = 0.0, p_y = 0.0, p_z = 0.0;

    double u_xx = 0.0, u_yy = 0.0, u_zz = 0.0;
    double v_xx = 0.0, v_yy = 0.0, v_zz = 0.0;
    double w_xx = 0.0, w_yy = 0.0, w_zz = 0.0;

    Vector_3 coor( 0.0, 0.0, 0.0 );

    element->get_3D_R_gradR_LaplacianR( qua, &R[0], &dR_dx[0], 
        &dR_dy[0], &dR_dz[0], &d2R_dxx[0], &d2R_dyy[0], &d2R_dzz[0] );

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

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    // Get the tau_m and tau_c
    const auto dxi_dx = element->get_invJacobian(qua);

    const std::array<double, 2> tau = get_tau(dt, dxi_dx, u, v, w);
    const double tau_m = tau[0];
    const double tau_c = tau[1];

    const double tau_m_2 = tau_m * tau_m;

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    // Get the body force
    const Vector_3 f_body = get_f(coor, curr);

    const double u_lap = u_xx + u_yy + u_zz;
    const double v_lap = v_xx + v_yy + v_zz;
    const double w_lap = w_xx + w_yy + w_zz;

    const double rx = rho0 * ( u_t + u_x * u + u_y * v + u_z * w - f_body.x() ) + p_x - vis_mu * u_lap;
    const double ry = rho0 * ( v_t + v_x * u + v_y * v + v_z * w - f_body.y() ) + p_y - vis_mu * v_lap;
    const double rz = rho0 * ( w_t + w_x * u + w_y * v + w_z * w - f_body.z() ) + p_z - vis_mu * w_lap;

    const double div_vel = u_x + v_y + w_z;

    const double u_prime = -1.0 * tau_m * rx;
    const double v_prime = -1.0 * tau_m * ry;
    const double w_prime = -1.0 * tau_m * rz;

    const double r_dot_gradu = u_x * rx + u_y * ry + u_z * rz;
    const double r_dot_gradv = v_x * rx + v_y * ry + v_z * rz;
    const double r_dot_gradw = w_x * rx + w_y * ry + w_z * rz;

    // Get the Discontinuity Capturing tau
    const double tau_dc = get_DC( dxi_dx, u_prime, v_prime, w_prime );

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
          - NA * rho0 * f_body.x() );

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
          - NA * rho0 * f_body.y() );

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
          - NA * rho0 * f_body.z() );
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

  const double two_mu = 2.0 * vis_mu;

  const double rho0_2 = rho0 * rho0;

  const double curr = time + alpha_f * dt;

  const double dd_dv = alpha_f * gamma * dt;

  Zero_Tangent_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);
  std::vector<double> d2R_dxx(nLocBas, 0.0), d2R_dyy(nLocBas, 0.0), d2R_dzz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    double u = 0.0, u_t = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v = 0.0, v_t = 0.0, v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w = 0.0, w_t = 0.0, w_x = 0.0, w_y = 0.0, w_z = 0.0;
    double p = 0.0, p_x = 0.0, p_y = 0.0, p_z = 0.0;

    double u_xx = 0.0, u_yy = 0.0, u_zz = 0.0;
    double v_xx = 0.0, v_yy = 0.0, v_zz = 0.0;
    double w_xx = 0.0, w_yy = 0.0, w_zz = 0.0;

    Vector_3 coor( 0.0, 0.0, 0.0 );

    element->get_3D_R_gradR_LaplacianR( qua, &R[0], &dR_dx[0], 
        &dR_dy[0], &dR_dz[0], &d2R_dxx[0], &d2R_dyy[0], &d2R_dzz[0] );

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

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    const auto dxi_dx = element->get_invJacobian(qua);

    const std::array<double, 2> tau = get_tau(dt, dxi_dx, u, v, w);
    const double tau_m = tau[0];
    const double tau_c = tau[1];

    const double tau_m_2 = tau_m * tau_m;

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua); 

    const Vector_3 f_body = get_f(coor, curr);

    const double u_lap = u_xx + u_yy + u_zz;
    const double v_lap = v_xx + v_yy + v_zz;
    const double w_lap = w_xx + w_yy + w_zz;

    const double rx = rho0 * ( u_t + u_x * u + u_y * v + u_z * w - f_body.x() ) + p_x - vis_mu * u_lap;
    const double ry = rho0 * ( v_t + v_x * u + v_y * v + v_z * w - f_body.y() ) + p_y - vis_mu * v_lap ;
    const double rz = rho0 * ( w_t + w_x * u + w_y * v + w_z * w - f_body.z() ) + p_z - vis_mu * w_lap;

    const double div_vel = u_x + v_y + w_z;

    const double u_prime = -1.0 * tau_m * rx;
    const double v_prime = -1.0 * tau_m * ry;
    const double w_prime = -1.0 * tau_m * rz;

    const double r_dot_gradu = u_x * rx + u_y * ry + u_z * rz;
    const double r_dot_gradv = v_x * rx + v_y * ry + v_z * rz;
    const double r_dot_gradw = w_x * rx + w_y * ry + w_z * rz;

    const double tau_dc = get_DC( dxi_dx, u_prime, v_prime, w_prime );

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
          - NA * rho0 * f_body.x() );

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
          - NA * rho0 * f_body.y() );

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
          - NA * rho0 * f_body.z() );

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

  const double two_mu = 2.0 * vis_mu;

  const double curr = 0.0;

  Zero_Tangent_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    double u = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v = 0.0, v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w = 0.0, w_x = 0.0, w_y = 0.0, w_z = 0.0;
    double p = 0.0;
    
    Vector_3 coor( 0.0, 0.0, 0.0 );

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

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const Vector_3 f_body = get_f(coor, curr);

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      Residual[4*A+1] += gwts * ( NA * rho0 * (u*u_x + v*u_y + w*u_z) 
          - NA_x * p
          + two_mu * NA_x * u_x
          + vis_mu * NA_y * (u_y + v_x)
          + vis_mu * NA_z * (u_z + w_x)
          - NA * rho0 * f_body.x() );

      Residual[4*A+2] += gwts * ( NA * rho0 * (u*v_x + v*v_y + w*v_z) 
          - NA_y * p
          + vis_mu * NA_x * (u_y + v_x)
          + two_mu * NA_y * v_y
          + vis_mu * NA_z * (v_z + w_y)
          - NA * rho0 * f_body.y() );

      Residual[4*A+3] += gwts * ( NA * rho0 * (u*w_x + v*w_y + w*w_z) 
          - NA_z * p
          + vis_mu * NA_x * (u_z + w_x)
          + vis_mu * NA_y * (w_y + v_z)
          + two_mu * NA_z * w_z
          - NA * rho0 * f_body.z() );

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

  const double curr = time + alpha_f * dt;

  Zero_sur_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);
  
    double surface_area;

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    Vector_3 coor( 0.0, 0.0, 0.0 );

    for(int ii=0; ii<snLocBas; ++ii)
    {
      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    const Vector_3 traction = get_ebc_fun( ebc_id, coor, curr, n_out );

    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual[4*A+1] -= surface_area * quad -> get_qw(qua) * R[A] * traction.x();
      sur_Residual[4*A+2] -= surface_area * quad -> get_qw(qua) * R[A] * traction.y();
      sur_Residual[4*A+3] -= surface_area * quad -> get_qw(qua) * R[A] * traction.z();
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

  Zero_sur_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);
  
    double surface_area;

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual[4*A+1] += surface_area * quad -> get_qw(qua) * R[A] * n_out.x() * val;
      sur_Residual[4*A+2] += surface_area * quad -> get_qw(qua) * R[A] * n_out.y() * val;
      sur_Residual[4*A+3] += surface_area * quad -> get_qw(qua) * R[A] * n_out.z() * val;
    }
  }
}

void PLocAssem_Tet_CMM_GenAlpha::Assem_Residual_BackFlowStab(
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  Zero_sur_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);
    
    double surface_area, factor;

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    double u = 0.0, v = 0.0, w = 0.0;;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      const int ii4 = ii * 4;
      u += sol[ii4+1] * R[ii];
      v += sol[ii4+2] * R[ii];
      w += sol[ii4+3] * R[ii];
    }

    const double temp = u * n_out.x() + v * n_out.y() + w * n_out.z();

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
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double dd_dv = alpha_f * gamma * dt;

  Zero_sur_Tangent_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);
    
    double surface_area, factor;

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    double u = 0.0, v = 0.0, w = 0.0;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      u += sol[ii*4+1] * R[ii];
      v += sol[ii*4+2] * R[ii];
      w += sol[ii*4+3] * R[ii];
    }

    const double temp = u * n_out.x() + v * n_out.y() + w * n_out.z();

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

  double flrate = 0.0;

  for(int qua =0; qua< face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);
  
    double surface_area;

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    double u = 0.0, v = 0.0, w = 0.0;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      u += sol[ii*4+1] * R[ii];
      v += sol[ii*4+2] * R[ii];
      w += sol[ii*4+3] * R[ii];
    }

    flrate += surface_area * quad->get_qw(qua) * ( u * n_out.x() + v * n_out.y() + w * n_out.z() );
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

  // Initialize the two variables to be passed out
  pres = 0.0; area = 0.0;

  for(int qua =0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    double pp = 0.0;
    for(int ii=0; ii<snLocBas; ++ii) pp += sol[4*ii+0] * R[ii];

    pres += element->get_detJac(qua) * quad->get_qw(qua) * pp;
    area += element->get_detJac(qua) * quad->get_qw(qua);
  }
}


void PLocAssem_Tet_CMM_GenAlpha::Assem_Residual_EBC_Wall(
    const double &time, const double &dt,
    const double * const &dot_sol,
    const double * const &sol,
    const double * const &sol_wall_disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &ele_thickness,
    const double * const &ele_youngsmod,
    const double * const &ele_springconst,
    const double * const &ele_dampingconst,
    const double * const &qua_prestress,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z ); 

  const double curr = time + alpha_f * dt;

  // Global Cauchy stress at all quadrature points
  std::vector<Tensor2_3D> sigma( face_nqp, Tensor2_3D() );
  get_Wall_CauchyStress(sol_wall_disp, element, ele_youngsmod, sigma );

  Zero_sur_Residual();

  for(int qua=0; qua<face_nqp; ++qua)
  {
    const std::vector<double> R = element -> get_R( qua );

    // For membrane elements, basis function gradients are computed
    // with respect to lamina coords
    const std::vector<double> dR_dxl = element -> get_dR_dx( qua );
    const std::vector<double> dR_dyl = element -> get_dR_dy( qua );

    // Global-to-local rotation matrix Q
    const Tensor2_3D Q = element->get_rotationMatrix(qua);

    double u_t = 0.0, v_t = 0.0, w_t = 0.0, u = 0.0, v = 0.0, w = 0.0;
    double disp_x = 0.0, disp_y = 0.0, disp_z = 0.0;
    double h_w = 0.0, ks_w = 0.0, cs_w = 0.0;

    Vector_3 coor( 0.0, 0.0, 0.0 );

    for(int ii=0; ii<snLocBas; ++ii)
    {
      u_t += dot_sol[ii*4+1] * R[ii];
      v_t += dot_sol[ii*4+2] * R[ii];
      w_t += dot_sol[ii*4+3] * R[ii];

      u += sol[ii*4+1] * R[ii];
      v += sol[ii*4+2] * R[ii];
      w += sol[ii*4+3] * R[ii];

      disp_x += sol_wall_disp[ii*3+0] * R[ii];
      disp_y += sol_wall_disp[ii*3+1] * R[ii];
      disp_z += sol_wall_disp[ii*3+2] * R[ii];

      h_w  += ele_thickness[ii] * R[ii];
      ks_w += ele_springconst[ii]  * R[ii];
      cs_w += ele_dampingconst[ii] * R[ii];

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    // Add prestress: convert from Voigt notation (comps 11, 22, 33, 23, 13, 12)
    sigma[qua].xx() += qua_prestress[qua*6];
    sigma[qua].xy() += qua_prestress[qua*6+5];
    sigma[qua].xz() += qua_prestress[qua*6+4];
    sigma[qua].yx() += qua_prestress[qua*6+5];
    sigma[qua].yy() += qua_prestress[qua*6+1];
    sigma[qua].yz() += qua_prestress[qua*6+3];
    sigma[qua].zx() += qua_prestress[qua*6+4];
    sigma[qua].zy() += qua_prestress[qua*6+3];
    sigma[qua].zz() += qua_prestress[qua*6+2];

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    // Body force acting on the wall
    const Vector_3 fw =  get_fw(coor, curr);

    // Basis function gradients with respect to global coords
    // dR/dx_{i} = Q_{ji} * dR/dxl_{j}. Note that dR/dzl = 0.0
    std::vector<double> dR_dx(snLocBas, 0.0), dR_dy(snLocBas, 0.0), dR_dz(snLocBas, 0.0);
    for(int ii=0; ii<snLocBas; ++ii)
    {
      dR_dx[ii] = Q.xx() * dR_dxl[ii] + Q.yx() * dR_dyl[ii];
      dR_dy[ii] = Q.xy() * dR_dxl[ii] + Q.yy() * dR_dyl[ii];
      dR_dz[ii] = Q.xz() * dR_dxl[ii] + Q.yz() * dR_dyl[ii];
    }

    for(int A=0; A<snLocBas; ++A)
    {
      const double NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      sur_Residual[4*A+1] += gwts * h_w * ( R[A] * rho_w * ( u_t - fw.x() )
          + NA_x * sigma[qua].xx() + NA_y * sigma[qua].xy() + NA_z * sigma[qua].xz() )
          + gwts * R[A] * ( ks_w * disp_x + cs_w * u ); 
      
      sur_Residual[4*A+2] += gwts * h_w * ( R[A] * rho_w * ( v_t - fw.y() )
          + NA_x * sigma[qua].yx() + NA_y * sigma[qua].yy() + NA_z * sigma[qua].yz() )
          + gwts * R[A] * ( ks_w * disp_y + cs_w * v ); 
      
      sur_Residual[4*A+3] += gwts * h_w * ( R[A] * rho_w * ( w_t - fw.z() )
          + NA_x * sigma[qua].zx() + NA_y * sigma[qua].zy() + NA_z * sigma[qua].zz() )
          + gwts * R[A] * ( ks_w * disp_z + cs_w * w ); 
    }

  } // end qua loop

}


void PLocAssem_Tet_CMM_GenAlpha::Assem_Tangent_Residual_EBC_Wall(
    const double &time, const double &dt,
    const double * const &dot_sol,
    const double * const &sol,
    const double * const &sol_wall_disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &ele_thickness,
    const double * const &ele_youngsmod,
    const double * const &ele_springconst,
    const double * const &ele_dampingconst,
    const double * const &qua_prestress,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z ); 

  const int dim = 3;

  const double dd_dv = alpha_f * gamma * dt;
  const double dd_du = dd_dv * dd_dv / alpha_m;

  const double curr = time + alpha_f * dt;

  // Global Cauchy stress at all quadrature points
  std::vector<Tensor2_3D> sigma( face_nqp, Tensor2_3D() );
  get_Wall_CauchyStress(sol_wall_disp, element, ele_youngsmod, sigma );

  Zero_sur_Tangent_Residual();

  for(int qua=0; qua<face_nqp; ++qua)
  {
    const std::vector<double> R = element -> get_R( qua );

    // For membrane elements, basis function gradients are computed
    // with respect to lamina coords
    const std::vector<double> dR_dxl = element -> get_dR_dx( qua );
    const std::vector<double> dR_dyl = element -> get_dR_dy( qua );

    // Lamina and global stiffness matrices
    double * Kl = new double [ (snLocBas*dim) * (snLocBas*dim) ] {};
    double * Kg = new double [ (snLocBas*dim) * (snLocBas*dim) ] {};

    // Global-to-local rotation matrix Q
    const Tensor2_3D Q = element->get_rotationMatrix(qua);

    double u_t = 0.0, v_t = 0.0, w_t = 0.0, u = 0.0, v = 0.0, w = 0.0; 
    double disp_x = 0.0, disp_y = 0.0, disp_z = 0.0;
    double h_w = 0.0, E_w = 0.0, ks_w = 0.0, cs_w = 0.0;

    Vector_3 coor( 0.0, 0.0, 0.0 );

    for(int ii=0; ii<snLocBas; ++ii)
    {
      u_t += dot_sol[ii*4+1] * R[ii];
      v_t += dot_sol[ii*4+2] * R[ii];
      w_t += dot_sol[ii*4+3] * R[ii];

      u += sol[ii*4+1] * R[ii];
      v += sol[ii*4+2] * R[ii];
      w += sol[ii*4+3] * R[ii];

      disp_x += sol_wall_disp[ii*3+0] * R[ii];
      disp_y += sol_wall_disp[ii*3+1] * R[ii];
      disp_z += sol_wall_disp[ii*3+2] * R[ii];

      h_w += ele_thickness[ii] * R[ii];
      E_w += ele_youngsmod[ii] * R[ii];

      ks_w += ele_springconst[ii]  * R[ii];
      cs_w += ele_dampingconst[ii] * R[ii];

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    // Body force acting on the wall
    const Vector_3 fw = get_fw(coor, curr);

    const double coef = E_w / (1.0 - nu_w * nu_w);

    // Add prestress: convert from Voigt notation (comps 11, 22, 33, 23, 13, 12)
    sigma[qua].xx() += qua_prestress[qua*6];
    sigma[qua].xy() += qua_prestress[qua*6+5];
    sigma[qua].xz() += qua_prestress[qua*6+4];
    sigma[qua].yx() += qua_prestress[qua*6+5];
    sigma[qua].yy() += qua_prestress[qua*6+1];
    sigma[qua].yz() += qua_prestress[qua*6+3];
    sigma[qua].zx() += qua_prestress[qua*6+4];
    sigma[qua].zy() += qua_prestress[qua*6+3];
    sigma[qua].zz() += qua_prestress[qua*6+2];

    // Basis function gradients with respect to global coords
    // dR/dx_{i} = Q_{ji} * dR/dxl_{j}. Note that dR/dzl = 0.0
    std::vector<double> dR_dx(snLocBas, 0.0), dR_dy(snLocBas, 0.0), dR_dz(snLocBas, 0.0);
    for(int ii=0; ii<snLocBas; ++ii)
    {
      dR_dx[ii] = Q.xx() * dR_dxl[ii] + Q.yx() * dR_dyl[ii];
      dR_dy[ii] = Q.xy() * dR_dxl[ii] + Q.yy() * dR_dyl[ii];
      dR_dz[ii] = Q.xz() * dR_dxl[ii] + Q.yz() * dR_dyl[ii];
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

      sur_Residual[4*A+1] += gwts * h_w * ( R[A] * rho_w * ( u_t - fw.x() )
          + NA_x * sigma[qua].xx() + NA_y * sigma[qua].xy() + NA_z * sigma[qua].xz() )
          + gwts * R[A] * ( ks_w * disp_x + cs_w * u ); 

      sur_Residual[4*A+2] += gwts * h_w * ( R[A] * rho_w * ( v_t - fw.y() )
          + NA_x * sigma[qua].yx() + NA_y * sigma[qua].yy() + NA_z * sigma[qua].yz() )
          + gwts * R[A] * ( ks_w * disp_y + cs_w * v ); 
      
      sur_Residual[4*A+3] += gwts * h_w * ( R[A] * rho_w * ( w_t - fw.z() )
          + NA_x * sigma[qua].zx() + NA_y * sigma[qua].zy() + NA_z * sigma[qua].zz() )
          + gwts * R[A] * ( ks_w * disp_z + cs_w * w ); 

      for(int B=0; B<snLocBas; ++B)
      {
        // Momentum-x with respect to u, v, w
        sur_Tangent[ 4*snLocBas*(4*A+1) + 4*B+1 ] += gwts * h_w * (
            alpha_m * rho_w * R[A] * R[B]
            + dd_du * Kg[ (snLocBas*dim)*(A*dim) + (B*dim) ] )
            + gwts * R[A] * R[B] * ( dd_du * ks_w + dd_dv * cs_w );

        sur_Tangent[ 4*snLocBas*(4*A+1) + 4*B+2 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim) + (B*dim+1) ] );

        sur_Tangent[ 4*snLocBas*(4*A+1) + 4*B+3 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim) + (B*dim+2) ] );

        // Momentum-y with respect to u, v, w
        sur_Tangent[ 4*snLocBas*(4*A+2) + 4*B+1 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim+1) + (B*dim) ] );

        sur_Tangent[ 4*snLocBas*(4*A+2) + 4*B+2 ] += gwts * h_w * (
            alpha_m * rho_w * R[A] * R[B]
            + dd_du * Kg[ (snLocBas*dim)*(A*dim+1) + (B*dim+1) ] )
            + gwts * R[A] * R[B] * ( dd_du * ks_w + dd_dv * cs_w );

        sur_Tangent[ 4*snLocBas*(4*A+2) + 4*B+3 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim+1) + (B*dim+2) ] );

        // Momentum-z with respect to u, v, w
        sur_Tangent[ 4*snLocBas*(4*A+3) + 4*B+1 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim+2) + (B*dim) ] );

        sur_Tangent[ 4*snLocBas*(4*A+3) + 4*B+2 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim+2) + (B*dim+1) ] );

        sur_Tangent[ 4*snLocBas*(4*A+3) + 4*B+3 ] += gwts * h_w * (
            alpha_m * rho_w * R[A] * R[B]
            + dd_du * Kg[ (snLocBas*dim)*(A*dim+2) + (B*dim+2) ] )
            + gwts * R[A] * R[B] * ( dd_du * ks_w + dd_dv * cs_w );

      } // end B loop
    } // end A loop

    delete [] Kl; delete [] Kg;
    Kl = nullptr; Kg = nullptr;

  } // end qua loop

}


void PLocAssem_Tet_CMM_GenAlpha::get_Wall_CauchyStress(
    const double * const &sol_wall_disp,
    const FEAElement * const &element,
    const double * const &ele_youngsmod,
    std::vector<Tensor2_3D> &sigma )
{
  SYS_T::print_fatal_if( element -> get_numQuapts() != face_nqp, "Error: The element's data structure is incompatible with the face_nqp in the local assembly.\n" );

  const int dim = 3;

  // Lamina displacements
  double * sol_wall_disp_l = new double [ snLocBas * dim ];

  for(int qua=0; qua<face_nqp; ++qua)
  {
    const std::vector<double> R = element -> get_R( qua );

    // For membrane elements, basis function gradients are computed
    // with respect to lamina coords
    const std::vector<double> dR_dxl = element -> get_dR_dx( qua );
    const std::vector<double> dR_dyl = element -> get_dR_dy( qua );

    // Global-to-local rotation matrix Q
    const Tensor2_3D Q = element->get_rotationMatrix(qua);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      sol_wall_disp_l[dim*ii]   = sol_wall_disp[dim*ii] * Q.xx() 
        + sol_wall_disp[dim*ii+1] * Q.xy() + sol_wall_disp[dim*ii+2] * Q.xz();

      sol_wall_disp_l[dim*ii+1] = sol_wall_disp[dim*ii] * Q.yx()
        + sol_wall_disp[dim*ii+1] * Q.yy() + sol_wall_disp[dim*ii+2] * Q.yz();

      sol_wall_disp_l[dim*ii+2] = sol_wall_disp[dim*ii] * Q.zx()
        + sol_wall_disp[dim*ii+1] * Q.zy() + sol_wall_disp[dim*ii+2] * Q.zz();
    }

    double E_w = 0.0;
    double u1l_xl = 0.0, u2l_xl = 0.0, u3l_xl = 0.0;
    double u1l_yl = 0.0, u2l_yl = 0.0, u3l_yl = 0.0;

    for(int ii=0; ii<snLocBas; ++ii)
    {
      E_w += ele_youngsmod[ii] * R[ii];

      u1l_xl += sol_wall_disp_l[dim*ii]   * dR_dxl[ii];
      u1l_yl += sol_wall_disp_l[dim*ii]   * dR_dyl[ii];
      u2l_xl += sol_wall_disp_l[dim*ii+1] * dR_dxl[ii];
      u2l_yl += sol_wall_disp_l[dim*ii+1] * dR_dyl[ii];
      u3l_xl += sol_wall_disp_l[dim*ii+2] * dR_dxl[ii];
      u3l_yl += sol_wall_disp_l[dim*ii+2] * dR_dyl[ii];
    }

    const double coef = E_w / (1.0 - nu_w * nu_w);

    // Lamina Cauchy stress
    sigma[qua] = Tensor2_3D(
        u1l_xl + nu_w*u2l_yl,               0.5*(1.0-nu_w) * (u1l_yl + u2l_xl), 0.5*kappa_w*(1.0-nu_w) * u3l_xl,
        0.5*(1.0-nu_w) * (u1l_yl + u2l_xl), nu_w*u1l_xl + u2l_yl,               0.5*kappa_w*(1.0-nu_w) * u3l_yl,
        0.5*kappa_w*(1.0-nu_w) * u3l_xl,    0.5*kappa_w*(1.0-nu_w) * u3l_yl,    0.0 );
    sigma[qua] *= coef;

    // Global Cauchy stress: Q^T * lamina_Cauchy * Q
    sigma[qua].MatRot(Q);
  }

  delete [] sol_wall_disp_l; sol_wall_disp_l = nullptr;
}

// EOF
