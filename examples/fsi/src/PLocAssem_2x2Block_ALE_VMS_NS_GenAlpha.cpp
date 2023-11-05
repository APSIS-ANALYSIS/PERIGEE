#include "PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha.hpp"

PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha(
    const TimeMethod_GenAlpha * const &tm_gAlpha,
    const int &in_nlocbas, const int &in_snlocbas,
    const double &in_rho, const double &in_vis_mu, const double &in_beta )
: rho0( in_rho ), vis_mu( in_vis_mu ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), beta(in_beta), CI(36.0), CT(4.0),
  nLocBas(in_nlocbas), snLocBas(in_snlocbas),
  vec_size_0( nLocBas * 3 ), vec_size_1( nLocBas ), 
  sur_size_0( snLocBas * 3 )
{
  Tangent00 = new PetscScalar[vec_size_0 * vec_size_0];
  Tangent01 = new PetscScalar[vec_size_0 * vec_size_1];
  Tangent10 = new PetscScalar[vec_size_1 * vec_size_0];
  Tangent11 = new PetscScalar[vec_size_1 * vec_size_1];

  Residual0 = new PetscScalar[vec_size_0];
  Residual1 = new PetscScalar[vec_size_1];

  sur_Tangent00 = new PetscScalar[sur_size_0 * sur_size_0];
  sur_Residual0 = new PetscScalar[sur_size_0];

  Zero_Tangent_Residual();
  Zero_sur_Tangent_Residual();

  // print info of this assembly routine
  print_info();
}

PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::~PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha()
{
  delete [] Tangent00; Tangent00 = nullptr;
  delete [] Tangent01; Tangent01 = nullptr;
  delete [] Tangent10; Tangent10 = nullptr;
  delete [] Tangent11; Tangent11 = nullptr;

  delete [] Residual0; Residual0 = nullptr;
  delete [] Residual1; Residual1 = nullptr;

  delete [] sur_Tangent00; sur_Tangent00 = nullptr;

  delete [] sur_Residual0; sur_Residual0 = nullptr;
}


void PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint("  Three-dimensional Incompressible Navier-Stokes equations: \n");
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
  SYS_T::commPrint("  Input solution vector should have 7 dofs including the mesh motion\n      information in the first three slots. \n");
  SYS_T::commPrint("  Density rho is in front of the du/dt and the dimension of the\n      equations is momentum time rate. \n");
  SYS_T::print_sep_line();
}

SymmTensor2_3D PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::get_metric(
    const std::array<double, 9> &f ) const
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

std::array<double, 2> PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::get_tau(
    const double &dt, const std::array<double, 9> &dxi_dx,
    const double &u, const double &v, const double &w ) const
{
  // Use K matrix to correct the metric
  const SymmTensor2_3D G = get_metric( dxi_dx );

  const Vector_3 velo_vec( u, v, w );

  const double temp_nu = vis_mu / rho0;

  const double denom_m = std::sqrt(CT / (dt*dt) 
  + G.VecMatVec( velo_vec, velo_vec ) 
  + CI * temp_nu * temp_nu * G.MatContraction( G ));

  return {{1.0 / (rho0 * denom_m ), rho0 * denom_m / G.tr()}};
}

double PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::get_DC(
    const std::array<double, 9> &dxi_dx,
    const double &u, const double &v, const double &w ) const
{
  // const SymmTensor2_3D G = get_metric( dxi_dx );
  // const Vector_3 velo_vec( u, v, w );
  // double dc_tau = G.VecMatVec( velo_vec, velo_vec );

  //if(dc_tau > 1.0e-15) dc_tau = rho0 * std::pow(dc_tau, -0.5);
  //else dc_tau = 0.0;
  const double dc_tau = 0.0;

  return dc_tau;
}

void PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::Assem_Residual(
    const double &time, const double &dt,
    const double * const &dot_disp,
    const double * const &dot_velo,
    const double * const &dot_pres,
    const double * const &disp,
    const double * const &velo,
    const double * const &pres,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  double R[nLocBas], dR_dx[nLocBas], dR_dy[nLocBas], dR_dz[nLocBas];
  double curPt_x[nLocBas], curPt_y[nLocBas], curPt_z[nLocBas];

  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp, nLocBas, curPt_x, curPt_y, curPt_z);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const double two_mu = 2.0 * vis_mu;

  const double curr = time + alpha_f * dt;

  Zero_Residual();

  const int nqp = quad -> get_num_quadPts();

  for(int qua=0; qua<nqp; ++qua)
  {
    double u = 0.0, u_t = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v = 0.0, v_t = 0.0, v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w = 0.0, w_t = 0.0, w_x = 0.0, w_y = 0.0, w_z = 0.0;
    double p = 0.0, p_x = 0.0, p_y = 0.0, p_z = 0.0;
    double mu = 0.0, mv = 0.0, mw = 0.0; // mesh velocity, i.e. hat-v
    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_R_gradR( qua, R, dR_dx, dR_dy, dR_dz );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      mu += dot_disp[ii*3+0] * R[ii];
      mv += dot_disp[ii*3+1] * R[ii];
      mw += dot_disp[ii*3+2] * R[ii];

      u_t += dot_velo[ii*3+0] * R[ii];
      v_t += dot_velo[ii*3+1] * R[ii];
      w_t += dot_velo[ii*3+2] * R[ii];

      u += velo[ii*3+0] * R[ii];
      v += velo[ii*3+1] * R[ii];
      w += velo[ii*3+2] * R[ii];
      p += pres[ii]     * R[ii];

      u_x += velo[ii*3+0] * dR_dx[ii];
      v_x += velo[ii*3+1] * dR_dx[ii];
      w_x += velo[ii*3+2] * dR_dx[ii];
      p_x += pres[ii]     * dR_dx[ii];

      u_y += velo[ii*3+0] * dR_dy[ii];
      v_y += velo[ii*3+1] * dR_dy[ii];
      w_y += velo[ii*3+2] * dR_dy[ii];
      p_y += pres[ii]     * dR_dy[ii];

      u_z += velo[ii*3+0] * dR_dz[ii];
      v_z += velo[ii*3+1] * dR_dz[ii];
      w_z += velo[ii*3+2] * dR_dz[ii];
      p_z += pres[ii]     * dR_dz[ii];

      coor.x() += curPt_x[ii] * R[ii];
      coor.y() += curPt_y[ii] * R[ii];
      coor.z() += curPt_z[ii] * R[ii];
    }

    // v - hat(v)
    const double cu = u - mu;
    const double cv = v - mv;
    const double cw = w - mw;

    const auto dxi_dx = element->get_invJacobian(qua);
    const std::array<double, 2> tau = get_tau(dt, dxi_dx, cu, cv, cw);
    const double tau_m = tau[0];
    const double tau_c = tau[1];

    const double tau_m_2 = tau_m * tau_m;

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const Vector_3 f_body = get_f(coor, curr);

    const double rx = rho0 * ( u_t + u_x * cu + u_y * cv + u_z * cw - f_body.x() ) + p_x;
    const double ry = rho0 * ( v_t + v_x * cu + v_y * cv + v_z * cw - f_body.y() ) + p_y;
    const double rz = rho0 * ( w_t + w_x * cu + w_y * cv + w_z * cw - f_body.z() ) + p_z;

    const double div_vel = u_x + v_y + w_z;

    const double u_prime = -1.0 * tau_m * rx;
    const double v_prime = -1.0 * tau_m * ry;
    const double w_prime = -1.0 * tau_m * rz;

    const double tau_dc = get_DC( dxi_dx, u_prime, v_prime, w_prime );

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      const double velo_dot_gradR = NA_x * cu + NA_y * cv + NA_z * cw;
      const double r_dot_gradR = NA_x * rx + NA_y * ry + NA_z * rz;
      const double r_dot_gradu = u_x * rx + u_y * ry + u_z * rz;
      const double r_dot_gradv = v_x * rx + v_y * ry + v_z * rz;
      const double r_dot_gradw = w_x * rx + w_y * ry + w_z * rz;
      const double velo_prime_dot_gradR = NA_x * u_prime + NA_y * v_prime + NA_z * w_prime;

      Residual1[A] += gwts * ( NA * div_vel + tau_m * r_dot_gradR );

      Residual0[3*A+0] += gwts * ( NA * rho0 * u_t
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
          - NA * rho0 * f_body.x() );

      Residual0[3*A+1] += gwts * ( NA * rho0 * v_t
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
          - NA * rho0 * f_body.y() );

      Residual0[3*A+2] += gwts * (NA * rho0 * w_t
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
          - NA * rho0 * f_body.z() );
    }
  }
}

void PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &dot_disp,
    const double * const &dot_velo,
    const double * const &dot_pres,
    const double * const &disp,
    const double * const &velo,
    const double * const &pres,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  double R[nLocBas], dR_dx[nLocBas], dR_dy[nLocBas], dR_dz[nLocBas];
  double curPt_x[nLocBas], curPt_y[nLocBas], curPt_z[nLocBas];

  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp, nLocBas, curPt_x, curPt_y, curPt_z);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const double two_mu = 2.0 * vis_mu;

  const double rho0_2 = rho0 * rho0;

  const double curr = time + alpha_f * dt;

  const double dd_dv = alpha_f * gamma * dt;

  Zero_Tangent_Residual();

  const int nqp = quad -> get_num_quadPts();

  for(int qua=0; qua<nqp; ++qua)
  {
    double u = 0.0, u_t = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v = 0.0, v_t = 0.0, v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w = 0.0, w_t = 0.0, w_x = 0.0, w_y = 0.0, w_z = 0.0;
    double p = 0.0, p_x = 0.0, p_y = 0.0, p_z = 0.0;
    double mu = 0.0, mv = 0.0, mw = 0.0;
    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_R_gradR( qua, R, dR_dx, dR_dy, dR_dz );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      mu += dot_disp[ii*3+0] * R[ii];
      mv += dot_disp[ii*3+1] * R[ii];
      mw += dot_disp[ii*3+2] * R[ii];

      u_t += dot_velo[ii*3+0] * R[ii];
      v_t += dot_velo[ii*3+1] * R[ii];
      w_t += dot_velo[ii*3+2] * R[ii];

      u += velo[ii*3+0] * R[ii];
      v += velo[ii*3+1] * R[ii];
      w += velo[ii*3+2] * R[ii];
      p += pres[ii]     * R[ii];

      u_x += velo[ii*3+0] * dR_dx[ii];
      v_x += velo[ii*3+1] * dR_dx[ii];
      w_x += velo[ii*3+2] * dR_dx[ii];
      p_x += pres[ii]     * dR_dx[ii];

      u_y += velo[ii*3+0] * dR_dy[ii];
      v_y += velo[ii*3+1] * dR_dy[ii];
      w_y += velo[ii*3+2] * dR_dy[ii];
      p_y += pres[ii]     * dR_dy[ii];

      u_z += velo[ii*3+0] * dR_dz[ii];
      v_z += velo[ii*3+1] * dR_dz[ii];
      w_z += velo[ii*3+2] * dR_dz[ii];
      p_z += pres[ii]     * dR_dz[ii];

      coor.x() += curPt_x[ii] * R[ii];
      coor.y() += curPt_y[ii] * R[ii];
      coor.z() += curPt_z[ii] * R[ii];
    }

    const double cu = u - mu;
    const double cv = v - mv;
    const double cw = w - mw;

    const auto dxi_dx = element->get_invJacobian(qua);
    const std::array<double, 2> tau = get_tau(dt, dxi_dx, cu, cv, cw);
    const double tau_m = tau[0];
    const double tau_c = tau[1];

    const double tau_m_2 = tau_m * tau_m;

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua); 

    const Vector_3 f_body = get_f(coor, curr);

    const double rx = rho0 * ( u_t + u_x * cu + u_y * cv + u_z * cw - f_body.x() ) + p_x;
    const double ry = rho0 * ( v_t + v_x * cu + v_y * cv + v_z * cw - f_body.y() ) + p_y;
    const double rz = rho0 * ( w_t + w_x * cu + w_y * cv + w_z * cw - f_body.z() ) + p_z;

    const double div_vel = u_x + v_y + w_z;

    const double u_prime = -1.0 * tau_m * rx;
    const double v_prime = -1.0 * tau_m * ry;
    const double w_prime = -1.0 * tau_m * rz;

    const double tau_dc = get_DC( dxi_dx, u_prime, v_prime, w_prime );

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      const double velo_dot_gradR = NA_x * cu + NA_y * cv + NA_z * cw;
      const double r_dot_gradR = NA_x * rx + NA_y * ry + NA_z * rz;
      const double r_dot_gradu = u_x * rx + u_y * ry + u_z * rz;
      const double r_dot_gradv = v_x * rx + v_y * ry + v_z * rz;
      const double r_dot_gradw = w_x * rx + w_y * ry + w_z * rz;
      const double velo_prime_dot_gradR = NA_x * u_prime + NA_y * v_prime + NA_z * w_prime;

      Residual1[A] += gwts * ( NA * div_vel + tau_m * r_dot_gradR );

      Residual0[3*A+0] += gwts * ( NA * rho0 * u_t
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
          - NA * rho0 * f_body.x() );

      Residual0[3*A+1] += gwts * ( NA * rho0 * v_t
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
          - NA * rho0 * f_body.y() );

      Residual0[3*A+2] += gwts * (NA * rho0 * w_t
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
          - NA * rho0 * f_body.z() );

      for(int B=0; B<nLocBas; ++B)
      {
        const double NB = R[B], NB_x = dR_dx[B], NB_y = dR_dy[B], NB_z = dR_dz[B];
        const double velo_dot_gradNB = cu * NB_x + cv * NB_y + cw * NB_z;
        const double velo_prime_dot_gradNB = u_prime * NB_x + v_prime * NB_y + w_prime * NB_z;

        const double NANB  = NA*NB, NANBx = NA*NB_x, NANBy = NA*NB_y, NANBz = NA*NB_z;
        const double NAxNB = NA_x*NB, NAxNBx = NA_x*NB_x, NAxNBy = NA_x*NB_y, NAxNBz = NA_x*NB_z;
        const double NAyNB = NA_y*NB, NAyNBx = NA_y*NB_x, NAyNBy = NA_y*NB_y, NAyNBz = NA_y*NB_z;
        const double NAzNB = NA_z*NB, NAzNBx = NA_z*NB_x, NAzNBy = NA_z*NB_y, NAzNBz = NA_z*NB_z;

        const double drx_du_B = rho0 * ( u_x * NB + velo_dot_gradNB ); 
        const double drx_dv_B = rho0 * u_y * NB;
        const double drx_dw_B = rho0 * u_z * NB;

        const double dry_du_B = rho0 * v_x * NB;
        const double dry_dv_B = rho0 * ( v_y * NB + velo_dot_gradNB );
        const double dry_dw_B = rho0 * v_z * NB;

        const double drz_du_B = rho0 * w_x * NB;
        const double drz_dv_B = rho0 * w_y * NB;
        const double drz_dw_B = rho0 * (w_z * NB + velo_dot_gradNB);

        // Continuity equation with respect to p, u, v, w
        Tangent11[nLocBas*A+B] += gwts * dd_dv * tau_m * (NAxNBx + NAyNBy + NAzNBz);

        Tangent10[3*nLocBas*A+3*B+0] += gwts * ( alpha_m * tau_m * rho0 * NAxNB
            + dd_dv * ( NANBx + tau_m * NA_x * drx_du_B
              + tau_m * NA_y * dry_du_B + tau_m * NA_z * drz_du_B ) );

        Tangent10[3*nLocBas*A+3*B+1] += gwts * ( alpha_m * tau_m * rho0 * NAyNB
            + dd_dv * ( NANBy + tau_m * NA_x * drx_dv_B
              + tau_m * NA_y * dry_dv_B + tau_m * NA_z * drz_dv_B ) );

        Tangent10[3*nLocBas*A+3*B+2] += gwts * ( alpha_m * tau_m * rho0 * NAzNB
            + dd_dv * ( NANBz + tau_m * NA_x * drx_dw_B
              + tau_m * NA_y * dry_dw_B + tau_m * NA_z * drz_dw_B ) );

        // Momentum-x with respect to p, u, v, w
        Tangent01[nLocBas*(3*A+0)+B] += gwts * dd_dv * ((-1.0) * NAxNB
            + velo_dot_gradR * tau_m * rho0 * NB_x
            - NA * tau_m * rho0 * (u_x * NB_x + u_y * NB_y + u_z * NB_z)
            - 2.0 * tau_m_2 * rho0 * rx * NAxNBx
            - tau_m_2 * rho0 * NA_y * (rx * NB_y + ry * NB_x)
            - tau_m_2 * rho0 * NA_z * (rx * NB_z + rz * NB_x) );

        Tangent00[3*nLocBas*(3*A+0)+3*B+0] += gwts * ( 
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

        Tangent00[3*nLocBas*(3*A+0)+3*B+1] += gwts * ( 
            alpha_m * (-1.0) * rho0_2 * (tau_m * u_y * NANB + tau_m_2 * rx * NAyNB)
            + dd_dv * ( NANB * rho0 * u_y + vis_mu * NAyNBx 
              + rho0 * tau_m * rx * NAyNB
              + velo_dot_gradR * rho0 * tau_m * drx_dv_B
              - rho0 * tau_m * NA * (u_x*drx_dv_B + u_y*dry_dv_B + u_z*drz_dv_B)
              + tau_c * NAxNBy
              - 2.0 * rho0 * tau_m_2 * rx * NA_x * drx_dv_B
              - rho0 * tau_m_2 * NA_y * (rx * dry_dv_B + ry * drx_dv_B)
              - rho0 * tau_m_2 * NA_z * (rx * drz_dv_B + rz * drx_dv_B) ) );

        Tangent00[3*nLocBas*(3*A+0)+3*B+2] += gwts * (
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
        Tangent01[nLocBas*(3*A+1)+B] += gwts * dd_dv * ( (-1.0) * NAyNB
            + velo_dot_gradR * tau_m * rho0 * NB_y
            - NA * tau_m * rho0 * (v_x * NB_x + v_y * NB_y + v_z * NB_z)
            - tau_m_2 * rho0 * NA_x * (rx * NB_y + ry * NB_x)
            - 2.0 * tau_m_2 * rho0 * ry * NAyNBy
            - tau_m_2 * rho0 * NA_z * (ry * NB_z + rz * NB_y) );

        Tangent00[3*nLocBas*(3*A+1)+3*B+0] += gwts * (
            alpha_m * (-1.0) * rho0_2 * (tau_m * v_x * NANB + tau_m_2 * ry * NAxNB)
            + dd_dv * ( NANB * rho0 * v_x + vis_mu * NAxNBy
              + rho0 * tau_m * ry * NAxNB
              + velo_dot_gradR * rho0 * tau_m * dry_du_B
              - rho0 * tau_m * NA * (v_x*drx_du_B + v_y*dry_du_B + v_z*drz_du_B)
              + tau_c * NAyNBx
              - rho0 * tau_m_2 * NA_x * (ry * drx_du_B + rx * dry_du_B)
              - 2.0 * rho0 * tau_m_2 * ry * NA_y * dry_du_B
              - rho0 * tau_m_2 * NA_z * (ry * drz_du_B + rz * dry_du_B) ) );

        Tangent00[3*nLocBas*(3*A+1)+3*B+1] += gwts * (
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

        Tangent00[3*nLocBas*(3*A+1)+3*B+2] += gwts * (
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
        Tangent01[nLocBas*(3*A+2)+B] += gwts * dd_dv * ( (-1.0) * NAzNB
            + velo_dot_gradR * tau_m * rho0 * NB_z
            - NA * tau_m * rho0 * (w_x * NB_x + w_y * NB_y + w_z * NB_z)
            - tau_m_2 * rho0 * NA_x * (rx * NB_z + rz * NB_x)
            - tau_m_2 * rho0 * NA_y * (ry * NB_z + rz * NB_y)
            - 2.0 * tau_m_2 * rho0 * rz * NAzNBz );

        Tangent00[3*nLocBas*(3*A+2)+3*B+0] += gwts * (
            alpha_m * (-1.0) * rho0_2 * (tau_m * w_x * NANB + tau_m_2 * rz * NAxNB)
            + dd_dv * ( NANB * rho0 * w_x + vis_mu * NAxNBz
              + rho0 * tau_m * rz * NAxNB
              + velo_dot_gradR * rho0 * tau_m * drz_du_B
              - rho0 * tau_m * NA * (w_x*drx_du_B + w_y*dry_du_B + w_z*drz_du_B)
              + tau_c * NAzNBx
              - rho0 * tau_m_2 * NA_x * (rx * drz_du_B + rz * drx_du_B)
              - rho0 * tau_m_2 * NA_y * (ry * drz_du_B + rz * dry_du_B)
              - 2.0 * rho0 * tau_m_2 * rz * NA_z * drz_du_B ) );

        Tangent00[3*nLocBas*(3*A+2)+3*B+1] += gwts * (
            alpha_m * (-1.0) * rho0_2 * (tau_m * w_y * NANB + tau_m_2 * rz * NAyNB)
            + dd_dv * ( NANB * rho0 * w_y + vis_mu * NAyNBz
              + rho0 * tau_m * rz * NAyNB
              + velo_dot_gradR * rho0 * tau_m * drz_dv_B
              - rho0 * tau_m * NA * (w_x*drx_dv_B + w_y*dry_dv_B + w_z*drz_dv_B)
              + tau_c * NAzNBy
              - rho0 * tau_m_2 * NA_x * (rx * drz_dv_B + rz * drx_dv_B)
              - rho0 * tau_m_2 * NA_y * (ry * drz_dv_B + rz * dry_dv_B)
              - 2.0 * rho0 * tau_m_2 * rz * NA_z * drz_dv_B ) );

        Tangent00[3*nLocBas*(3*A+2)+3*B+2] += gwts * (
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
}

void PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::Assem_Mass_Residual(
    const double * const &disp,
    const double * const &velo,
    const double * const &pres,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  double R[nLocBas], dR_dx[nLocBas], dR_dy[nLocBas], dR_dz[nLocBas];
  double curPt_x[nLocBas], curPt_y[nLocBas], curPt_z[nLocBas];

  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp, nLocBas, curPt_x, curPt_y, curPt_z);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const double two_mu = 2.0 * vis_mu;

  const double curr = 0.0;

  Zero_Tangent_Residual();

  const int nqp = quad -> get_num_quadPts();

  for(int qua=0; qua<nqp; ++qua)
  {
    double u = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v = 0.0, v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w = 0.0, w_x = 0.0, w_y = 0.0, w_z = 0.0;
    double p = 0.0;
    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_R_gradR( qua, R, dR_dx, dR_dy, dR_dz );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      u += velo[ii*3+0] * R[ii];
      v += velo[ii*3+1] * R[ii];
      w += velo[ii*3+2] * R[ii];
      p += pres[ii]     * R[ii];

      u_x += velo[ii*3+0] * dR_dx[ii];
      v_x += velo[ii*3+1] * dR_dx[ii];
      w_x += velo[ii*3+2] * dR_dx[ii];

      u_y += velo[ii*3+0] * dR_dy[ii];
      v_y += velo[ii*3+1] * dR_dy[ii];
      w_y += velo[ii*3+2] * dR_dy[ii];

      u_z += velo[ii*3+0] * dR_dz[ii];
      v_z += velo[ii*3+1] * dR_dz[ii];
      w_z += velo[ii*3+2] * dR_dz[ii];

      coor.x() += curPt_x[ii] * R[ii];
      coor.y() += curPt_y[ii] * R[ii];
      coor.z() += curPt_z[ii] * R[ii];
    }

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const Vector_3 f_body = get_f(coor, curr);

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      Residual0[3*A+0] += gwts * ( NA * rho0 * (u*u_x + v*u_y + w*u_z) 
          - NA_x * p
          + two_mu * NA_x * u_x
          + vis_mu * NA_y * (u_y + v_x)
          + vis_mu * NA_z * (u_z + w_x)
          - NA * rho0 * f_body.x() );

      Residual0[3*A+1] += gwts * ( NA * rho0 * (u*v_x + v*v_y + w*v_z) 
          - NA_y * p
          + vis_mu * NA_x * (u_y + v_x)
          + two_mu * NA_y * v_y
          + vis_mu * NA_z * (v_z + w_y)
          - NA * rho0 * f_body.y() );

      Residual0[3*A+2] += gwts * ( NA * rho0 * (u*w_x + v*w_y + w*w_z) 
          - NA_z * p
          + vis_mu * NA_x * (u_z + w_x)
          + vis_mu * NA_y * (w_y + v_z)
          + two_mu * NA_z * w_z
          - NA * rho0 * f_body.z() );

      for(int B=0; B<nLocBas; ++B)
      {
        Tangent11[nLocBas*A + B] += gwts * rho0 * NA * R[B];

        Tangent00[3*nLocBas*(3*A+0) + 3*B+0] += gwts * rho0 * NA * R[B];
        Tangent00[3*nLocBas*(3*A+1) + 3*B+1] += gwts * rho0 * NA * R[B];
        Tangent00[3*nLocBas*(3*A+2) + 3*B+2] += gwts * rho0 * NA * R[B];
      }
    }
  }
}

void PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::Assem_Residual_EBC(
    const int &ebc_id,
    const double &time, const double &dt,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  double curPt_x[snLocBas], curPt_y[snLocBas], curPt_z[snLocBas];

  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp, snLocBas, curPt_x, curPt_y, curPt_z);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const int face_nqp = quad -> get_num_quadPts();

  double surface_area;

  const double curr = time + alpha_f * dt;

  Zero_sur_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);
    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    Vector_3 coor(0.0, 0.0, 0.0);
    for(int ii=0; ii<snLocBas; ++ii)
    {
      coor.x() += curPt_x[ii] * R[ii];
      coor.y() += curPt_y[ii] * R[ii];
      coor.z() += curPt_z[ii] * R[ii];
    }

    const Vector_3 traction = get_ebc_fun( ebc_id, coor, curr, n_out );

    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual0[3*A+0] -= surface_area * quad -> get_qw(qua) * R[A] * traction.x();
      sur_Residual0[3*A+1] -= surface_area * quad -> get_qw(qua) * R[A] * traction.y();
      sur_Residual0[3*A+2] -= surface_area * quad -> get_qw(qua) * R[A] * traction.z();
    }
  }
}

double PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::get_flowrate( 
    const double * const &disp,
    const double * const &velo,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  double curPt_x[snLocBas], curPt_y[snLocBas], curPt_z[snLocBas];

  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp, snLocBas, curPt_x, curPt_y, curPt_z);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const int face_nqp = quad -> get_num_quadPts();

  double flrate = 0.0;

  for(int qua =0; qua< face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    double surface_area;

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    double u = 0.0, v = 0.0, w = 0.0;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      u += velo[3*ii+0] * R[ii];
      v += velo[3*ii+1] * R[ii];
      w += velo[3*ii+2] * R[ii];
    }

    flrate += surface_area * quad->get_qw(qua) * ( u * n_out.x() + v * n_out.y() + w * n_out.z() );
  }

  return flrate;
}

void PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::get_pressure_area( 
    const double * const &disp,
    const double * const &pres,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad,
    double &pressure, double &area )
{
  double curPt_x[snLocBas], curPt_y[snLocBas], curPt_z[snLocBas];

  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp, snLocBas, curPt_x, curPt_y, curPt_z);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const int face_nqp = quad -> get_num_quadPts();

  // Initialize the two variables to be passed out
  pressure = 0.0; area = 0.0;

  for(int qua =0; qua< face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    double pp = 0.0;
    for(int ii=0; ii<snLocBas; ++ii) pp += pres[ii] * R[ii];

    pressure += element->get_detJac(qua) * quad->get_qw(qua) * pp;
    area     += element->get_detJac(qua) * quad->get_qw(qua);
  }
}

void PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::Assem_Residual_EBC_Resistance(
    const double &val,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  double curPt_x[snLocBas], curPt_y[snLocBas], curPt_z[snLocBas];

  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp, snLocBas, curPt_x, curPt_y, curPt_z);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const int face_nqp = quad -> get_num_quadPts();

  Zero_sur_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    double surface_area;

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual0[3*A+0] += surface_area * quad -> get_qw(qua) * R[A] * n_out.x() * val;
      sur_Residual0[3*A+1] += surface_area * quad -> get_qw(qua) * R[A] * n_out.y() * val;
      sur_Residual0[3*A+2] += surface_area * quad -> get_qw(qua) * R[A] * n_out.z() * val;
    }
  }
}

void PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::Assem_Residual_BackFlowStab(
    const double * const &dot_disp,
    const double * const &disp,
    const double * const &velo,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  double curPt_x[snLocBas], curPt_y[snLocBas], curPt_z[snLocBas];

  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp, snLocBas, curPt_x, curPt_y, curPt_z);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const int face_nqp = quad -> get_num_quadPts();

  double surface_area, factor;

  Zero_sur_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    double u = 0.0, v = 0.0, w = 0.0, mu = 0.0, mv = 0.0, mw = 0.0;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      mu += dot_disp[ii*3+0] * R[ii];
      mv += dot_disp[ii*3+1] * R[ii];
      mw += dot_disp[ii*3+2] * R[ii];

      u += velo[ii*3+0] * R[ii];
      v += velo[ii*3+1] * R[ii];
      w += velo[ii*3+2] * R[ii];
    }

    const double cu = u - mu;
    const double cv = v - mv;
    const double cw = w - mw;

    const double temp = cu * n_out.x() + cv * n_out.y() + cw * n_out.z();

    if(temp < 0.0) factor = temp * rho0 * beta;
    else factor = 0.0;

    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual0[3*A+0] -= surface_area * quad -> get_qw(qua) * R[A] * factor * u;
      sur_Residual0[3*A+1] -= surface_area * quad -> get_qw(qua) * R[A] * factor * v;
      sur_Residual0[3*A+2] -= surface_area * quad -> get_qw(qua) * R[A] * factor * w;
    }
  }
}

void PLocAssem_2x2Block_ALE_VMS_NS_GenAlpha::Assem_Tangent_Residual_BackFlowStab(
    const double &dt,
    const double * const &dot_disp,
    const double * const &disp,
    const double * const &velo,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  double curPt_x[snLocBas], curPt_y[snLocBas], curPt_z[snLocBas];

  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp, snLocBas, curPt_x, curPt_y, curPt_z);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const int face_nqp = quad -> get_num_quadPts();

  double surface_area, factor;

  const double dd_dv = alpha_f * gamma * dt;

  Zero_sur_Tangent_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    double u = 0.0, v = 0.0, w = 0.0, mu = 0.0, mv = 0.0, mw = 0.0;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      mu += dot_disp[ii*3+0] * R[ii];
      mv += dot_disp[ii*3+1] * R[ii];
      mw += dot_disp[ii*3+2] * R[ii];

      u += velo[ii*3+0] * R[ii];
      v += velo[ii*3+1] * R[ii];
      w += velo[ii*3+2] * R[ii];
    }

    const double cu = u - mu;
    const double cv = v - mv;
    const double cw = w - mw;

    const double temp = cu * n_out.x() + cv * n_out.y() + cw * n_out.z();

    if(temp < 0.0) factor = temp * rho0 * beta;
    else factor = 0.0;

    const double gwts = surface_area * quad -> get_qw(qua);

    // snLocBas = 3 for linear tet/tri element
    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual0[3*A+0] -= gwts * R[A] * factor * u;
      sur_Residual0[3*A+1] -= gwts * R[A] * factor * v;
      sur_Residual0[3*A+2] -= gwts * R[A] * factor * w;

      for(int B=0; B<snLocBas; ++B)
      {
        sur_Tangent00[ 3*snLocBas*(3*A+0) + 3*B+0 ] -= gwts * dd_dv * R[A] * factor * R[B];
        sur_Tangent00[ 3*snLocBas*(3*A+1) + 3*B+1 ] -= gwts * dd_dv * R[A] * factor * R[B];
        sur_Tangent00[ 3*snLocBas*(3*A+2) + 3*B+2 ] -= gwts * dd_dv * R[A] * factor * R[B];
      }
    }
  }
}

// EOF
