#include "PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha.hpp"

PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha::PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha(
    const TimeMethod_GenAlpha * const &tm_gAlpha, const int &in_nqp, 
    const double &in_rho, const double &in_vis_mu, const double &in_beta )
: rho0( in_rho ), vis_mu( in_vis_mu ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), beta(in_beta), CI(36.0), CT(4.0),
  nLocBas_v(4), nLocBas_p(4), snLocBas_v(3), snLocBas_p(3),
  vec_size_0( nLocBas_v * 3 ), vec_size_1( nLocBas_p ), 
  sur_size_0( snLocBas_v * 3 ), sur_size_1( snLocBas_p ), nqp( in_nqp )
{
  Tangent00 = new PetscScalar[vec_size_0 * vec_size_0];
  Tangent01 = new PetscScalar[vec_size_0 * vec_size_1];
  Tangent10 = new PetscScalar[vec_size_1 * vec_size_0];
  Tangent11 = new PetscScalar[vec_size_1 * vec_size_1];

  Residual0 = new PetscScalar[vec_size_0];
  Residual1 = new PetscScalar[vec_size_1];
  
  sur_Tangent00 = new PetscScalar[sur_size_0 * sur_size_0];
  sur_Tangent01 = new PetscScalar[sur_size_0 * sur_size_1];
  sur_Tangent10 = new PetscScalar[sur_size_1 * sur_size_0];
  sur_Tangent11 = new PetscScalar[sur_size_1 * sur_size_1];

  sur_Residual0 = new PetscScalar[sur_size_0];
  sur_Residual1 = new PetscScalar[sur_size_1];

  Zero_Tangent_Residual();
  Zero_sur_Tangent_Residual();

  // print info of this assembly routine
  print_info();
}


PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha::~PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha()
{}


void PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha::print_info() const
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

void PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha::get_metric(
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

void PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha::get_tau(
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


double PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha::get_DC(
    const double * const &dxi_dx,
    const double &u, const double &v, const double &w ) const
{
  //double G11, G12, G13, G22, G23, G33;
  //get_metric( dxi_dx, G11, G12, G13, G22, G23, G33 );

  //dc_tau = G11 * u * u + 2.0 * G12 * u * v + 2.0 * G13 * u * w + G22 * v * v
  //  + 2.0 * G23 * v * w + G33 * w * w;

  //if(dc_tau > 1.0e-15) dc_tau = rho0 * std::pow(dc_tau, -0.5);
  //else dc_tau = 0.0;
  // return dc_tau;

  return 0.0;
}

void PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha::Assem_Residual(
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
  double R[4], dR_dx[4], dR_dy[4], dR_dz[4], dxi_dx[9];
  double curPt_x[4], curPt_y[4], curPt_z[4];

  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, disp, nLocBas, curPt_x, curPt_y, curPt_z);

  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

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
    double mu = 0.0, mv = 0.0, mw = 0.0; // mesh velocity, i.e. hat-v

    element->get_R_gradR( qua, R, dR_dx, dR_dy, dR_dz );
    element->get_invJacobian( qua, dxi_dx );

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

      coor_x += curPt_x[ii] * R[ii];
      coor_y += curPt_y[ii] * R[ii];
      coor_z += curPt_z[ii] * R[ii];
    }

    // v - hat(v)
    const double cu = u - mu;
    const double cv = v - mv;
    const double cw = w - mw;

    double tau_m, tau_c;
    get_tau(tau_m, tau_c, dt, dxi_dx, cu, cv, cw);

    const double tau_m_2 = tau_m * tau_m;

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const Vector_3 f_body = get_f(coor_x, coor_y, coor_z, curr);

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





























// EOF
