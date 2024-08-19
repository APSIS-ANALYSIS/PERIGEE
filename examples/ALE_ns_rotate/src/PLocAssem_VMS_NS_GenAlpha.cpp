#include "PLocAssem_VMS_NS_GenAlpha.hpp"

PLocAssem_VMS_NS_GenAlpha::PLocAssem_VMS_NS_GenAlpha(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas,
        const double &in_rho, const double &in_vis_mu,
        const double &in_beta, const int &elemtype, 
        const double &angular,
        const Vector_3 &point_xyz, const Vector_3 &angular_direc,
        const double &in_ct, const double &in_ctauc )
: rho0( in_rho ), vis_mu( in_vis_mu ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), beta(in_beta),
  CI( (elemtype == 501 || elemtype == 601) ? 36.0 : 60.0 ),
  CT( in_ct ), Ctauc( in_ctauc ),
  nqp(in_nqp), nLocBas( in_nlocbas ), snLocBas( in_snlocbas ),
  vec_size( in_nlocbas * 4 ), sur_size ( in_snlocbas * 4 ), angular_velo(angular),
  point_rotated(point_xyz), direction_rotated(angular_direc),
  coef( (elemtype == 501 || elemtype == 502) ? 0.6299605249474365 : 1.0 ),
  mm( (elemtype == 501 || elemtype == 502) ? std::array<double, 9>{2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0} :
                                             std::array<double, 9>{1.0, 0.0, 0.0, 0.0, 1.0 ,0.0, 0.0, 0.0 ,1.0} )
{
  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  sur_Tangent = new PetscScalar[sur_size * sur_size];
  sur_Residual = new PetscScalar[sur_size];

  // const int disp_size = nLocBas * 3;
  // disp_mesh = new PetscScalar[disp_size];

  // for(int ii=0; ii<disp_size; ++ii) disp_mesh[ii] = 0.0;

  Zero_Tangent_Residual();

  Zero_sur_Tangent_Residual();

  flist = new locassem_vms_ns_funs[2];
  flist[0] = &PLocAssem_VMS_NS_GenAlpha::get_Poiseuille_traction;
  flist[1] = &PLocAssem_VMS_NS_GenAlpha::get_Poiseuille_traction;

  print_info();
}

PLocAssem_VMS_NS_GenAlpha::~PLocAssem_VMS_NS_GenAlpha()
{
  delete [] Tangent; Tangent = nullptr; 
  delete [] Residual; Residual = nullptr;
  delete [] sur_Tangent; sur_Tangent = nullptr;
  delete [] sur_Residual; sur_Residual = nullptr;

  delete [] flist;
}

void PLocAssem_VMS_NS_GenAlpha::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("  Three-dimensional Incompressible Navier-Stokes equations: \n");
  if(nLocBas == 4)
    SYS_T::commPrint("  FEM: 4-node Tetrahedral element \n");
  else if(nLocBas == 10)
    SYS_T::commPrint("  FEM: 10-node Tetrahedral element \n");
  else if(nLocBas == 8)
    SYS_T::commPrint("  FEM: 8-node Hexahedral element \n");
  else if(nLocBas == 27)
    SYS_T::commPrint("  FEM: 27-node Hexahedral element \n");
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

SymmTensor2_3D PLocAssem_VMS_NS_GenAlpha::get_metric(
    const std::array<double, 9> &f ) const
{
  const double fk0 = mm[0] * f[0] + (mm[1] * f[3] + mm[2] * f[6]);
  const double fk1 = mm[4] * f[3] + (mm[3] * f[0] + mm[5] * f[6]);
  const double fk2 = mm[8] * f[6] + (mm[6] * f[0] + mm[7] * f[3]);
  const double fk3 = mm[0] * f[1] + (mm[1] * f[4] + mm[2] * f[7]);
  const double fk4 = mm[4] * f[4] + (mm[3] * f[1] + mm[5] * f[7]);
  const double fk5 = mm[8] * f[7] + (mm[6] * f[1] + mm[7] * f[4]);
  const double fk6 = mm[0] * f[2] + (mm[1] * f[5] + mm[2] * f[8]);
  const double fk7 = mm[4] * f[5] + (mm[3] * f[2] + mm[5] * f[8]);
  const double fk8 = mm[8] * f[8] + (mm[6] * f[2] + mm[7] * f[5]);

  return SymmTensor2_3D( coef * ( fk0 * f[0] + fk1 * f[3] + fk2 * f[6] ),
  coef * ( fk3 * f[1] + fk4 * f[4] + fk5 * f[7] ),
  coef * ( fk6 * f[2] + fk7 * f[5] + fk8 * f[8] ),
  coef * ( fk3 * f[2] + fk4 * f[5] + fk5 * f[8] ),
  coef * ( fk0 * f[2] + fk1 * f[5] + fk2 * f[8] ),
  coef * ( fk0 * f[1] + fk1 * f[4] + fk2 * f[7] ) );
}

std::array<double, 2> PLocAssem_VMS_NS_GenAlpha::get_tau(
    const double &dt, const std::array<double, 9> &dxi_dx,
    const double &u, const double &v, const double &w ) const
{
  const SymmTensor2_3D G = get_metric( dxi_dx );

  const Vector_3 velo_vec( u, v, w );

  const double temp_nu = vis_mu / rho0;

  const double denom_m = std::sqrt( CT / (dt*dt) + G.VecMatVec( velo_vec, velo_vec) + CI * temp_nu * temp_nu * G.MatContraction( G ) );

  // return tau_m followed by tau_c
  return {{1.0 / ( rho0 * denom_m ), Ctauc * rho0 * denom_m / G.tr()}};
}

double PLocAssem_VMS_NS_GenAlpha::get_DC(
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

void PLocAssem_VMS_NS_GenAlpha::Assem_Residual(
    const double &time, const double &dt,
    const double * const &dot_sol,
    const double * const &sol,
    const double * const &mvelo,
    const double * const &mdisp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  std::vector<double> curPt_x(nLocBas, 0.0), curPt_y(nLocBas, 0.0), curPt_z(nLocBas, 0.0);

  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, mdisp, nLocBas, &curPt_x[0], &curPt_y[0], &curPt_z[0]);

  element->buildBasis( quad, &curPt_x[0], &curPt_y[0], &curPt_z[0] );

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
    double mu = 0.0, mv = 0.0, mw = 0.0; // mesh velocity, i.e. hat-v

    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_3D_R_gradR_LaplacianR( qua, &R[0], &dR_dx[0], 
        &dR_dy[0], &dR_dz[0], &d2R_dxx[0], &d2R_dyy[0], &d2R_dzz[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int ii3 = 3 * ii;
      const int ii4 = 4 * ii;

      mu += mvelo[ii3+0] * R[ii];
      mv += mvelo[ii3+1] * R[ii];
      mw += mvelo[ii3+2] * R[ii];

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

      coor.x() += curPt_x[ii] * R[ii];
      coor.y() += curPt_y[ii] * R[ii];
      coor.z() += curPt_z[ii] * R[ii];
    }

    // v - hat(v)
    const double cu = u - mu;
    const double cv = v - mv;
    const double cw = w - mw;

    // Get the tau_m and tau_c
    const auto dxi_dx = element->get_invJacobian(qua);

    const std::array<double, 2> tau = get_tau( dt, dxi_dx, cu, cv, cw );
    const double tau_m = tau[0];
    const double tau_c = tau[1];

    const double tau_m_2 = tau_m * tau_m;

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    // Get the body force
    const Vector_3 f_body = get_f( coor, curr );

    const double u_lap = u_xx + u_yy + u_zz;
    const double v_lap = v_xx + v_yy + v_zz;
    const double w_lap = w_xx + w_yy + w_zz;

    const double rx = rho0 * ( u_t + u_x * cu + u_y * cv + u_z * cw - f_body.x() ) + p_x - vis_mu * u_lap;
    const double ry = rho0 * ( v_t + v_x * cu + v_y * cv + v_z * cw - f_body.y() ) + p_y - vis_mu * v_lap;
    const double rz = rho0 * ( w_t + w_x * cu + w_y * cv + w_z * cw - f_body.z() ) + p_z - vis_mu * w_lap;

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
      const double velo_dot_gradR = NA_x * cu + NA_y * cv + NA_z * cw;
      const double r_dot_gradR = NA_x * rx + NA_y * ry + NA_z * rz;
      const double velo_prime_dot_gradR = NA_x * u_prime + NA_y * v_prime + NA_z * w_prime;

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
          - NA * rho0 * f_body.x() );

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
          - NA * rho0 * f_body.y() );

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
          - NA * rho0 * f_body.z() );
    }
  }
}

void PLocAssem_VMS_NS_GenAlpha::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &dot_sol,
    const double * const &sol,
    const double * const &mvelo,
    const double * const &mdisp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  std::vector<double> curPt_x(nLocBas, 0.0), curPt_y(nLocBas, 0.0), curPt_z(nLocBas, 0.0);

  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, mdisp, nLocBas, &curPt_x[0], &curPt_y[0], &curPt_z[0]);

  element->buildBasis( quad, &curPt_x[0], &curPt_y[0], &curPt_z[0] );

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
    double mu = 0.0, mv = 0.0, mw = 0.0; // mesh velocity, i.e. hat-v

    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_3D_R_gradR_LaplacianR( qua, &R[0], &dR_dx[0], 
        &dR_dy[0], &dR_dz[0], &d2R_dxx[0], &d2R_dyy[0], &d2R_dzz[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int ii3 = 3 * ii;      
      const int ii4 = 4 * ii;

      mu += mvelo[ii3+0] * R[ii];
      mv += mvelo[ii3+1] * R[ii];
      mw += mvelo[ii3+2] * R[ii];

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

      coor.x() += curPt_x[ii] * R[ii];
      coor.y() += curPt_y[ii] * R[ii];
      coor.z() += curPt_z[ii] * R[ii];
    }

    // v - hat(v)
    const double cu = u - mu;
    const double cv = v - mv;
    const double cw = w - mw;

    const auto dxi_dx = element->get_invJacobian(qua);

    const std::array<double, 2> tau = get_tau( dt, dxi_dx, cu, cv, cw );
    const double tau_m = tau[0];
    const double tau_c = tau[1];

    const double tau_m_2 = tau_m * tau_m;

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua); 

    const Vector_3 f_body = get_f( coor, curr );

    const double u_lap = u_xx + u_yy + u_zz;
    const double v_lap = v_xx + v_yy + v_zz;
    const double w_lap = w_xx + w_yy + w_zz;

    const double rx = rho0 * ( u_t + u_x * cu + u_y * cv + u_z * cw - f_body.x() ) + p_x - vis_mu * u_lap;
    const double ry = rho0 * ( v_t + v_x * cu + v_y * cv + v_z * cw - f_body.y() ) + p_y - vis_mu * v_lap ;
    const double rz = rho0 * ( w_t + w_x * cu + w_y * cv + w_z * cw - f_body.z() ) + p_z - vis_mu * w_lap;

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

      const double velo_dot_gradR = NA_x * cu + NA_y * cv + NA_z * cw;
      const double r_dot_gradR = NA_x * rx + NA_y * ry + NA_z * rz;
      const double velo_prime_dot_gradR = NA_x * u_prime + NA_y * v_prime + NA_z * w_prime;

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
          - NA * rho0 * f_body.x() );

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
          - NA * rho0 * f_body.y() );

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
          - NA * rho0 * f_body.z() );

      for(int B=0; B<nLocBas; ++B)
      {
        const double NB = R[B], NB_x = dR_dx[B], NB_y = dR_dy[B], NB_z = dR_dz[B];
        const double NB_xx = d2R_dxx[B], NB_yy = d2R_dyy[B], NB_zz = d2R_dzz[B];
        const double NB_lap = NB_xx + NB_yy + NB_zz;
        const double velo_dot_gradNB = cu * NB_x + cv * NB_y + cw * NB_z;
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

void PLocAssem_VMS_NS_GenAlpha::Assem_Mass_Residual(
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

    Vector_3 coor(0.0, 0.0, 0.0);

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

    const Vector_3 f_body = get_f( coor, curr );

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

void PLocAssem_VMS_NS_GenAlpha::Assem_Residual_EBC(
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

  Zero_sur_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    double surface_area;

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    Vector_3 coor(0.0, 0.0, 0.0);
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

void PLocAssem_VMS_NS_GenAlpha::Assem_Residual_EBC_Resistance(
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

  double surface_area;

  Zero_sur_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual[4*A+1] += surface_area * quad -> get_qw(qua) * R[A] * n_out.x() * val;
      sur_Residual[4*A+2] += surface_area * quad -> get_qw(qua) * R[A] * n_out.y() * val;
      sur_Residual[4*A+3] += surface_area * quad -> get_qw(qua) * R[A] * n_out.z() * val;
    }
  }
}

void PLocAssem_VMS_NS_GenAlpha::Assem_Residual_BackFlowStab(
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  Zero_sur_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    double surface_area;

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    Vector_3 velo(0.0, 0.0, 0.0);
    for(int ii=0; ii<snLocBas; ++ii)
    {
      const int ii4 = ii * 4;
      velo.x() += sol[ii4+1] * R[ii];
      velo.y() += sol[ii4+2] * R[ii];
      velo.z() += sol[ii4+3] * R[ii];
    }

    const double temp = Vec3::dot_product( velo, n_out );

    const double factor = temp < 0.0 ? temp * rho0 * beta : 0.0;

    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual[4*A+1] -= surface_area * quad -> get_qw(qua) * R[A] * factor * velo.x();
      sur_Residual[4*A+2] -= surface_area * quad -> get_qw(qua) * R[A] * factor * velo.y();
      sur_Residual[4*A+3] -= surface_area * quad -> get_qw(qua) * R[A] * factor * velo.z();
    }
  }
}

void PLocAssem_VMS_NS_GenAlpha::Assem_Tangent_Residual_BackFlowStab(
    const double &dt,
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

  Zero_sur_Tangent_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    double surface_area;

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    Vector_3 velo(0.0, 0.0, 0.0);
    for(int ii=0; ii<snLocBas; ++ii)
    {
      velo.x() += sol[ii*4+1] * R[ii];
      velo.y() += sol[ii*4+2] * R[ii];
      velo.z() += sol[ii*4+3] * R[ii];
    }

    const double temp = Vec3::dot_product( velo, n_out );

    const double factor = temp < 0.0 ? temp * rho0 * beta : 0.0;

    const double gwts = surface_area * quad -> get_qw(qua);

    // snLocBas = 3 for linear tri element
    //            6 for quadratic tri element
    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual[4*A+1] -= gwts * R[A] * factor * velo.x();
      sur_Residual[4*A+2] -= gwts * R[A] * factor * velo.y();
      sur_Residual[4*A+3] -= gwts * R[A] * factor * velo.z();

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

double PLocAssem_VMS_NS_GenAlpha::get_flowrate( const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  double flrate = 0.0;

  for(int qua =0; qua< face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    double surface_area;
    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    Vector_3 velo(0.0, 0.0, 0.0);
    for(int ii=0; ii<snLocBas; ++ii)
    {
      velo.x() += sol[ii*4+1] * R[ii];
      velo.y() += sol[ii*4+2] * R[ii];
      velo.z() += sol[ii*4+3] * R[ii];
    }

    flrate += surface_area * quad->get_qw(qua) * Vec3::dot_product( velo, n_out ); 
  }

  return flrate;
}

void PLocAssem_VMS_NS_GenAlpha::get_pressure_area( 
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

// EOF
