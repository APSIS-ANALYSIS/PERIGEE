#include "PLocAssem_VMS_NS_SemiBDF1.hpp"

PLocAssem_VMS_NS_SemiBDF1::PLocAssem_VMS_NS_SemiBDF1(
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas,
        const double &in_rho, const double &in_vis_mu,
        const double &in_beta, const int &elemtype, 
        const double &in_ct, const double &in_ctauc )
: rho0( in_rho ), vis_mu( in_vis_mu ),
  CI( (elemtype == 501 || elemtype == 601) ? 36.0 : 60.0 ),
  CT( in_ct ), Ctauc( in_ctauc ),
  nqp(in_nqp), nLocBas( in_nlocbas ), snLocBas( in_snlocbas ),
  vec_size( in_nlocbas * 4 ), sur_size ( in_snlocbas * 4 ),
  coef( (elemtype == 501 || elemtype == 502) ? 0.6299605249474365 : 1.0 ),
  mm( (elemtype == 501 || elemtype == 502) ? std::array<double, 9>{2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0} :
                                             std::array<double, 9>{1.0, 0.0, 0.0, 0.0, 1.0 ,0.0, 0.0, 0.0 ,1.0} )
{
  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  sur_Tangent = new PetscScalar[sur_size * sur_size];
  sur_Residual = new PetscScalar[sur_size];

  Zero_Tangent_Residual();

  Zero_sur_Tangent_Residual();

  print_info();
}

PLocAssem_VMS_NS_GenAlpha::~PLocAssem_VMS_NS_SemiBDF1()
{
  delete [] Tangent; Tangent = nullptr; 
  delete [] Residual; Residual = nullptr;
  delete [] sur_Tangent; sur_Tangent = nullptr;
  delete [] sur_Residual; sur_Residual = nullptr;
}

void PLocAssem_VMS_NS_SemiBDF1::print_info() const
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
  SYS_T::commPrint("  Temporal: Semi-BDF-1st Method \n");
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
  SYS_T::commPrint("----------------------------------------------------------- \n");
}

SymmTensor2_3D PLocAssem_VMS_NS_SemiBDF1::get_metric(
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

std::array<double, 2> PLocAssem_VMS_NS_SemiBDF1::get_tau(
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

double PLocAssem_VMS_NS_SemiBDF1::get_DC(
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

void PLocAssem_VMS_NS_SemiBDF1::Assem_Residual(
    const double &time, const double &dt,
    const double * const &sol_n,
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double two_mu = 2.0 * vis_mu;

  const double curr = time + dt;

  Zero_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);
  std::vector<double> d2R_dxx(nLocBas, 0.0), d2R_dyy(nLocBas, 0.0), d2R_dzz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    // variables at t_n
    double u_n = 0.0, u_x_n = 0.0, u_y_n = 0.0, u_z_n = 0.0;
    double v_n = 0.0, v_x_n = 0.0, v_y_n = 0.0, v_z_n = 0.0;
    double w_n = 0.0, w_x_n = 0.0, w_y_n = 0.0, w_z_n = 0.0;
    double p_n = 0.0, p_x_n = 0.0, p_y_n = 0.0, p_z_n = 0.0;
    double u_xx_n = 0.0, u_yy_n = 0.0, u_zz_n = 0.0;
    double v_xx_n = 0.0, v_yy_n = 0.0, v_zz_n = 0.0;
    double w_xx_n = 0.0, w_yy_n = 0.0, w_zz_n = 0.0;

    // variables at t_n+1
    double u = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v = 0.0, v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w = 0.0, w_x = 0.0, w_y = 0.0, w_z = 0.0;
    double p = 0.0, p_x = 0.0, p_y = 0.0, p_z = 0.0;
    double u_xx = 0.0, u_yy = 0.0, u_zz = 0.0;
    double v_xx = 0.0, v_yy = 0.0, v_zz = 0.0;
    double w_xx = 0.0, w_yy = 0.0, w_zz = 0.0;


    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_3D_R_gradR_LaplacianR( qua, &R[0], &dR_dx[0], 
        &dR_dy[0], &dR_dz[0], &d2R_dxx[0], &d2R_dyy[0], &d2R_dzz[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int ii4 = 4 * ii;

      // variables at t_n
      u_n += sol_n[ii4+1] * R[ii];
      v_n += sol_n[ii4+2] * R[ii];
      w_n += sol_n[ii4+3] * R[ii];
      p_n += sol_n[ii4+0] * R[ii];

      u_x_n += sol_n[ii4+1] * dR_dx[ii];
      v_x_n += sol_n[ii4+2] * dR_dx[ii];
      w_x_n += sol_n[ii4+3] * dR_dx[ii];
      p_x_n += sol_n[ii4+0] * dR_dx[ii];

      u_y_n += sol_n[ii4+1] * dR_dy[ii];
      v_y_n += sol_n[ii4+2] * dR_dy[ii];
      w_y_n += sol_n[ii4+3] * dR_dy[ii];
      p_y_n += sol_n[ii4+0] * dR_dy[ii];

      u_z_n += sol_n[ii4+1] * dR_dz[ii];
      v_z_n += sol_n[ii4+2] * dR_dz[ii];
      w_z_n += sol_n[ii4+3] * dR_dz[ii];
      p_z_n += sol_n[ii4+0] * dR_dz[ii];

      u_xx_n += sol_n[ii4+1] * d2R_dxx[ii];
      u_yy_n += sol_n[ii4+1] * d2R_dyy[ii];
      u_zz_n += sol_n[ii4+1] * d2R_dzz[ii];

      v_xx_n += sol_n[ii4+2] * d2R_dxx[ii];
      v_yy_n += sol_n[ii4+2] * d2R_dyy[ii];
      v_zz_n += sol_n[ii4+2] * d2R_dzz[ii];

      w_xx_n += sol_n[ii4+3] * d2R_dxx[ii];
      w_yy_n += sol_n[ii4+3] * d2R_dyy[ii];
      w_zz_n += sol_n[ii4+3] * d2R_dzz[ii];

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];

      // variables at t_n+1 
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

    }
    // 1st-order BDF scheme  
    const double alpha_sig = 1;

    const double u_n_BDF = u_n;
    const double v_n_BDF = v_n;
    const double w_n_BDF = z_n;

    // Newton-Gregory Polynomial
    const double u_NG = u_n;
    const double v_NG = v_n;
    const double w_NG = w_n;
    const double p_NG = p_n;

    const double u_x_NG = u_x_n;
    const double v_x_NG = v_x_n;
    const double w_x_NG = w_x_n;
    const double p_x_NG = p_x_n;

    const double u_y_NG = u_y_n;
    const double v_y_NG = v_y_n;
    const double w_y_NG = w_y_n;
    const double p_y_NG = p_y_n;

    const double u_z_NG = u_z_n;
    const double v_z_NG = v_z_n;
    const double w_z_NG = w_z_n;
    const double p_z_NG = p_z_n;

    const double u_xx_NG = u_xx_n;
    const double v_xx_NG = v_xx_n;
    const double w_xx_NG = w_xx_n;

    const double u_yy_NG = u_yy_n;
    const double v_yy_NG = v_yy_n;
    const double w_yy_NG = w_yy_n;

    const double u_zz_NG = u_zz_n;
    const double v_zz_NG = v_zz_n;
    const double w_zz_NG = w_zz_n;

    // Compute the laplacian of vector
    const double u_lap_NG = u_xx_n + u_yy_n + u_zz_n;
    const double v_lap_NG = v_xx_n + v_yy_n + v_zz_n;
    const double w_lap_NG = w_xx_n + w_yy_n + w_zz_n;

    const double u_lap = u_xx + u_yy + u_zz;
    const double v_lap = v_xx + v_yy + v_zz;
    const double w_lap = w_xx + w_yy + w_zz;

    // Get the tau_m and tau_c
    const auto dxi_dx = element->get_invJacobian(qua);

    const std::array<double, 2> tau = get_tau( dt, dxi_dx, u_NG, v_NG, w_NG );
    const double tau_m = tau[0];
    const double tau_c = tau[1];

    const double tau_m_2 = tau_m * tau_m;

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    // Get the body force
    const Vector_3 f_body_p0 = get_f( coor, time );
    const Vector_3 f_body_p1 = get_f( coor, curr );


    const double w_lap_p1 = w_xx_p1 + w_yy_p1 + w_zz_p1;

    // Compute the residue 
    const double rx_hat = rho0 * ( (alpha_sig * u_NG - u_n_BDF)/dt + u_NG * u_x_NG + v_NG * u_y_NG + w_NG * u_z_NG - f_body.x() ) + p_x_NG - vis_mu * u_lap_NG; 
    const double ry_hat = rho0 * ( (alpha_sig * v_NG - v_n_BDF)/dt + u_NG * v_x_NG + v_NG * v_y_NG + w_NG * v_z_NG - f_body.y() ) + p_y_NG - vis_mu * v_lap_NG; 
    const double rz_hat = rho0 * ( (alpha_sig * w_NG - w_n_BDF)/dt + u_NG * w_x_NG + v_NG * w_y_NG + w_NG * w_z_NG - f_body.z() ) + p_z_NG - vis_mu * w_lap_NG;  

    const double rx_tilda = rho0 * ( u_NG * u_x + v_NG * u_y + w_NG * u_z - f_body.x() ) + p_x - vis_mu * u_lap;
    const double ry_tilda = rho0 * ( u_NG * v_x + v_NG * v_y + w_NG * v_z - f_body.y() ) + p_y - vis_mu * v_lap;
    const double rz_tilda = rho0 * ( u_NG * w_x + v_NG * w_y + w_NG * w_z - f_body.z() ) + p_z - vis_mu * w_lap;

    const double rx = rho0 * ( (alpha_sig * u - u_n_BDF)/dt + u_NG * u_x + v_NG * u_y + w_NG * u_z - f_body.x() ) + p_x - vis_mu * u_lap;
    const double ry = rho0 * ( (alpha_sig * v - v_n_BDF)/dt + u_NG * v_x + v_NG * v_y + w_NG * v_z - f_body.y() ) + p_y - vis_mu * v_lap;
    const double rz = rho0 * ( (alpha_sig * w - w_n_BDF)/dt + u_NG * w_x + v_NG * w_y + w_NG * w_z - f_body.z() ) + p_z - vis_mu * w_lap;

    const double div_vel = u_x + v_y + w_z;

    const double u_prime = -1.0 * tau_m * rx;
    const double v_prime = -1.0 * tau_m * ry;
    const double w_prime = -1.0 * tau_m * rz;

    const double r_dot_gradu_NG = u_x_NG * rx + u_y_NG * ry + u_z_NG * rz;
    const double r_dot_gradv_NG = v_x_NG * rx + v_y_NG * ry + v_z_NG * rz;
    const double r_dot_gradw_NG = w_x_NG * rx + w_y_NG * ry + w_z_NG * rz;
    
    // Get the Discontinuity Capturing tau
    const double tau_dc = get_DC( dxi_dx, u_prime, v_prime, w_prime );

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];
      
      const double velo_NG_dot_gradR    = NA_x * u_NG    + NA_y * v_NG    + NA_z * w_NG;
      const double r_hat_dot_gradR      = NA_x * rx_hat  + NA_y * ry_hat  + NA_z * rz_hat;
      const double r_dot_gradR          = NA_x * rx      + NA_y * ry.     + NA_z * rz;
      const double velo_prime_dot_gradR = NA_x * u_prime + NA_y * v_prime + NA_z * w_prime;

      Residual[4*A] += gwts * ( NA * div_vel + tau_m * r_dot_gradR );

      Residual[4*A+1] += gwts * ( NA * rho0 * (alpha_sig * u - u_n_BDF)/dt 
          + NA * rho0 * (u_NG * u_x + v_NG * u_y + w_NG * u_z)
          - NA_x * p
          + NA_x * two_mu * u_x
          + NA_y * vis_mu * (u_y + v_x)
          + NA_z * vis_mu * (u_z + w_x)
          - NA * rho0 * f_body.x()
          + velo_NG_dot_gradR * rho0 * tau_m * rx
          - NA * rho0 * tau_m * r_dot_gradu_NG
          - r_hat_dot_gradR * rho0 * tau_m_2 * rx_tilda
          - r_hat_dot_gradR * rho0 * tau_m_2 * rho0 * (alpha_sig * u)/dt
          + r_dot_gradR * rho0 * tau_m_2 * rho0 * u_n_BDF/dt
          + NA_x * tau_c * div_vel );

      Residual[4*A+2] += gwts * ( NA * rho0 * (alpha_sig * v - v_n_BDF)/dt 
          + NA * rho0 * (u_NG * v_x + v_NG * v_y + w_NG * v_z)
          - NA_y * p
          + NA_x * vis_mu * (v_x + u_y)
          + NA_y * two_mu * v_y 
          + NA_z * vis_mu * (v_z + w_y)
          - NA * rho0 * f_body.y()
          + velo_NG_dot_gradR * rho0 * tau_m * ry
          - NA * rho0 * tau_m * r_dot_gradv_NG
          - r_hat_dot_gradR * rho0 * tau_m_2 * ry_tilda
          - r_hat_dot_gradR * rho0 * tau_m_2 * rho0 * (alpha_sig * v)/dt
          + r_dot_gradR * rho0 * tau_m_2 * rho0 * v_n_BDF/dt
          + NA_y * tau_c * div_vel );

      Residual[4*A+3] += gwts * ( NA * rho0 * (alpha_sig * v - v_n_BDF)/dt 
          + NA * rho0 * (u_NG * v_x + v_NG * v_y + w_NG * v_z)
          - NA_y * p
          + NA_x * vis_mu * (v_x + u_y)
          + NA_y * two_mu * v_y 
          + NA_z * vis_mu * (v_z + w_y)
          - NA * rho0 * f_body.y()
          + velo_NG_dot_gradR * rho0 * tau_m * ry
          - NA * rho0 * tau_m * r_dot_gradv_NG
          - r_hat_dot_gradR * rho0 * tau_m_2 * ry_tilda
          - r_hat_dot_gradR * rho0 * tau_m_2 * rho0 * (alpha_sig * v)/dt
          + r_dot_gradR * rho0 * tau_m_2 * rho0 * v_n_BDF/dt
          + NA_y * tau_c * div_vel );
 
    } // A-loop
  } //qua-loop
}

void PLocAssem_VMS_NS_SemiBDF1::Assem_Tanget_Residual(
    const double &time, const double &dt,
    const double * const &sol_n,
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double two_mu = 2.0 * vis_mu;

  const double curr = time + dt;

  Zero_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);
  std::vector<double> d2R_dxx(nLocBas, 0.0), d2R_dyy(nLocBas, 0.0), d2R_dzz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    // variables at t_n
    double u_n = 0.0, u_x_n = 0.0, u_y_n = 0.0, u_z_n = 0.0;
    double v_n = 0.0, v_x_n = 0.0, v_y_n = 0.0, v_z_n = 0.0;
    double w_n = 0.0, w_x_n = 0.0, w_y_n = 0.0, w_z_n = 0.0;
    double p_n = 0.0, p_x_n = 0.0, p_y_n = 0.0, p_z_n = 0.0;
    double u_xx_n = 0.0, u_yy_n = 0.0, u_zz_n = 0.0;
    double v_xx_n = 0.0, v_yy_n = 0.0, v_zz_n = 0.0;
    double w_xx_n = 0.0, w_yy_n = 0.0, w_zz_n = 0.0;

    // variables at t_n+1
    double u = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v = 0.0, v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w = 0.0, w_x = 0.0, w_y = 0.0, w_z = 0.0;
    double p = 0.0, p_x = 0.0, p_y = 0.0, p_z = 0.0;
    double u_xx = 0.0, u_yy = 0.0, u_zz = 0.0;
    double v_xx = 0.0, v_yy = 0.0, v_zz = 0.0;
    double w_xx = 0.0, w_yy = 0.0, w_zz = 0.0;


    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_3D_R_gradR_LaplacianR( qua, &R[0], &dR_dx[0], 
        &dR_dy[0], &dR_dz[0], &d2R_dxx[0], &d2R_dyy[0], &d2R_dzz[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int ii4 = 4 * ii;

      // variables at t_n
      u_n += sol_n[ii4+1] * R[ii];
      v_n += sol_n[ii4+2] * R[ii];
      w_n += sol_n[ii4+3] * R[ii];
      p_n += sol_n[ii4+0] * R[ii];

      u_x_n += sol_n[ii4+1] * dR_dx[ii];
      v_x_n += sol_n[ii4+2] * dR_dx[ii];
      w_x_n += sol_n[ii4+3] * dR_dx[ii];
      p_x_n += sol_n[ii4+0] * dR_dx[ii];

      u_y_n += sol_n[ii4+1] * dR_dy[ii];
      v_y_n += sol_n[ii4+2] * dR_dy[ii];
      w_y_n += sol_n[ii4+3] * dR_dy[ii];
      p_y_n += sol_n[ii4+0] * dR_dy[ii];

      u_z_n += sol_n[ii4+1] * dR_dz[ii];
      v_z_n += sol_n[ii4+2] * dR_dz[ii];
      w_z_n += sol_n[ii4+3] * dR_dz[ii];
      p_z_n += sol_n[ii4+0] * dR_dz[ii];

      u_xx_n += sol_n[ii4+1] * d2R_dxx[ii];
      u_yy_n += sol_n[ii4+1] * d2R_dyy[ii];
      u_zz_n += sol_n[ii4+1] * d2R_dzz[ii];

      v_xx_n += sol_n[ii4+2] * d2R_dxx[ii];
      v_yy_n += sol_n[ii4+2] * d2R_dyy[ii];
      v_zz_n += sol_n[ii4+2] * d2R_dzz[ii];

      w_xx_n += sol_n[ii4+3] * d2R_dxx[ii];
      w_yy_n += sol_n[ii4+3] * d2R_dyy[ii];
      w_zz_n += sol_n[ii4+3] * d2R_dzz[ii];

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];

      // variables at t_n+1 
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

    }
    // 1st-order BDF scheme  
    const double alpha_sig = 1;

    const double u_n_BDF = u_n;
    const double v_n_BDF = v_n;
    const double w_n_BDF = z_n;

    // Newton-Gregory Polynominal
    const double u_NG = u_n;
    const double v_NG = v_n;
    const double w_NG = w_n;
    const double p_NG = p_n;

    const double u_x_NG = u_x_n;
    const double v_x_NG = v_x_n;
    const double w_x_NG = w_x_n;
    const double p_x_NG = p_x_n;

    const double u_y_NG = u_y_n;
    const double v_y_NG = v_y_n;
    const double w_y_NG = w_y_n;
    const double p_y_NG = p_y_n;

    const double u_z_NG = u_z_n;
    const double v_z_NG = v_z_n;
    const double w_z_NG = w_z_n;
    const double p_z_NG = p_z_n;

    const double u_xx_NG = u_xx_n;
    const double v_xx_NG = v_xx_n;
    const double w_xx_NG = w_xx_n;

    const double u_yy_NG = u_yy_n;
    const double v_yy_NG = v_yy_n;
    const double w_yy_NG = w_yy_n;

    const double u_zz_NG = u_zz_n;
    const double v_zz_NG = v_zz_n;
    const double w_zz_NG = w_zz_n;

    // Compute the laplacian of vector
    const double u_lap_NG = u_xx_n + u_yy_n + u_zz_n;
    const double v_lap_NG = v_xx_n + v_yy_n + v_zz_n;
    const double w_lap_NG = w_xx_n + w_yy_n + w_zz_n;

    const double u_lap = u_xx + u_yy + u_zz;
    const double v_lap = v_xx + v_yy + v_zz;
    const double w_lap = w_xx + w_yy + w_zz;

    // Get the tau_m and tau_c
    const auto dxi_dx = element->get_invJacobian(qua);

    const std::array<double, 2> tau = get_tau( dt, dxi_dx, u_NG, v_NG, w_NG );
    const double tau_m = tau[0];
    const double tau_c = tau[1];

    const double tau_m_2 = tau_m * tau_m;

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    // Get the body force
    const Vector_3 f_body_p0 = get_f( coor, time );
    const Vector_3 f_body_p1 = get_f( coor, curr );


    const double w_lap_p1 = w_xx_p1 + w_yy_p1 + w_zz_p1;

    // Compute the residue 
    const double rx_hat = rho0 * ( (alpha_sig * u_NG - u_n_BDF)/dt + u_NG * u_x_NG + v_NG * u_y_NG + w_NG * u_z_NG - f_body.x() ) + p_x_NG - vis_mu * u_lap_NG; 
    const double ry_hat = rho0 * ( (alpha_sig * v_NG - v_n_BDF)/dt + u_NG * v_x_NG + v_NG * v_y_NG + w_NG * v_z_NG - f_body.y() ) + p_y_NG - vis_mu * v_lap_NG; 
    const double rz_hat = rho0 * ( (alpha_sig * w_NG - w_n_BDF)/dt + u_NG * w_x_NG + v_NG * w_y_NG + w_NG * w_z_NG - f_body.z() ) + p_z_NG - vis_mu * w_lap_NG;  

    const double rx_tilda = rho0 * ( u_NG * u_x + v_NG * u_y + w_NG * u_z - f_body.x() ) + p_x - vis_mu * u_lap;
    const double ry_tilda = rho0 * ( u_NG * v_x + v_NG * v_y + w_NG * v_z - f_body.y() ) + p_y - vis_mu * v_lap;
    const double rz_tilda = rho0 * ( u_NG * w_x + v_NG * w_y + w_NG * w_z - f_body.z() ) + p_z - vis_mu * w_lap;

    const double rx = rho0 * ( (alpha_sig * u - u_n_BDF)/dt + u_NG * u_x + v_NG * u_y + w_NG * u_z - f_body.x() ) + p_x - vis_mu * u_lap;
    const double ry = rho0 * ( (alpha_sig * v - v_n_BDF)/dt + u_NG * v_x + v_NG * v_y + w_NG * v_z - f_body.y() ) + p_y - vis_mu * v_lap;
    const double rz = rho0 * ( (alpha_sig * w - w_n_BDF)/dt + u_NG * w_x + v_NG * w_y + w_NG * w_z - f_body.z() ) + p_z - vis_mu * w_lap;

    const double div_vel = u_x + v_y + w_z;

    const double u_prime = -1.0 * tau_m * rx;
    const double v_prime = -1.0 * tau_m * ry;
    const double w_prime = -1.0 * tau_m * rz;

    const double r_dot_gradu_NG = u_x_NG * rx + u_y_NG * ry + u_z_NG * rz;
    const double r_dot_gradv_NG = v_x_NG * rx + v_y_NG * ry + v_z_NG * rz;
    const double r_dot_gradw_NG = w_x_NG * rx + w_y_NG * ry + w_z_NG * rz;
    
    // Get the Discontinuity Capturing tau
    const double tau_dc = get_DC( dxi_dx, u_prime, v_prime, w_prime );

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];
      
      const double velo_NG_dot_gradR    = NA_x * u_NG    + NA_y * v_NG    + NA_z * w_NG;
      const double r_hat_dot_gradR      = NA_x * rx_hat  + NA_y * ry_hat  + NA_z * rz_hat;
      const double r_dot_gradR          = NA_x * rx      + NA_y * ry.     + NA_z * rz;
      const double velo_prime_dot_gradR = NA_x * u_prime + NA_y * v_prime + NA_z * w_prime;

      Residual[4*A] += gwts * ( NA * div_vel + tau_m * r_dot_gradR );

      Residual[4*A+1] += gwts * ( NA * rho0 * (alpha_sig * u - u_n_BDF)/dt 
          + NA * rho0 * (u_NG * u_x + v_NG * u_y + w_NG * u_z)
          - NA_x * p
          + NA_x * two_mu * u_x
          + NA_y * vis_mu * (u_y + v_x)
          + NA_z * vis_mu * (u_z + w_x)
          - NA * rho0 * f_body.x()
          + velo_NG_dot_gradR * rho0 * tau_m * rx
          - NA * rho0 * tau_m * r_dot_gradu_NG
          - r_hat_dot_gradR * rho0 * tau_m_2 * rx_tilda
          - r_hat_dot_gradR * rho0 * tau_m_2 * rho0 * (alpha_sig * u)/dt
          + r_dot_gradR * rho0 * tau_m_2 * rho0 * u_n_BDF/dt
          + NA_x * tau_c * div_vel );

      Residual[4*A+2] += gwts * ( NA * rho0 * (alpha_sig * v - v_n_BDF)/dt 
          + NA * rho0 * (u_NG * v_x + v_NG * v_y + w_NG * v_z)
          - NA_y * p
          + NA_x * vis_mu * (v_x + u_y)
          + NA_y * two_mu * v_y 
          + NA_z * vis_mu * (v_z + w_y)
          - NA * rho0 * f_body.y()
          + velo_NG_dot_gradR * rho0 * tau_m * ry
          - NA * rho0 * tau_m * r_dot_gradv_NG
          - r_hat_dot_gradR * rho0 * tau_m_2 * ry_tilda
          - r_hat_dot_gradR * rho0 * tau_m_2 * rho0 * (alpha_sig * v)/dt
          + r_dot_gradR * rho0 * tau_m_2 * rho0 * v_n_BDF/dt
          + NA_y * tau_c * div_vel );

      Residual[4*A+3] += gwts * ( NA * rho0 * (alpha_sig * v - v_n_BDF)/dt 
          + NA * rho0 * (u_NG * v_x + v_NG * v_y + w_NG * v_z)
          - NA_y * p
          + NA_x * vis_mu * (v_x + u_y)
          + NA_y * two_mu * v_y 
          + NA_z * vis_mu * (v_z + w_y)
          - NA * rho0 * f_body.y()
          + velo_NG_dot_gradR * rho0 * tau_m * ry
          - NA * rho0 * tau_m * r_dot_gradv_NG
          - r_hat_dot_gradR * rho0 * tau_m_2 * ry_tilda
          - r_hat_dot_gradR * rho0 * tau_m_2 * rho0 * (alpha_sig * v)/dt
          + r_dot_gradR * rho0 * tau_m_2 * rho0 * v_n_BDF/dt
          + NA_y * tau_c * div_vel );
      
      for(int B=0; B<nLocBas; ++B)
      {
        const double NB = R[B], NB_x = dR_dx[B], NB_y = dR_dy[B], NB_z = dR_dz[B];
        const double NB_xx = d2R_dxx[B], NB_yy = d2R_dyy[B], NB_zz = d2R_dzz[B];
        const double NB_lap = NB_xx + NB_yy + NB_zz;

        const double velo_dot_gradNB = u * NB_x + v * NB_y + w * NB_z;
        const double velo_NG_dot_gradNB = u_NG * NB_x + v_NG * NB_y + w_NG * NB_z;

        const double NANB  = NA*NB, NANBx = NA*NB_x, NANBy = NA*NB_y, NANBz = NA*NB_z;
        const double NAxNB = NA_x*NB, NAxNBx = NA_x*NB_x, NAxNBy = NA_x*NB_y, NAxNBz = NA_x*NB_z;
        const double NAyNB = NA_y*NB, NAyNBx = NA_y*NB_x, NAyNBy = NA_y*NB_y, NAyNBz = NA_y*NB_z;
        const double NAzNB = NA_z*NB, NAzNBx = NA_z*NB_x, NAzNBy = NA_z*NB_y, NAzNBz = NA_z*NB_z;

        const double drx_du_B = rho0 * ( alpha_sig * NB / dt + velo_NG_dot_gradNB ) - vis_mu * NB_lap; 
        const double drx_dv_B = 0;
        const double drx_dw_B = 0;

        const double dry_du_B = 0;
        const double dry_dv_B = rho0 * ( alpha_sig * NB / dt + velo_NG_dot_gradNB ) - vis_mu * NB_lap; 
        const double dry_dw_B = 0;

        const double drz_du_B = 0;
        const double drz_dv_B = 0;
        const double drz_dw_B = rho0 * ( alpha_sig * NB / dt + velo_NG_dot_gradNB ) - vis_mu * NB_lap; 

        const double drx_NG_du_B = rho0 * velo_NG_dot_gradNB - vis_mu * NB_lap; 
        const double drx_NG_dv_B = 0;
        const double drx_NG_dw_B = 0;

        const double dry_NG_du_B = 0;
        const double dry_NG_dv_B = rho0 * velo_NG_dot_gradNB - vis_mu * NB_lap; 
        const double dry_NG_dw_B = 0;

        const double drz_NG_du_B = 0;
        const double drz_NG_dv_B = 0;
        const double drz_NG_dw_B = rho0 * velo_NG_dot_gradNB - vis_mu * NB_lap; 

        // Continuity equation with respect to p, u, v, w
        Tangent[16*nLocBas*A+4*B] += gwts *  tau_m * (NAxNBx + NAyNBy + NAzNBz);
        Tangent[16*nLocBas*A+4*B+1] += gwts * tau_m * (NA_x * drx_du_B + NA_y * dry_du_B + NA_z * drz_du_B );
        Tangent[16*nLocBas*A+4*B+2] += gwts * tau_m * (NA_x * drx_dv_B + NA_y * dry_dv_B + NA_z * drz_dv_B );
        Tangent[16*nLocBas*A+4*B+3] += gwts * tau_m * (NA_x * drx_dw_B + NA_y * dry_dw_B + NA_z * drz_dw_B );

        // Momentum-x with respect to p, u, v, w
        Tangent[4*nLocBas*(4*A+1)+4*B] += gwts * ( (-1.0) * NAxNB
            + velo_NG_dot_gradR * rho0 * tau_m * NB_x
            - NA * rho0 * tau_m * (u_x_NG * NB_x + u_y_NG * NB_y + u_z_NG * NB_z)
            - r_hat_dot_gradR * rho0 * tau_m_2 *  * NB_x
            + (NAxNBx + NAyNBy + NAzNBz) * rho0 * tau_m_2 * rho0 /dt * u_n_BDF );










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
      } // B-loop
    } // A-loop
  } //qua-loop
}
// EOF
