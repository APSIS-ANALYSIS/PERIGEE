#include "PLocAssem_VMS_NS_ImexBDF1.hpp"

PLocAssem_VMS_NS_ImexBDF1::PLocAssem_VMS_NS_ImexBDF1(
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas,
        const double &in_rho, const double &in_vis_mu,
        const double &in_beta, const int &elemtype, 
        const double &in_ct, const double &in_ctauc )
: rho0( in_rho ), vis_mu( in_vis_mu ), beta(in_beta),
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

  flist = new locassem_vms_ns_funs[1];
  flist[0] = &PLocAssem_VMS_NS_ImexBDF1::get_cubic_velo_traction;

  print_info();
}

PLocAssem_VMS_NS_ImexBDF1::~PLocAssem_VMS_NS_ImexBDF1()
{
  delete [] Tangent; Tangent = nullptr; 
  delete [] Residual; Residual = nullptr;
  delete [] sur_Tangent; sur_Tangent = nullptr;
  delete [] sur_Residual; sur_Residual = nullptr;
}

void PLocAssem_VMS_NS_ImexBDF1::print_info() const
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
  SYS_T::commPrint("  Temporal: Imex-BDF-1st Method \n");
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

SymmTensor2_3D PLocAssem_VMS_NS_ImexBDF1::get_metric(
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

std::array<double, 2> PLocAssem_VMS_NS_ImexBDF1::get_tau(
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

double PLocAssem_VMS_NS_ImexBDF1::get_DC(
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

void PLocAssem_VMS_NS_ImexBDF1::Assem_Residual(
    const double &time, const double &dt,
    const double * const &sol_0,
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
    // variables at t_n+1
    double u = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v = 0.0, v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w = 0.0, w_x = 0.0, w_y = 0.0, w_z = 0.0;
    double p = 0.0;

    // variables at t_n
    double u_0 = 0.0, u_x_0 = 0.0, u_y_0 = 0.0, u_z_0 = 0.0;
    double v_0 = 0.0, v_x_0 = 0.0, v_y_0 = 0.0, v_z_0 = 0.0;
    double w_0 = 0.0, w_x_0 = 0.0, w_y_0 = 0.0, w_z_0 = 0.0;
    double p_x_0 = 0.0, p_y_0 = 0.0, p_z_0 = 0.0;
    double u_xx_0 = 0.0, u_yy_0 = 0.0, u_zz_0 = 0.0;
    double v_xx_0 = 0.0, v_yy_0 = 0.0, v_zz_0 = 0.0;
    double w_xx_0 = 0.0, w_yy_0 = 0.0, w_zz_0 = 0.0;

    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_3D_R_gradR_LaplacianR( qua, &R[0], &dR_dx[0], 
        &dR_dy[0], &dR_dz[0], &d2R_dxx[0], &d2R_dyy[0], &d2R_dzz[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int ii4 = 4 * ii;
      // variables at t_n+1 
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

      // variables at t_n
      u_0 += sol_0[ii4+1] * R[ii];
      v_0 += sol_0[ii4+2] * R[ii];
      w_0 += sol_0[ii4+3] * R[ii];

      u_x_0 += sol_0[ii4+1] * dR_dx[ii];
      v_x_0 += sol_0[ii4+2] * dR_dx[ii];
      w_x_0 += sol_0[ii4+3] * dR_dx[ii];
      p_x_0 += sol_0[ii4+0] * dR_dx[ii];

      u_y_0 += sol_0[ii4+1] * dR_dy[ii];
      v_y_0 += sol_0[ii4+2] * dR_dy[ii];
      w_y_0 += sol_0[ii4+3] * dR_dy[ii];
      p_y_0 += sol_0[ii4+0] * dR_dy[ii];

      u_z_0 += sol_0[ii4+1] * dR_dz[ii];
      v_z_0 += sol_0[ii4+2] * dR_dz[ii];
      w_z_0 += sol_0[ii4+3] * dR_dz[ii];
      p_z_0 += sol_0[ii4+0] * dR_dz[ii];

      u_xx_0 += sol_0[ii4+1] * d2R_dxx[ii];
      v_xx_0 += sol_0[ii4+2] * d2R_dxx[ii];
      w_xx_0 += sol_0[ii4+3] * d2R_dxx[ii];
      
      u_yy_0 += sol_0[ii4+1] * d2R_dyy[ii];
      v_yy_0 += sol_0[ii4+2] * d2R_dyy[ii];
      w_yy_0 += sol_0[ii4+3] * d2R_dyy[ii];
      
      u_zz_0 += sol_0[ii4+1] * d2R_dzz[ii];
      v_zz_0 += sol_0[ii4+2] * d2R_dzz[ii];
      w_zz_0 += sol_0[ii4+3] * d2R_dzz[ii];

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];

    }
    // 1st-order BDF scheme  
    const double alpha_sig = 1;

    const double u_BDF = u_0;
    const double v_BDF = v_0;
    const double w_BDF = w_0;

    // Newton-Gregory Polynominal
    const double u_NG = u_0;
    const double v_NG = v_0;
    const double w_NG = w_0;
    //const double p_NG = p_0;

    const double u_x_NG = u_x_0;
    const double v_x_NG = v_x_0;
    const double w_x_NG = w_x_0;
    const double p_x_NG = p_x_0;

    const double u_y_NG = u_y_0;
    const double v_y_NG = v_y_0;
    const double w_y_NG = w_y_0;
    const double p_y_NG = p_y_0;

    const double u_z_NG = u_z_0;
    const double v_z_NG = v_z_0;
    const double w_z_NG = w_z_0;
    const double p_z_NG = p_z_0;

    const double u_xx_NG = u_xx_0;
    const double v_xx_NG = v_xx_0;
    const double w_xx_NG = w_xx_0;

    const double u_yy_NG = u_yy_0;
    const double v_yy_NG = v_yy_0;
    const double w_yy_NG = w_yy_0;

    const double u_zz_NG = u_zz_0;
    const double v_zz_NG = v_zz_0;
    const double w_zz_NG = w_zz_0;

    // Compute the laplacian of vector
    const double u_lap_NG = u_xx_NG + u_yy_NG + u_zz_NG;
    const double v_lap_NG = v_xx_NG + v_yy_NG + v_zz_NG;
    const double w_lap_NG = w_xx_NG + w_yy_NG + w_zz_NG;

    // Get the tau_m and tau_c
    const auto dxi_dx = element->get_invJacobian(qua);

    const std::array<double, 2> tau = get_tau( dt, dxi_dx, u_NG, v_NG, w_NG );
    const double tau_m = tau[0];
    const double tau_c = tau[1];

    const double tau_m_2 = tau_m * tau_m;

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    // Get the body force
    const Vector_3 f_body = get_f( coor, curr );

    // Compute the residue 
    const double rx_hat = rho0 * ( (alpha_sig * u_NG - u_BDF)/dt 
                        + u_NG * u_x_NG + v_NG * u_y_NG + w_NG * u_z_NG - f_body.x() ) 
                        + p_x_NG - vis_mu * u_lap_NG; 

    const double ry_hat = rho0 * ( (alpha_sig * v_NG - v_BDF)/dt 
                        + u_NG * v_x_NG + v_NG * v_y_NG + w_NG * v_z_NG - f_body.y() ) 
                        + p_y_NG - vis_mu * v_lap_NG; 

    const double rz_hat = rho0 * ( (alpha_sig * w_NG - w_BDF)/dt 
                        + u_NG * w_x_NG + v_NG * w_y_NG + w_NG * w_z_NG - f_body.z() ) 
                        + p_z_NG - vis_mu * w_lap_NG;  

    const double div_vel = u_x + v_y + w_z;

    const double r_hat_dot_gradu_NG = rx_hat * u_x_NG + ry_hat * u_y_NG + rz_hat * u_z_NG;
    const double r_hat_dot_gradv_NG = rx_hat * v_x_NG + ry_hat * v_y_NG + rz_hat * v_z_NG;
    const double r_hat_dot_gradw_NG = rx_hat * w_x_NG + ry_hat * w_y_NG + rz_hat * w_z_NG;
    
    // Get the Discontinuity Capturing tau
    //const double tau_dc = get_DC( dxi_dx, u_prime, v_prime, w_prime );

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];
      
      const double velo_NG_dot_gradR  =  u_NG * NA_x    + v_NG * NA_y    + w_NG * NA_z;
      const double r_hat_dot_gradR    =  rx_hat * NA_x  + ry_hat * NA_y  + rz_hat * NA_z;
      
      Residual[4*A] += gwts * ( NA * div_vel + tau_m * r_hat_dot_gradR );

      Residual[4*A+1] += gwts * ( NA * rho0 / dt * (alpha_sig * u - u_BDF) 
          + NA * rho0 * (u_NG * u_x_NG + v_NG * u_y_NG + w_NG * u_z_NG)
          - NA_x * p
          + NA_x * two_mu * u_x
          + NA_y * vis_mu * (u_y + v_x)
          + NA_z * vis_mu * (u_z + w_x)
          - NA * rho0 * f_body.x()
          + velo_NG_dot_gradR * rho0 * tau_m * rx_hat
          - NA * rho0 * tau_m * r_hat_dot_gradu_NG
          - r_hat_dot_gradR * rho0 * tau_m_2 * rx_hat
          + NA_x * tau_c * div_vel );

      Residual[4*A+2] += gwts * ( NA * rho0 / dt * (alpha_sig * v - v_BDF) 
          + NA * rho0 * (u_NG * v_x_NG + v_NG * v_y_NG + w_NG * v_z_NG)
          - NA_y * p
          + NA_x * vis_mu * (v_x + u_y)
          + NA_y * two_mu * v_y 
          + NA_z * vis_mu * (v_z + w_y)
          - NA * rho0 * f_body.y()
          + velo_NG_dot_gradR * rho0 * tau_m * ry_hat
          - NA * rho0 * tau_m * r_hat_dot_gradv_NG
          - r_hat_dot_gradR * rho0 * tau_m_2 * ry_hat 
          + NA_y * tau_c * div_vel );

      Residual[4*A+3] += gwts * ( NA * rho0 / dt * (alpha_sig * w - w_BDF)
          + NA * rho0 * (u_NG * w_x_NG + v_NG * w_y_NG + w_NG * w_z_NG)
          - NA_z * p
          + NA_x * vis_mu * (w_x + u_z)
          + NA_y * vis_mu * (w_y + v_z)
          + NA_z * two_mu * w_z
          - NA * rho0 * f_body.z()
          + velo_NG_dot_gradR * rho0 * tau_m * rz_hat
          - NA * rho0 * tau_m * r_hat_dot_gradw_NG
          - r_hat_dot_gradR * rho0 * tau_m_2 * rz_hat 
          + NA_z * tau_c * div_vel );
    } // A-loop
  } //qua-loop
}

void PLocAssem_VMS_NS_ImexBDF1::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &sol_0,
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

  Zero_Tangent_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);
  std::vector<double> d2R_dxx(nLocBas, 0.0), d2R_dyy(nLocBas, 0.0), d2R_dzz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    // variables at t_n+1
    double u = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v = 0.0, v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w = 0.0, w_x = 0.0, w_y = 0.0, w_z = 0.0;
    double p = 0.0;

    // variables at t_n
    double u_0 = 0.0, u_x_0 = 0.0, u_y_0 = 0.0, u_z_0 = 0.0;
    double v_0 = 0.0, v_x_0 = 0.0, v_y_0 = 0.0, v_z_0 = 0.0;
    double w_0 = 0.0, w_x_0 = 0.0, w_y_0 = 0.0, w_z_0 = 0.0;
    double p_x_0 = 0.0, p_y_0 = 0.0, p_z_0 = 0.0;
    double u_xx_0 = 0.0, u_yy_0 = 0.0, u_zz_0 = 0.0;
    double v_xx_0 = 0.0, v_yy_0 = 0.0, v_zz_0 = 0.0;
    double w_xx_0 = 0.0, w_yy_0 = 0.0, w_zz_0 = 0.0;

    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_3D_R_gradR_LaplacianR( qua, &R[0], &dR_dx[0], 
        &dR_dy[0], &dR_dz[0], &d2R_dxx[0], &d2R_dyy[0], &d2R_dzz[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int ii4 = 4 * ii;
      // variables at t_n+1 
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

      // variables at t_n
      u_0 += sol_0[ii4+1] * R[ii];
      v_0 += sol_0[ii4+2] * R[ii];
      w_0 += sol_0[ii4+3] * R[ii];

      u_x_0 += sol_0[ii4+1] * dR_dx[ii];
      v_x_0 += sol_0[ii4+2] * dR_dx[ii];
      w_x_0 += sol_0[ii4+3] * dR_dx[ii];
      p_x_0 += sol_0[ii4+0] * dR_dx[ii];

      u_y_0 += sol_0[ii4+1] * dR_dy[ii];
      v_y_0 += sol_0[ii4+2] * dR_dy[ii];
      w_y_0 += sol_0[ii4+3] * dR_dy[ii];
      p_y_0 += sol_0[ii4+0] * dR_dy[ii];

      u_z_0 += sol_0[ii4+1] * dR_dz[ii];
      v_z_0 += sol_0[ii4+2] * dR_dz[ii];
      w_z_0 += sol_0[ii4+3] * dR_dz[ii];
      p_z_0 += sol_0[ii4+0] * dR_dz[ii];

      u_xx_0 += sol_0[ii4+1] * d2R_dxx[ii];
      v_xx_0 += sol_0[ii4+2] * d2R_dxx[ii];
      w_xx_0 += sol_0[ii4+3] * d2R_dxx[ii];

      u_yy_0 += sol_0[ii4+1] * d2R_dyy[ii];
      v_yy_0 += sol_0[ii4+2] * d2R_dyy[ii];
      w_yy_0 += sol_0[ii4+3] * d2R_dyy[ii];

      u_zz_0 += sol_0[ii4+1] * d2R_dzz[ii];
      v_zz_0 += sol_0[ii4+2] * d2R_dzz[ii];
      w_zz_0 += sol_0[ii4+3] * d2R_dzz[ii];

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];

    }
    // 1st-order BDF scheme  
    const double alpha_sig = 1.0;

    const double u_BDF = u_0;
    const double v_BDF = v_0;
    const double w_BDF = w_0;

    // Newton-Gregory Polynominal
    const double u_NG = u_0;
    const double v_NG = v_0;
    const double w_NG = w_0;
    //const double p_NG = p_0;

    const double u_x_NG = u_x_0;
    const double v_x_NG = v_x_0;
    const double w_x_NG = w_x_0;
    const double p_x_NG = p_x_0;

    const double u_y_NG = u_y_0;
    const double v_y_NG = v_y_0;
    const double w_y_NG = w_y_0;
    const double p_y_NG = p_y_0;

    const double u_z_NG = u_z_0;
    const double v_z_NG = v_z_0;
    const double w_z_NG = w_z_0;
    const double p_z_NG = p_z_0;

    const double u_xx_NG = u_xx_0;
    const double v_xx_NG = v_xx_0;
    const double w_xx_NG = w_xx_0;

    const double u_yy_NG = u_yy_0;
    const double v_yy_NG = v_yy_0;
    const double w_yy_NG = w_yy_0;

    const double u_zz_NG = u_zz_0;
    const double v_zz_NG = v_zz_0;
    const double w_zz_NG = w_zz_0;

    // Compute the laplacian of vector
    const double u_lap_NG = u_xx_NG + u_yy_NG + u_zz_NG;
    const double v_lap_NG = v_xx_NG + v_yy_NG + v_zz_NG;
    const double w_lap_NG = w_xx_NG + w_yy_NG + w_zz_NG;

    // Get the tau_m and tau_c
    const auto dxi_dx = element->get_invJacobian(qua);

    const std::array<double, 2> tau = get_tau( dt, dxi_dx, u_NG, v_NG, w_NG );
    const double tau_m = tau[0];
    const double tau_c = tau[1];

    const double tau_m_2 = tau_m * tau_m;

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    // Get the body force
    const Vector_3 f_body = get_f( coor, curr );

    // Compute the residue 
    const double rx_hat = rho0 * ( (alpha_sig * u_NG - u_BDF)/dt 
                        + u_NG * u_x_NG + v_NG * u_y_NG + w_NG * u_z_NG - f_body.x() ) 
                        + p_x_NG - vis_mu * u_lap_NG; 

    const double ry_hat = rho0 * ( (alpha_sig * v_NG - v_BDF)/dt 
                        + u_NG * v_x_NG + v_NG * v_y_NG + w_NG * v_z_NG - f_body.y() ) 
                        + p_y_NG - vis_mu * v_lap_NG; 

    const double rz_hat = rho0 * ( (alpha_sig * w_NG - w_BDF)/dt 
                        + u_NG * w_x_NG + v_NG * w_y_NG + w_NG * w_z_NG - f_body.z() ) 
                        + p_z_NG - vis_mu * w_lap_NG;  

    const double div_vel = u_x + v_y + w_z;

    const double r_hat_dot_gradu_NG = rx_hat * u_x_NG + ry_hat * u_y_NG + rz_hat * u_z_NG;
    const double r_hat_dot_gradv_NG = rx_hat * v_x_NG + ry_hat * v_y_NG + rz_hat * v_z_NG;
    const double r_hat_dot_gradw_NG = rx_hat * w_x_NG + ry_hat * w_y_NG + rz_hat * w_z_NG;
    
    // Get the Discontinuity Capturing tau
    //const double tau_dc = get_DC( dxi_dx, u_prime, v_prime, w_prime );

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];
      
      const double velo_NG_dot_gradR  =  u_NG * NA_x    + v_NG * NA_y    + w_NG * NA_z;
      const double r_hat_dot_gradR    =  rx_hat * NA_x  + ry_hat * NA_y  + rz_hat * NA_z;
      
      Residual[4*A] += gwts * ( NA * div_vel + tau_m * r_hat_dot_gradR );

      Residual[4*A+1] += gwts * ( NA * rho0 / dt * (alpha_sig * u - u_BDF) 
          + NA * rho0 * (u_NG * u_x_NG + v_NG * u_y_NG + w_NG * u_z_NG)
          - NA_x * p
          + NA_x * two_mu * u_x
          + NA_y * vis_mu * (u_y + v_x)
          + NA_z * vis_mu * (u_z + w_x)
          - NA * rho0 * f_body.x()
          + velo_NG_dot_gradR * rho0 * tau_m * rx_hat
          - NA * rho0 * tau_m * r_hat_dot_gradu_NG
          - r_hat_dot_gradR * rho0 * tau_m_2 * rx_hat
          + NA_x * tau_c * div_vel );

      Residual[4*A+2] += gwts * ( NA * rho0 / dt * (alpha_sig * v - v_BDF) 
          + NA * rho0 * (u_NG * v_x_NG + v_NG * v_y_NG + w_NG * v_z_NG)
          - NA_y * p
          + NA_x * vis_mu * (v_x + u_y)
          + NA_y * two_mu * v_y 
          + NA_z * vis_mu * (v_z + w_y)
          - NA * rho0 * f_body.y()
          + velo_NG_dot_gradR * rho0 * tau_m * ry_hat
          - NA * rho0 * tau_m * r_hat_dot_gradv_NG
          - r_hat_dot_gradR * rho0 * tau_m_2 * ry_hat 
          + NA_y * tau_c * div_vel );

      Residual[4*A+3] += gwts * ( NA * rho0 / dt * (alpha_sig * w - w_BDF)
          + NA * rho0 * (u_NG * w_x_NG + v_NG * w_y_NG + w_NG * w_z_NG)
          - NA_z * p
          + NA_x * vis_mu * (w_x + u_z)
          + NA_y * vis_mu * (w_y + v_z)
          + NA_z * two_mu * w_z
          - NA * rho0 * f_body.z()
          + velo_NG_dot_gradR * rho0 * tau_m * rz_hat
          - NA * rho0 * tau_m * r_hat_dot_gradw_NG
          - r_hat_dot_gradR * rho0 * tau_m_2 * rz_hat 
          + NA_z * tau_c * div_vel );
      
      for(int B=0; B<nLocBas; ++B)
      {
        const double NB = R[B], NB_x = dR_dx[B], NB_y = dR_dy[B], NB_z = dR_dz[B];

        const double NANB  = NA*NB;
        const double NAxNB = NA_x*NB, NAxNBx = NA_x*NB_x, NAxNBy = NA_x*NB_y, NAxNBz = NA_x*NB_z;
        const double NAyNB = NA_y*NB, NAyNBx = NA_y*NB_x, NAyNBy = NA_y*NB_y, NAyNBz = NA_y*NB_z;
        const double NAzNB = NA_z*NB, NAzNBx = NA_z*NB_x, NAzNBy = NA_z*NB_y, NAzNBz = NA_z*NB_z;

        // Continuity equation with respect to p, u, v, w
        Tangent[16*nLocBas*A+4*B]   += 0.0;
        Tangent[16*nLocBas*A+4*B+1] += gwts *  NA * NB_x;
        Tangent[16*nLocBas*A+4*B+2] += gwts *  NA * NB_y;
        Tangent[16*nLocBas*A+4*B+3] += gwts *  NA * NB_z;

        // Momentum-x with respect to p, u, v, w
        Tangent[4*nLocBas*(4*A+1)+4*B] += gwts * ( (-1.0) * NAxNB);
        
        Tangent[4*nLocBas*(4*A+1)+4*B+1] += gwts * ( rho0 /dt * alpha_sig  * NANB
            + vis_mu * (2.0*NAxNBx + NAyNBy + NAzNBz)
            + tau_c * NAxNBx );

        Tangent[4*nLocBas*(4*A+1)+4*B+2] += gwts * ( vis_mu * NAyNBx 
            + tau_c * NAxNBy );

        Tangent[4*nLocBas*(4*A+1)+4*B+3] += gwts * ( vis_mu * NAzNBx 
            + tau_c * NAxNBz );

        // Momentum-y with respect to p, u, v, w
        Tangent[4*nLocBas*(4*A+2)+4*B] += gwts * (-1.0) * NAyNB;
          
        Tangent[4*nLocBas*(4*A+2)+4*B+1] += gwts * ( vis_mu * NAxNBy
            + tau_c * NAyNBx );

         Tangent[4*nLocBas*(4*A+2)+4*B+2] += gwts * ( rho0 /dt * alpha_sig  * NANB
            + vis_mu * (NAxNBx + 2.0*NAyNBy + NAzNBz)
            + tau_c * NAyNBy );  

        Tangent[4*nLocBas*(4*A+2)+4*B+3] += gwts * ( vis_mu * NAzNBy
            + tau_c * NAyNBz );

        // Momentum-z with respect to p, u, v, w
        Tangent[4*nLocBas*(4*A+3)+4*B] += gwts *  (-1.0) * NAzNB;

        Tangent[4*nLocBas*(4*A+3)+4*B+1] += gwts * ( vis_mu * NAxNBz
            + tau_c * NAzNBx );

        Tangent[4*nLocBas*(4*A+3)+4*B+2] += gwts * ( vis_mu * NAyNBz
            + tau_c * NAzNBy );

        Tangent[4*nLocBas*(4*A+3)+4*B+3] += gwts * ( rho0 /dt * alpha_sig  * NANB
            + vis_mu * (NAxNBx + NAyNBy + 2.0*NAzNBz)
            + tau_c * NAzNBz );

      } // B-loop
    } // A-loop
  } //qua-loop
    // ----------------------------------------------------------------
    // The local `stiffness' matrix 
    //            K[p][q] = Sub_Tan[4*ii + jj][A*nLocBas+B],
    // where p = 4*A+ii, q = 4*B+jj, and K has 4*nLocBas rows/columns.
    // Tangent is a 1D vector storing K by rows:
    // Tangent[4*nLocBas*p + q] = K[p][q] = Sub_Tan[4*ii+jj][A*nLocBas+B]
  // ----------------------------------------------------------------  
}

void PLocAssem_VMS_NS_ImexBDF1::Assem_Mass_Residual(
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

void PLocAssem_VMS_NS_ImexBDF1::Assem_Residual_EBC(
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

  const double curr = time + dt;

  Zero_Residual();

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
      Residual[4*A+1] -= surface_area * quad -> get_qw(qua) * R[A] * traction.x();
      Residual[4*A+2] -= surface_area * quad -> get_qw(qua) * R[A] * traction.y();
      Residual[4*A+3] -= surface_area * quad -> get_qw(qua) * R[A] * traction.z();
    }
  }
}

void PLocAssem_VMS_NS_ImexBDF1::Assem_Residual_EBC_Resistance(
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

  Zero_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    for(int A=0; A<snLocBas; ++A)
    {
      Residual[4*A+1] += surface_area * quad -> get_qw(qua) * R[A] * n_out.x() * val;
      Residual[4*A+2] += surface_area * quad -> get_qw(qua) * R[A] * n_out.y() * val;
      Residual[4*A+3] += surface_area * quad -> get_qw(qua) * R[A] * n_out.z() * val;
    }
  }
}

void PLocAssem_VMS_NS_ImexBDF1::Assem_Residual_BackFlowStab(
    const double * const &sol_0,
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

    double surface_area, factor;

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    double u_0 = 0.0, v_0 = 0.0, w_0 = 0.0;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      const int ii4 = ii * 4;
      u_0 += sol_0[ii4+1] * R[ii];
      v_0 += sol_0[ii4+2] * R[ii];
      w_0 += sol_0[ii4+3] * R[ii];
    }

    // Newton-Gregory Polynominal
    const double u_NG = u_0;
    const double v_NG = v_0;
    const double w_NG = w_0;

    const double temp = u_NG * n_out.x() + v_NG * n_out.y() + w_NG * n_out.z();

    if(temp < 0.0) factor = temp * rho0 * beta;
    else factor = 0.0;

    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual[4*A+1] -= surface_area * quad -> get_qw(qua) * R[A] * factor * u_NG;
      sur_Residual[4*A+2] -= surface_area * quad -> get_qw(qua) * R[A] * factor * v_NG;
      sur_Residual[4*A+3] -= surface_area * quad -> get_qw(qua) * R[A] * factor * w_NG;
    }
  }
}

double PLocAssem_VMS_NS_ImexBDF1::get_flowrate( const double * const &sol,
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

    double u = 0.0, v = 0.0, w = 0.0;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      const int ii4 = ii*4;
      u += sol[ii4+1] * R[ii];
      v += sol[ii4+2] * R[ii];
      w += sol[ii4+3] * R[ii];
    }

    flrate += surface_area * quad->get_qw(qua) * ( u * n_out.x() + v * n_out.y() + w * n_out.z() );
  }

  return flrate;
}

void PLocAssem_VMS_NS_ImexBDF1::get_pressure_area( 
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
