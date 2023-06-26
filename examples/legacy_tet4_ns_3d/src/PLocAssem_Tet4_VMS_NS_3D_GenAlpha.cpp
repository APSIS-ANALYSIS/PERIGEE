#include "PLocAssem_Tet4_VMS_NS_3D_GenAlpha.hpp"

PLocAssem_Tet4_VMS_NS_3D_GenAlpha::PLocAssem_Tet4_VMS_NS_3D_GenAlpha(
    const TimeMethod_GenAlpha * const &tm_gAlpha,
    const int &in_nlocbas, const int &in_nqp,
    const int &in_snlocbas )
: rho0( 1.056 ), nu( 0.035 / rho0 ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()),
  num_ebc_fun(0), nLocBas(4), dof_per_node(4), vec_size(16),
  nqp(in_nqp), snLocBas(in_snlocbas), CI(36.0), CT(4.0)
{
  // Allocate the memory for Tangent and Residual
  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  Zero_Tangent_Residual();

  if( num_ebc_fun == 0 ) flist = NULL;
  else flist = new locassem_tet4_vms_ns_funs [num_ebc_fun];

  // Users have to specify the functions for ebc integration in correct order
  //flist[0] = &PLocAssem_Tet4_VMS_NS_3D_GenAlpha::get_H1;
  //flist[1] = &PLocAssem_Tet4_VMS_NS_3D_GenAlpha::get_H2;

  print_info();
}



PLocAssem_Tet4_VMS_NS_3D_GenAlpha::~PLocAssem_Tet4_VMS_NS_3D_GenAlpha()
{
  delete [] Tangent; Tangent = NULL; delete [] Residual; Residual = NULL;
  if(num_ebc_fun > 0) delete [] flist;
}


void PLocAssem_Tet4_VMS_NS_3D_GenAlpha::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Three-dimensional Incompressible Navier-Stokes: \n");
  PetscPrintf(PETSC_COMM_WORLD, "  FEM: 4-node Tetrahedral \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Spatial: Variational Multiscale Methods \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Temporal: Generalized-alpha Method \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Density rho = %e \n", rho0);
  PetscPrintf(PETSC_COMM_WORLD, "  Dynamic Viscosity mu = %e \n", nu * rho0);
  PetscPrintf(PETSC_COMM_WORLD, "  Kienmatic Viscosity nu = %e \n", nu);
  PetscPrintf(PETSC_COMM_WORLD, "  Stabilization para CI = %e \n", CI);
  PetscPrintf(PETSC_COMM_WORLD, "  Stabilization para CT = %e \n", CT);
  PetscPrintf(PETSC_COMM_WORLD, "  Consistent tangent matrix used. \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Nonlinear quadratic term is in advective form. \n");
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PLocAssem_Tet4_VMS_NS_3D_GenAlpha::get_metric( 
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


void PLocAssem_Tet4_VMS_NS_3D_GenAlpha::get_tau( 
    double &tau_m_qua, double &tau_c_qua,
    const double &dt, const double * const &dxi_dx,
    const double &u, const double &v, const double &w ) const
{
  double G11, G12, G13, G22, G23, G33;
  get_metric( dxi_dx, G11, G12, G13, G22, G23, G33 );
  
  const double GdG = G11 * G11 + 2.0 * G12 * G12 + 2.0 * G13 * G13
    + G22 * G22 + 2.0 * G23 * G23 + G33 * G33;

  const double uGu = G11 * u * u + 2.0 * G12 * u * v + 2.0 * G13 * u * w
    + G22 * v * v + 2.0 * G23 * v * w + G33 * w * w;

  const double g_dot_g = G11 + G22 + G33;

  const double denom_m = CT / (dt*dt) + uGu + CI * nu * nu * GdG;

  tau_m_qua = 1.0 / sqrt(denom_m);

  const double denom_c = tau_m_qua * g_dot_g;

  tau_c_qua = 1.0 / denom_c;
}


void PLocAssem_Tet4_VMS_NS_3D_GenAlpha::get_DC( double &dc_tau,
    const double * const &dxi_dx,
    const double &u, const double &v, const double &w ) const
{
  //double G11, G12, G13, G22, G23, G33;
  //get_metric( dxi_dx, G11, G12, G13, G22, G23, G33 );

  //dc_tau = G11 * u * u + 2.0 * G12 * u * v + 2.0 * G13 * u * w + G22 * v * v
  //  + 2.0 * G23 * v * w + G33 * w * w;

  //if(dc_tau > 1.0e-10) dc_tau = std::pow(dc_tau, -0.5);
  
  dc_tau = 0.0;
}


void PLocAssem_Tet4_VMS_NS_3D_GenAlpha::Assem_Residual(
    const double &time, const double &dt,
    const double * const &velo,
    const double * const &disp,
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

  // Momentum residual and their derivatives
  double rx, ry, rz;

  double u_prime, v_prime, w_prime;

  double tau_m, tau_c, tau_dc;

  const double two_nu = 2.0 * nu;
  double NA, NA_x, NA_y, NA_z;
  double velo_dot_gradR, div_vel, r_dot_gradR;
  double tau_m_2;

  double velo_prime_dot_gradR; // v' dot gradNA

  const double curr = time + alpha_f * dt;

  Zero_Residual();

  for(qua=0; qua<nqp; ++qua)
  {
    u = 0.0; u_t = 0.0; u_x = 0.0; u_y = 0.0; u_z = 0.0;
    v = 0.0; v_t = 0.0; v_x = 0.0; v_y = 0.0; v_z = 0.0;
    w = 0.0; w_t = 0.0; w_x = 0.0; w_y = 0.0; w_z = 0.0;
    p = 0.0; coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    p_x = 0.0; p_y = 0.0; p_z = 0.0;

    element->get_R_gradR( qua, R, dR_dx, dR_dy, dR_dz );
    element->get_invJacobian( qua, dxi_dx );

    for(ii=0; ii<nLocBas; ++ii)
    {
      ii4 = 4 * ii;
      u_t += velo[ii4+1] * R[ii];
      v_t += velo[ii4+2] * R[ii];
      w_t += velo[ii4+3] * R[ii];

      u += disp[ii4+1] * R[ii];
      v += disp[ii4+2] * R[ii];
      w += disp[ii4+3] * R[ii];
      p += disp[ii4] * R[ii];

      u_x += disp[ii4+1] * dR_dx[ii];
      v_x += disp[ii4+2] * dR_dx[ii];
      w_x += disp[ii4+3] * dR_dx[ii];
      p_x += disp[ii4]   * dR_dx[ii];

      u_y += disp[ii4+1] * dR_dy[ii];
      v_y += disp[ii4+2] * dR_dy[ii];
      w_y += disp[ii4+3] * dR_dy[ii];
      p_y += disp[ii4]   * dR_dy[ii];

      u_z += disp[ii4+1] * dR_dz[ii];
      v_z += disp[ii4+2] * dR_dz[ii];
      w_z += disp[ii4+3] * dR_dz[ii];
      p_z += disp[ii4]   * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    get_tau(tau_m, tau_c, dt, dxi_dx, u, v, w);

    tau_m_2 = tau_m * tau_m;

    gwts = element->get_detJac(qua) * quad->get_qw(qua); 

    get_f(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    rx = u_t + u_x * u + u_y * v + u_z * w + p_x - f1;
    ry = v_t + v_x * u + v_y * v + v_z * w + p_y - f2;
    rz = w_t + w_x * u + w_y * v + w_z * w + p_z - f3;

    div_vel = u_x + v_y + w_z;

    // DC terms
    u_prime = -1.0 * tau_m * rx;
    v_prime = -1.0 * tau_m * ry;
    w_prime = -1.0 * tau_m * rz;

    get_DC( tau_dc, dxi_dx, u_prime, v_prime, w_prime );

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];
      velo_dot_gradR = NA_x * u + NA_y * v + NA_z * w;
      r_dot_gradR = NA_x * rx + NA_y * ry + NA_z * rz;

      velo_prime_dot_gradR = NA_x * u_prime + NA_y * v_prime + NA_z * w_prime;

      Residual[4*A] += gwts * ( NA * div_vel + tau_m * r_dot_gradR );

      Residual[4*A+1] += gwts * (NA * u_t
          + NA * (u * u_x + v * u_y + w * u_z)
          - NA_x * p
          + NA_x * two_nu * u_x
          + NA_y * nu * (u_y + v_x)
          + NA_z * nu * (u_z + w_x)
          + velo_dot_gradR * tau_m * rx
          + r_dot_gradR * u * tau_m
          + NA_x * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * rx
          + velo_prime_dot_gradR * tau_dc
          * (u_prime * u_x + v_prime * u_y + w_prime * u_z)
          - NA * f1 );

      Residual[4*A+2] += gwts * (NA * v_t
          + NA * (u * v_x + v * v_y + w * v_z)
          - NA_y * p
          + NA_x * nu * (u_y + v_x)
          + NA_y * two_nu * v_y
          + NA_z * nu * (v_z + w_y)
          + velo_dot_gradR * tau_m * ry
          + r_dot_gradR * v * tau_m
          + NA_y * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * ry
          + velo_prime_dot_gradR * tau_dc
          * (u_prime * v_x + v_prime * v_y + w_prime * v_z)
          - NA * f2 );

      Residual[4*A+3] += gwts * (NA * w_t
          + NA * (u * w_x + v * w_y + w * w_z)
          - NA_z * p
          + NA_x * nu * (u_z + w_x)
          + NA_y * nu * (w_y + v_z)
          + NA_z * two_nu * w_z
          + velo_dot_gradR * tau_m * rz
          + r_dot_gradR * w * tau_m
          + NA_z * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * rz
          + velo_prime_dot_gradR * tau_dc
          * (u_prime * w_x + v_prime * w_y + w_prime * w_z)
          - NA * f3 );
    }
  }
}


void PLocAssem_Tet4_VMS_NS_3D_GenAlpha::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &velo,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  int ii, qua, A, ii4, B, jj, index;
  double u, u_t, u_x, u_y, u_z;
  double v, v_t, v_x, v_y, v_z;
  double w, w_t, w_x, w_y, w_z;
  double p, f1, f2, f3, p_x, p_y, p_z;
  double gwts, coor_x, coor_y, coor_z;

  double rx, ry, rz; // Momentum residual

  double tau_m, tau_c, tau_dc;

  double u_prime, v_prime, w_prime;
  double velo_prime_dot_gradR, velo_prime_dot_gradNB;

  const double two_nu = 2.0 * nu;
  double NA, NA_x, NA_y, NA_z;
  double velo_dot_gradR, div_vel, r_dot_gradR;
  double tau_m_2;

  double NB, NB_x, NB_y, NB_z;
  double NANB, NAxNB, NAyNB, NAzNB;
  double NANBx, NAxNBx, NAyNBx, NAzNBx;
  double NANBy, NAxNBy, NAyNBy, NAzNBy;
  double NANBz, NAxNBz, NAyNBz, NAzNBz;

  double drx_du_B, drx_dv_B, drx_dw_B, drx_dp_B;
  double dry_du_B, dry_dv_B, dry_dw_B, dry_dp_B;
  double drz_du_B, drz_dv_B, drz_dw_B, drz_dp_B;

  double velo_dot_gradNB;

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

    element->get_R_gradR( qua, R, dR_dx, dR_dy, dR_dz );
    element->get_invJacobian( qua, dxi_dx );

    for(ii=0; ii<nLocBas; ++ii)
    {
      ii4 = 4 * ii;
      u_t += velo[ii4+1]   * R[ii];
      v_t += velo[ii4+2] * R[ii];
      w_t += velo[ii4+3] * R[ii];

      u += disp[ii4+1] * R[ii];
      v += disp[ii4+2] * R[ii];
      w += disp[ii4+3] * R[ii];
      p += disp[ii4] * R[ii];

      u_x += disp[ii4+1] * dR_dx[ii];
      v_x += disp[ii4+2] * dR_dx[ii];
      w_x += disp[ii4+3] * dR_dx[ii];
      p_x += disp[ii4]   * dR_dx[ii];

      u_y += disp[ii4+1] * dR_dy[ii];
      v_y += disp[ii4+2] * dR_dy[ii];
      w_y += disp[ii4+3] * dR_dy[ii];
      p_y += disp[ii4]   * dR_dy[ii];

      u_z += disp[ii4+1] * dR_dz[ii];
      v_z += disp[ii4+2] * dR_dz[ii];
      w_z += disp[ii4+3] * dR_dz[ii];
      p_z += disp[ii4]   * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    get_tau(tau_m, tau_c, dt, dxi_dx, u, v, w);

    tau_m_2 = tau_m * tau_m;

    gwts = element->get_detJac(qua) * quad->get_qw(qua); 

    get_f(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    rx = u_t + u_x * u + u_y * v + u_z * w + p_x - f1;
    ry = v_t + v_x * u + v_y * v + v_z * w + p_y - f2;
    rz = w_t + w_x * u + w_y * v + w_z * w + p_z - f3;

    div_vel = u_x + v_y + w_z;

    // Prepare for DC terms
    u_prime = -1.0 * tau_m * rx;
    v_prime = -1.0 * tau_m * ry;
    w_prime = -1.0 * tau_m * rz;

    get_DC( tau_dc, dxi_dx, u_prime, v_prime, w_prime );

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];
      velo_dot_gradR = NA_x * u + NA_y * v + NA_z * w;
      r_dot_gradR = NA_x * rx + NA_y * ry + NA_z * rz;

      velo_prime_dot_gradR = NA_x * u_prime + NA_y * v_prime + NA_z * w_prime;

      Residual[4*A] += gwts * ( NA * div_vel + tau_m * r_dot_gradR );

      Residual[4*A+1] += gwts * (NA * u_t
          + NA * (u * u_x + v * u_y + w * u_z)
          - NA_x * p
          + NA_x * two_nu * u_x
          + NA_y * nu * (u_y + v_x)
          + NA_z * nu * (u_z + w_x)
          + velo_dot_gradR * tau_m * rx
          + r_dot_gradR * u * tau_m
          + NA_x * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * rx
          + velo_prime_dot_gradR * tau_dc
          * (u_prime * u_x + v_prime * u_y + w_prime * u_z)
          - NA * f1 );

      Residual[4*A+2] += gwts * (NA * v_t
          + NA * (u * v_x + v * v_y + w * v_z)
          - NA_y * p
          + NA_x * nu * (u_y + v_x)
          + NA_y * two_nu * v_y
          + NA_z * nu * (v_z + w_y)
          + velo_dot_gradR * tau_m * ry
          + r_dot_gradR * v * tau_m
          + NA_y * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * ry
          + velo_prime_dot_gradR * tau_dc
          * (u_prime * v_x + v_prime * v_y + w_prime * v_z)
          - NA * f2 );

      Residual[4*A+3] += gwts * (NA * w_t
          + NA * (u * w_x + v * w_y + w * w_z)
          - NA_z * p
          + NA_x * nu * (u_z + w_x)
          + NA_y * nu * (w_y + v_z)
          + NA_z * two_nu * w_z
          + velo_dot_gradR * tau_m * rz
          + r_dot_gradR * w * tau_m
          + NA_z * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * rz
          + velo_prime_dot_gradR * tau_dc
          * (u_prime * w_x + v_prime * w_y + w_prime * w_z)
          - NA * f3 );

      for(B=0; B<nLocBas; ++B)
      {
        index = B + A * nLocBas;
        NB = R[B]; NB_x = dR_dx[B]; NB_y = dR_dy[B]; NB_z = dR_dz[B];
        velo_dot_gradNB = u * NB_x + v * NB_y + w * NB_z;

        velo_prime_dot_gradNB = u_prime * NB_x + v_prime * NB_y + w_prime * NB_z;

        NANB  = NA*NB; NANBx = NA*NB_x; NANBy = NA*NB_y; NANBz = NA*NB_z;
        NAxNB = NA_x*NB; NAxNBx = NA_x*NB_x; NAxNBy = NA_x*NB_y; NAxNBz = NA_x*NB_z;
        NAyNB = NA_y*NB; NAyNBx = NA_y*NB_x; NAyNBy = NA_y*NB_y; NAyNBz = NA_y*NB_z;
        NAzNB = NA_z*NB; NAzNBx = NA_z*NB_x; NAzNBy = NA_z*NB_y; NAzNBz = NA_z*NB_z;

        drx_du_B = u_x * NB + velo_dot_gradNB;
        drx_dv_B = u_y * NB;
        drx_dw_B = u_z * NB;
        drx_dp_B = NB_x;

        dry_du_B = v_x * NB;
        dry_dv_B = v_y * NB + velo_dot_gradNB;
        dry_dw_B = v_z * NB;
        dry_dp_B = NB_y;

        drz_du_B = w_x * NB;
        drz_dv_B = w_y * NB;
        drz_dw_B = w_z * NB + velo_dot_gradNB;
        drz_dp_B = NB_z;

        // Continuity equation with respect to p, u, v, w
        Sub_Tan[0][index] += gwts * tau_m * (NAxNBx + NAyNBy + NAzNBz);

        Sub_Tan[1][index] += gwts * ( alpha_m * tau_m * NAxNB
            + dd_dv * ( NANBx + tau_m * NA_x * drx_du_B
              + tau_m * NA_y * dry_du_B + tau_m * NA_z * drz_du_B ) );

        Sub_Tan[2][index] += gwts * ( alpha_m * tau_m * NAyNB
            + dd_dv * ( NANBy + tau_m * NA_x * drx_dv_B
              + tau_m * NA_y * dry_dv_B + tau_m * NA_z * drz_dv_B ) );

        Sub_Tan[3][index] += gwts * ( alpha_m * tau_m * NAzNB
            + dd_dv * ( NANBz + tau_m * NA_x * drx_dw_B
              + tau_m * NA_y * dry_dw_B + tau_m * NA_z * drz_dw_B ) );

        // Momentum-x with respect to p, u, v, w
        Sub_Tan[4][index] += gwts * ((-1.0) * NAxNB
            + (velo_dot_gradR + u * NA_x) * tau_m * NB_x
            + tau_m * u * (NAyNBy + NAzNBz)
            - 2.0 * tau_m_2 * rx * NA_x * drx_dp_B
            - tau_m_2 * NA_y * (rx * dry_dp_B + ry * drx_dp_B)
            - tau_m_2 * NA_z * (rx * drz_dp_B + rz * drx_dp_B) );

        Sub_Tan[5][index] += gwts * ( alpha_m * ( NANB
              + tau_m * (2.0 * u * NAxNB + v * NAyNB + w * NAzNB)
              - tau_m_2 * (2.0 * rx * NAxNB + ry * NAyNB + rz * NAzNB ) )
            + dd_dv *  ( NA * velo_dot_gradNB + NANB * u_x
              + nu * (2.0*NAxNBx + NAyNBy + NAzNBz)
              + tau_m * (2.0 * rx * NAxNB + ry * NAyNB + rz * NAzNB)
              + (velo_dot_gradR + u * NA_x) * tau_m * drx_du_B
              + tau_m * u * v_x * NAyNB + tau_m * u * w_x * NAzNB + tau_c * NAxNBx
              - 2.0 * tau_m_2 * rx  * NA_x * drx_du_B
              - tau_m_2 * ry * NA_y * drx_du_B
              - tau_m_2 * rz * NA_z * drx_du_B
              - tau_m_2 * rx * (v_x * NAyNB + w_x * NAzNB)
              + velo_prime_dot_gradR * tau_dc * velo_prime_dot_gradNB ) );

        Sub_Tan[6][index] += gwts * ( 
            alpha_m * ( tau_m * u - tau_m_2 * rx ) * NAyNB
            + dd_dv * ( NANB * u_y + nu * NAyNBx + tau_m * rx * NAyNB
              + (velo_dot_gradR + u * NA_x) * tau_m * drx_dv_B
              + tau_m * u * ( NA_y * dry_dv_B + NA_z * drz_dv_B ) + tau_c * NAxNBy
              - 2.0 * tau_m_2 * rx * NA_x * drx_dv_B
              - tau_m_2 * NA_y * (rx * dry_dv_B + ry * drx_dv_B)
              - tau_m_2 * NA_z * (rx * drz_dv_B + rz * drx_dv_B) ) );

        Sub_Tan[7][index] += gwts * (
            alpha_m * ( tau_m * u - tau_m_2 * rx ) * NAzNB
            + dd_dv * ( NANB * u_z + nu * NAzNBx + tau_m * rx * NAzNB
              + (velo_dot_gradR + u * NA_x) * tau_m * drx_dw_B
              + tau_m * u * ( NA_y * dry_dw_B + NA_z * drz_dw_B ) + tau_c * NAxNBz
              - 2.0 * tau_m_2 * rx * NA_x * drx_dw_B
              - tau_m_2 * NA_y * (rx * dry_dw_B + ry * drx_dw_B)
              - tau_m_2 * NA_z * (rx * drz_dw_B + rz * drx_dw_B) ) );

        // Momentum-y with respect to p u v w
        Sub_Tan[8][index] += gwts * ( (-1.0) * NAyNB
            + tau_m * v * (NAxNBx + NAzNBz)
            + (velo_dot_gradR + NA_y * v) * tau_m * NB_y
            - tau_m_2 * NA_x * (rx * NB_y + ry * NB_x)
            - 2.0 * tau_m_2 * ry * NAyNBy
            - tau_m_2 * NA_z * (ry * NB_z + rz * NB_y) );

        Sub_Tan[9][index] += gwts * (
            alpha_m * (tau_m * v - tau_m_2 * ry) * NAxNB
            + dd_dv * ( NANB * v_x + nu * NAxNBy
              + tau_m * ry * NAxNB
              + NA_x * v * tau_m * drx_du_B
              + (velo_dot_gradR + NA_y * v) * tau_m * dry_du_B
              + NA_z * v * tau_m * drz_du_B
              + tau_c * NAyNBx
              - tau_m_2 * NA_x * (ry * drx_du_B + rx * dry_du_B)
              - 2.0 * tau_m_2 * ry * NA_y * dry_du_B
              - tau_m_2 * NA_z * (ry * drz_du_B + rz * dry_du_B) ) );

        Sub_Tan[10][index] += gwts * (
            alpha_m * ( NANB + (velo_dot_gradR + NA_y * v) * tau_m * NB
              - tau_m_2 * (rx * NAxNB + 2.0 * ry * NAyNB + rz * NAzNB) )
            + dd_dv * ( NA * velo_dot_gradNB + NANB * v_y
              + nu * (NAxNBx + 2.0 * NAyNBy + NAzNBz)
              + tau_m * (rx * NAxNB + 2.0 * ry * NAyNB + rz * NAzNB)
              + tau_m * v * NA_x * drx_dv_B
              + (velo_dot_gradR + v * NA_y) * tau_m * dry_dv_B
              + tau_m * v * NA_z * drz_dv_B
              + tau_c * NAyNBy
              - tau_m_2 * NA_x * (rx * dry_dv_B + ry * drx_dv_B)
              - 2.0 * tau_m_2 * ry * NA_y * dry_dv_B
              - tau_m_2 * NA_z * (ry * drz_dv_B + rz * dry_dv_B)
              + velo_prime_dot_gradR * tau_dc * velo_prime_dot_gradNB ) );

        Sub_Tan[11][index] += gwts * (
            alpha_m * (tau_m * v - tau_m_2 * ry) * NAzNB
            + dd_dv * ( NANB * v_z + nu * NAzNBy
              + tau_m * ry * NAzNB
              + NA_x * v * tau_m * drx_dw_B
              + (velo_dot_gradR  + NA_y * v) * tau_m * dry_dw_B
              + NA_z * v * tau_m * drz_dw_B
              + tau_c * NAyNBz
              - tau_m_2 * NA_x * (rx * dry_dw_B + ry * drx_dw_B)
              - tau_m_2 * 2.0 * ry * NA_y * dry_dw_B
              - tau_m_2 * NA_z * (ry * drz_dw_B + rz * dry_dw_B) ) );

        // Momentum-z with respect to p u v w
        Sub_Tan[12][index] += gwts * ( (-1.0) * NAzNB
            + w * tau_m * NAxNBx
            + w * tau_m * NAyNBy
            + tau_m * (velo_dot_gradR + w * NA_z) * NB_z
            - tau_m_2 * NA_x * (rx * NB_z + rz * NB_x)
            - tau_m_2 * NA_y * (ry * NB_z + rz * NB_y)
            - 2.0 * tau_m_2 * rz * NAzNBz );

        Sub_Tan[13][index] += gwts * (
            alpha_m * ( (tau_m * w - tau_m_2 * rz) * NAxNB )
            + dd_dv * ( NANB * w_x + nu * NAxNBz
              + tau_m * rz * NAxNB
              + tau_m * w * NA_x * drx_du_B
              + tau_m * w * NA_y * dry_du_B
              + (velo_dot_gradR + w * NA_z) * tau_m * drz_du_B
              + tau_c * NAzNBx
              - tau_m_2 * NA_x * (rx * drz_du_B + rz * drx_du_B)
              - tau_m_2 * NA_y * (ry * drz_du_B + rz * dry_du_B)
              - 2.0 * tau_m_2 * rz * NA_z * drz_du_B ) );

        Sub_Tan[14][index] += gwts * (
            alpha_m * (tau_m * w - tau_m_2 * rz) * NAyNB
            + dd_dv * ( NANB * w_y + nu * NAyNBz
              + tau_m * rz * NAyNB
              + tau_m * w * (NA_x * drx_dv_B + NA_y * dry_dv_B)
              + (velo_dot_gradR + w * NA_z) * tau_m * drz_dv_B
              + tau_c * NAzNBy
              - tau_m_2 * NA_x * (rx * drz_dv_B + rz * drx_dv_B)
              - tau_m_2 * NA_y * (ry * drz_dv_B + rz * dry_dv_B)
              - 2.0 * tau_m_2 * rz * NA_z * drz_dv_B ) );

        Sub_Tan[15][index] += gwts * (
            alpha_m * ( NANB + (velo_dot_gradR + w * NA_z) * tau_m * NB
              - tau_m_2 * (rx * NAxNB + ry * NAyNB + 2.0 * rz * NAzNB) )
            + dd_dv * ( NA * velo_dot_gradNB + NANB * w_z
              + nu * (NAxNBx + NAyNBy + 2.0 * NAzNBz)
              + tau_m * (rx * NAxNB + ry * NAyNB + 2.0 * rz * NAzNB)
              + tau_m * w * NA_x * drx_dw_B
              + tau_m * w * NA_y * dry_dw_B
              + tau_m * (velo_dot_gradR + w * NA_z ) * drz_dw_B
              + tau_c * NAzNBz
              - tau_m_2 * NA_x * (rx * drz_dw_B + rz * drx_dw_B)
              - tau_m_2 * NA_y * (ry * drz_dw_B + rz * dry_dw_B)
              - 2.0 * tau_m_2 * NA_z * rz * drz_dw_B
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


void PLocAssem_Tet4_VMS_NS_3D_GenAlpha::Assem_Mass_Residual(
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  int ii, jj, qua, A, B, index;
  double u, u_x, u_y, u_z;
  double v, v_x, v_y, v_z;
  double w, w_x, w_y, w_z;
  double p, f1, f2, f3;
  double gwts, coor_x, coor_y, coor_z;

  double NA, NA_x, NA_y, NA_z;

  const double two_nu = 2.0 * nu;

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
      u += disp[4*ii+1]   * R[ii];
      v += disp[4*ii+2] * R[ii];
      w += disp[4*ii+3] * R[ii];
      p += disp[4*ii] * R[ii];

      u_x += disp[4*ii+1] * dR_dx[ii];
      v_x += disp[4*ii+2] * dR_dx[ii];
      w_x += disp[4*ii+3] * dR_dx[ii];

      u_y += disp[4*ii+1] * dR_dy[ii];
      v_y += disp[4*ii+2] * dR_dy[ii];
      w_y += disp[4*ii+3] * dR_dy[ii];

      u_z += disp[4*ii+1] * dR_dz[ii];
      v_z += disp[4*ii+2] * dR_dz[ii];
      w_z += disp[4*ii+3] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    gwts = element->get_detJac(qua) * quad->get_qw(qua);

    get_f(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      Residual[4*A+1] += gwts * ( NA*(u*u_x + v*u_y + w*u_z) 
          - NA_x * p
          + two_nu * NA_x * u_x
          + nu * NA_y * (u_y + v_x)
          + nu * NA_z * (u_z + w_x)
          - NA * f1 );

      Residual[4*A+2] += gwts * ( NA*(u*v_x + v*v_y + w*v_z) 
          - NA_y * p
          + nu * NA_x * (u_y + v_x)
          + two_nu * NA_y * v_y
          + nu * NA_z * (v_z + w_y)
          - NA * f2 );

      Residual[4*A+3] += gwts * ( NA*(u*w_x + v*w_y + w*w_z) 
          - NA_z * p
          + nu * NA_x * (u_z + w_x)
          + nu * NA_y * (w_y + v_z)
          + two_nu * NA_z * w_z
          - NA * f3 );

      for(B=0; B<nLocBas; ++B)
      {
        index = A * nLocBas + B;
        Sub_Tan[0][index]  += gwts * NA * R[B];
        Sub_Tan[5][index]  += gwts * NA * R[B];
        Sub_Tan[10][index] += gwts * NA * R[B];
        Sub_Tan[15][index] += gwts * NA * R[B];
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


void PLocAssem_Tet4_VMS_NS_3D_GenAlpha::Assem_Residual_EBC(
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
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

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
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
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


// EOF