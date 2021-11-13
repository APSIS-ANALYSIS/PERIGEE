#include "PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha.hpp"

PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha::PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha(
    IMaterialModel * const &in_matmodel,
    const TimeMethod_GenAlpha * const &tm_gAlpha,
    const int &in_nqp )
: rho0( in_matmodel->get_elastic_rho0() ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()),
  num_ebc_fun(0), nLocBas(4), dof_per_node(7), vec_size(16),
  nqp(in_nqp), snLocBas(3)
{
  matmodel = in_matmodel;

  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size]; 
  
  Zero_Tangent_Residual();

  if( num_ebc_fun == 0 ) flist = nullptr;
  else flist = new locassem_vms_seg_ela_fem_funs [num_ebc_fun];

  print_info();
}


PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha::~PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha()
{
  delete [] Tangent; delete [] Residual; Tangent = nullptr; Residual = nullptr;
  if(num_ebc_fun > 0) delete [] flist;
}


void PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("  Three-dimensional Hyper-elastic solid model with VMS mixed formulation, U-kinematic relation, Segregated, FEM formulation: \n");
  SYS_T::commPrint("\t  Spatial: Finite element with VMS stabilization \n");
  SYS_T::commPrint("\t  Temporal: Generalized-alpha method \n");
  SYS_T::commPrint("\t  Solid density rho0 = %e g/cm3\n\n", rho0);
  matmodel->print_info();
  SYS_T::commPrint("----------------------------------------------------------- \n");
}


void PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 0.0;
}


void PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha::Zero_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
}


void PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha::Assem_Estimate()
{
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 1.0;
}


void PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha::get_tau( 
    double &tau_m_qua, double &tau_c_qua,
    const double &dt, const double &Jin,
    const double &dx ) const
{
  const double mu = matmodel->get_elastic_mu();
  const double ka = matmodel->get_elastic_kappa();
  const double c_max = std::pow( rho0 / (ka + 4*mu/3.0), -0.5);
  const double dt_ka = dx / c_max;
  tau_m_qua = 0.000 * dt_ka * Jin / rho0;
  tau_c_qua = 0.000 * dx * c_max * rho0 / Jin;
}


void PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha::Assem_Residual(
    const double &time, const double &dt,
    const double * const &velo,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  double R[4], dR_dx[4], dR_dy[4], dR_dz[4];

  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double h_e = element->get_h( eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  int ii, qua, A, ii7;
  double p, p_t, p_x, p_y, p_z;
  double ux_t, uy_t, uz_t, vx, vy, vz;
  double vx_t, vy_t, vz_t;
  double ux_x, uy_x, uz_x, ux_y, uy_y, uz_y, ux_z, uy_z, uz_z;
  double vx_x, vy_x, vz_x, vx_y, vy_y, vz_y, vx_z, vy_z, vz_z;
  double fx, fy, fz, gwts, coor_x, coor_y, coor_z;
  double NA, NA_x, NA_y, NA_z;

  double Res_Mom[3];
  double Res_Mas;
  double GradNA_invF_ResMom;
  double invFDV_t;

  const double curr = time + alpha_f * dt;

  Zero_Residual();

  for(qua=0; qua < nqp; ++qua)
  {
    p = 0.0; p_t = 0.0; p_x = 0.0; p_y = 0.0; p_z = 0.0;

    ux_t = 0.0; uy_t = 0.0; uz_t = 0.0;
    vx   = 0.0; vy   = 0.0; vz   = 0.0;
    vx_t = 0.0; vy_t = 0.0; vz_t = 0.0;

    ux_x = 0.0; uy_x = 0.0; uz_x = 0.0;
    ux_y = 0.0; uy_y = 0.0; uz_y = 0.0;
    ux_z = 0.0; uy_z = 0.0; uz_z = 0.0;

    vx_x = 0.0; vy_x = 0.0; vz_x = 0.0;
    vx_y = 0.0; vy_y = 0.0; vz_y = 0.0;
    vx_z = 0.0; vy_z = 0.0; vz_z = 0.0;

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    element->get_R_gradR( qua, R, dR_dx, dR_dy, dR_dz );

    for(ii=0; ii<nLocBas; ++ii)
    {
      ii7 = 7*ii;
      p   += disp[ii7+3] * R[ii];
      p_x += disp[ii7+3] * dR_dx[ii];
      p_y += disp[ii7+3] * dR_dy[ii];
      p_z += disp[ii7+3] * dR_dz[ii];

      ux_t += velo[ii7  ] * R[ii];
      uy_t += velo[ii7+1] * R[ii];
      uz_t += velo[ii7+2] * R[ii];
      p_t  += velo[ii7+3] * R[ii];
      vx_t += velo[ii7+4] * R[ii];
      vy_t += velo[ii7+5] * R[ii];
      vz_t += velo[ii7+6] * R[ii];

      vx   += disp[ii7+4] * R[ii];
      vy   += disp[ii7+5] * R[ii];
      vz   += disp[ii7+6] * R[ii];

      ux_x += disp[ii7+0] * dR_dx[ii];
      uy_x += disp[ii7+1] * dR_dx[ii];
      uz_x += disp[ii7+2] * dR_dx[ii];

      ux_y += disp[ii7+0] * dR_dy[ii];
      uy_y += disp[ii7+1] * dR_dy[ii];
      uz_y += disp[ii7+2] * dR_dy[ii];

      ux_z += disp[ii7+0] * dR_dz[ii];
      uy_z += disp[ii7+1] * dR_dz[ii];
      uz_z += disp[ii7+2] * dR_dz[ii];

      vx_x += disp[ii7+4] * dR_dx[ii];
      vy_x += disp[ii7+5] * dR_dx[ii];
      vz_x += disp[ii7+6] * dR_dx[ii];

      vx_y += disp[ii7+4] * dR_dy[ii];
      vy_y += disp[ii7+5] * dR_dy[ii];
      vz_y += disp[ii7+6] * dR_dy[ii];

      vx_z += disp[ii7+4] * dR_dz[ii];
      vy_z += disp[ii7+5] * dR_dz[ii];
      vz_z += disp[ii7+6] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    gwts = element->get_detJac(qua) * quad->get_qw(qua);

    get_f(coor_x, coor_y, coor_z, curr, fx, fy, fz);

    const Matrix_3x3 F( ux_x + 1.0, ux_y, ux_z, uy_x, uy_y + 1.0, uy_z, uz_x, uz_y, uz_z + 1.0 );

    Matrix_3x3 invF( F ); invF.inverse();

    const Matrix_3x3 DVelo( vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z );

    invFDV_t = invF.MatTContraction(DVelo); // invF_Ii V_i,I

    Matrix_3x3 P_iso, S_iso;
    Tensor4_3D AA_iso;

    matmodel->get_PK_FFStiffness(F, P_iso, S_iso, AA_iso);

    const double rho = matmodel->get_rho(p);
    const double mbeta = matmodel->get_beta(p);
    const double detF = F.det();

    // Get stabilization parameters
    double tau_m, tau_c;
    get_tau(tau_m, tau_c, dt, detF, h_e);

    // Residual of momentum equation
    Res_Mom[0] = rho * detF * vx_t;
    Res_Mom[1] = rho * detF * vy_t;
    Res_Mom[2] = rho * detF * vz_t;

    Res_Mom[0] += detF * ( invF(0) * p_x + invF(3) * p_y + invF(6) * p_z );
    Res_Mom[1] += detF * ( invF(1) * p_x + invF(4) * p_y + invF(7) * p_z );
    Res_Mom[2] += detF * ( invF(2) * p_x + invF(5) * p_y + invF(8) * p_z );

    Res_Mom[0] -= rho * detF * fx;
    Res_Mom[1] -= rho * detF * fy;
    Res_Mom[2] -= rho * detF * fz;

    // Residual of mass equation
    Res_Mas = detF * ( mbeta * p_t + invFDV_t );

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      // NA_I invF_Ii 
      double GradNA_invF[3];
      invF.VecMultT( NA_x, NA_y, NA_z, GradNA_invF[0], GradNA_invF[1], GradNA_invF[2] );

      // tau_m stabilization term
      GradNA_invF_ResMom = tau_m * ( GradNA_invF[0] * Res_Mom[0] + GradNA_invF[1] * Res_Mom[1] + GradNA_invF[2] * Res_Mom[2] );

      Residual[4*A  ] += gwts * ( NA * Res_Mas + GradNA_invF_ResMom );

      Residual[4*A+1] += gwts * ( NA * rho * detF * (vx_t - fx)
          + NA_x * P_iso(0) + NA_y * P_iso(1) + NA_z * P_iso(2)
          - GradNA_invF[0] * (detF * p - tau_c * Res_Mas) );

      Residual[4*A+2] += gwts * ( NA * rho * detF * (vy_t - fy)
          + NA_x * P_iso(3) + NA_y * P_iso(4) + NA_z * P_iso(5)
          - GradNA_invF[1] * (detF * p - tau_c * Res_Mas) );

      Residual[4*A+3] += gwts * ( NA * rho * detF * (vz_t - fz)
          + NA_x * P_iso(6) + NA_y * P_iso(7) + NA_z * P_iso(8)
          - GradNA_invF[2] * (detF * p - tau_c * Res_Mas) );
    }
  }
}


void PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &velo,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  double R[4], dR_dx[4], dR_dy[4], dR_dz[4];
  
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double h_e = element->get_h( eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  int ii, jj, qua, A, B, ii7;
  double p, p_t, p_x, p_y, p_z;
  double ux_t, uy_t, uz_t, vx, vy, vz;
  double vx_t, vy_t, vz_t;
  double ux_x, uy_x, uz_x, ux_y, uy_y, uz_y, ux_z, uy_z, uz_z;
  double vx_x, vy_x, vz_x, vx_y, vy_y, vz_y, vx_z, vy_z, vz_z;
  double fx, fy, fz, gwts, coor_x, coor_y, coor_z;
  double NA, NA_x, NA_y, NA_z;
  double NB, NB_x, NB_y, NB_z;

  double Res_Mom[3];
  double Res_Mas;
  double GradNA_invF_ResMom;
  double invFDV_t;

  double GradNB_invF_dot_Res_Mom;
  double GradNA_invF_dot_GradP_invF;
  double GradNA_invF_dot_GradNB_invF;
  double GradNA_invF_dot_part_Mom;

  Matrix_3x3 Dvelo_invF;

  const double curr = time + alpha_f * dt;

  const double dd_dv = alpha_f * gamma * dt;

  const double ddvm = dd_dv * dd_dv / alpha_m;

  Zero_Tangent_Residual();

  for(qua=0; qua < nqp; ++qua)
  {
    p = 0.0; p_t = 0.0; p_x = 0.0; p_y = 0.0; p_z = 0.0;

    ux_t = 0.0; uy_t = 0.0; uz_t = 0.0;
    vx   = 0.0; vy   = 0.0; vz   = 0.0;
    vx_t = 0.0; vy_t = 0.0; vz_t = 0.0;

    ux_x = 0.0; uy_x = 0.0; uz_x = 0.0;
    ux_y = 0.0; uy_y = 0.0; uz_y = 0.0;
    ux_z = 0.0; uy_z = 0.0; uz_z = 0.0;

    vx_x = 0.0; vy_x = 0.0; vz_x = 0.0;
    vx_y = 0.0; vy_y = 0.0; vz_y = 0.0;
    vx_z = 0.0; vy_z = 0.0; vz_z = 0.0;

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    element->get_R_gradR(qua, R, dR_dx, dR_dy, dR_dz);

    for(ii=0; ii<nLocBas; ++ii)
    {
      ii7 = 7*ii;
      p   += disp[ii7+3] * R[ii];
      p_x += disp[ii7+3] * dR_dx[ii];
      p_y += disp[ii7+3] * dR_dy[ii];
      p_z += disp[ii7+3] * dR_dz[ii];

      ux_t += velo[ii7  ] * R[ii];
      uy_t += velo[ii7+1] * R[ii];
      uz_t += velo[ii7+2] * R[ii];
      p_t  += velo[ii7+3] * R[ii];
      vx_t += velo[ii7+4] * R[ii];
      vy_t += velo[ii7+5] * R[ii];
      vz_t += velo[ii7+6] * R[ii];

      vx   += disp[ii7+4] * R[ii];
      vy   += disp[ii7+5] * R[ii];
      vz   += disp[ii7+6] * R[ii];

      ux_x += disp[ii7+0] * dR_dx[ii];
      uy_x += disp[ii7+1] * dR_dx[ii];
      uz_x += disp[ii7+2] * dR_dx[ii];

      ux_y += disp[ii7+0] * dR_dy[ii];
      uy_y += disp[ii7+1] * dR_dy[ii];
      uz_y += disp[ii7+2] * dR_dy[ii];

      ux_z += disp[ii7+0] * dR_dz[ii];
      uy_z += disp[ii7+1] * dR_dz[ii];
      uz_z += disp[ii7+2] * dR_dz[ii];

      vx_x += disp[ii7+4] * dR_dx[ii];
      vy_x += disp[ii7+5] * dR_dx[ii];
      vz_x += disp[ii7+6] * dR_dx[ii];

      vx_y += disp[ii7+4] * dR_dy[ii];
      vy_y += disp[ii7+5] * dR_dy[ii];
      vz_y += disp[ii7+6] * dR_dy[ii];

      vx_z += disp[ii7+4] * dR_dz[ii];
      vy_z += disp[ii7+5] * dR_dz[ii];
      vz_z += disp[ii7+6] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    gwts = element->get_detJac(qua) * quad->get_qw(qua);

    get_f(coor_x, coor_y, coor_z, curr, fx, fy, fz);

    const Matrix_3x3 F( ux_x + 1.0, ux_y, ux_z, uy_x, uy_y + 1.0, uy_z, uz_x, uz_y, uz_z + 1.0 );

    Matrix_3x3 invF( F ); invF.inverse();

    const Matrix_3x3 DVelo( vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z );

    Dvelo_invF.MatMult(DVelo, invF); // v_i,I invF_Ij = v_i,j

    double GradP_invF[3];
    invF.VecMultT( p_x, p_y, p_z, GradP_invF[0], GradP_invF[1], GradP_invF[2] ); // p_I invF_ii = p,i

    invFDV_t = invF.MatTContraction(DVelo); // invF_Ii V_i,I

    Matrix_3x3 P_iso, S_iso;
    Tensor4_3D AA_iso;
    matmodel->get_PK_FFStiffness(F, P_iso, S_iso, AA_iso);

    const double rho = matmodel->get_rho(p);
    const double drho = matmodel->get_drho_dp(p);

    const double mbeta = matmodel->get_beta(p);
    const double dmbeta = matmodel->get_dbeta_dp(p); 

    const double detF = F.det();

    // Get stabilization parameters
    double tau_m, tau_c;
    get_tau(tau_m, tau_c, dt, detF, h_e);

    // Residual of momentum equation
    Res_Mom[0] = rho * detF * vx_t;
    Res_Mom[1] = rho * detF * vy_t;
    Res_Mom[2] = rho * detF * vz_t;

    Res_Mom[0] += detF * ( invF(0) * p_x + invF(3) * p_y + invF(6) * p_z );
    Res_Mom[1] += detF * ( invF(1) * p_x + invF(4) * p_y + invF(7) * p_z );
    Res_Mom[2] += detF * ( invF(2) * p_x + invF(5) * p_y + invF(8) * p_z );

    Res_Mom[0] -= rho * detF * fx;
    Res_Mom[1] -= rho * detF * fy;
    Res_Mom[2] -= rho * detF * fz;

    // Residual of mass equation
    Res_Mas = detF * ( mbeta * p_t + invFDV_t );

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      // NA_I invF_Ii 
      double GradNA_invF[3];
      invF.VecMultT( NA_x, NA_y, NA_z, GradNA_invF[0], GradNA_invF[1], GradNA_invF[2] );

      // tau_m stabilization term
      GradNA_invF_ResMom = tau_m * ( GradNA_invF[0] * Res_Mom[0] + GradNA_invF[1] * Res_Mom[1] + GradNA_invF[2] * Res_Mom[2] );

      GradNA_invF_dot_GradP_invF = GradNA_invF[0] * GradP_invF[0] +
        GradNA_invF[1] * GradP_invF[1] + GradNA_invF[2] * GradP_invF[2];

      GradNA_invF_dot_part_Mom = GradNA_invF[0] * (vx_t - fx) 
        + GradNA_invF[1] * (vy_t - fy) + GradNA_invF[2] * (vz_t - fz);

      Residual[4*A  ] += gwts * ( NA * Res_Mas + GradNA_invF_ResMom );

      Residual[4*A+1] += gwts * ( NA * rho * detF * (vx_t - fx)
          + NA_x * P_iso(0) + NA_y * P_iso(1) + NA_z * P_iso(2)
          - GradNA_invF[0] * (detF * p - tau_c * Res_Mas) );

      Residual[4*A+2] += gwts * ( NA * rho * detF * (vy_t - fy)
          + NA_x * P_iso(3) + NA_y * P_iso(4) + NA_z * P_iso(5)
          - GradNA_invF[1] * (detF * p - tau_c * Res_Mas) );

      Residual[4*A+3] += gwts * ( NA * rho * detF * (vz_t - fz)
          + NA_x * P_iso(6) + NA_y * P_iso(7) + NA_z * P_iso(8)
          - GradNA_invF[2] * (detF * p - tau_c * Res_Mas) );

      for(B=0; B<nLocBas; ++B)
      {
        NB = R[B]; NB_x = dR_dx[B]; NB_y = dR_dy[B]; NB_z = dR_dz[B];

        const double NANBJ = detF * NA * NB;

        double GradNB_invF[3];
        invF.VecMultT( NB_x, NB_y, NB_z, GradNB_invF[0], GradNB_invF[1], GradNB_invF[2] );

        GradNB_invF_dot_Res_Mom = GradNB_invF[0] * Res_Mom[0]
          + GradNB_invF[1] * Res_Mom[1] + GradNB_invF[2] * Res_Mom[2];

        GradNA_invF_dot_GradNB_invF = GradNA_invF[0] * GradNB_invF[0] +
          GradNA_invF[1] * GradNB_invF[1] + GradNA_invF[2] * GradNB_invF[2];

        Tangent[ 16*nLocBas*A + 4*B ] += gwts * (alpha_m * mbeta * NANBJ
            + dd_dv * ( dmbeta * p_t * NANBJ 
              + GradNA_invF_dot_GradNB_invF * tau_m * detF 
              + tau_m * drho * detF * NB * GradNA_invF_dot_part_Mom ) );

        Tangent[ 16*nLocBas*A + 4*B + 1 ] += gwts * ( 
            alpha_m * GradNA_invF[0] * tau_m * rho * detF * NB
            + dd_dv * NA * detF * GradNB_invF[0]
            + ddvm * ( NA * detF * mbeta * p_t * GradNB_invF[0]
              + NA * detF * ( invFDV_t*GradNB_invF[0] 
                - Dvelo_invF(0) * GradNB_invF[0]
                - Dvelo_invF(3) * GradNB_invF[1] 
                - Dvelo_invF(6) * GradNB_invF[2] )
              - GradNA_invF[0] * tau_m * GradNB_invF_dot_Res_Mom
              + tau_m * detF * ( GradNA_invF_dot_GradP_invF * GradNB_invF[0]
                - GradNA_invF_dot_GradNB_invF * GradP_invF[0] )
              + tau_m * rho * detF * GradNA_invF_dot_part_Mom * GradNB_invF[0]
              ) );

        Tangent[ 16*nLocBas*A + 4*B + 2 ] += gwts * (
            alpha_m * GradNA_invF[1] * tau_m * rho * detF * NB
            + dd_dv * NA * detF * GradNB_invF[1]
            + ddvm * ( NA * detF * mbeta * p_t * GradNB_invF[1]
              + NA * detF * ( invFDV_t*GradNB_invF[1] 
                - Dvelo_invF(1) * GradNB_invF[0]
                - Dvelo_invF(4) * GradNB_invF[1] 
                - Dvelo_invF(7) * GradNB_invF[2] )
              - GradNA_invF[1] * tau_m * GradNB_invF_dot_Res_Mom
              + tau_m * detF * ( GradNA_invF_dot_GradP_invF * GradNB_invF[1]
                - GradNA_invF_dot_GradNB_invF * GradP_invF[1] )
              + tau_m * rho * detF * GradNA_invF_dot_part_Mom * GradNB_invF[1]
              ) );

        Tangent[ 16*nLocBas*A + 4*B + 3 ] += gwts * (
            alpha_m * GradNA_invF[2] * tau_m * rho * detF * NB
            + dd_dv * NA * detF * GradNB_invF[2]
            + ddvm * ( NA * detF * mbeta * p_t * GradNB_invF[2]
              + NA * detF * ( invFDV_t*GradNB_invF[2] 
                - Dvelo_invF(2) * GradNB_invF[0]
                - Dvelo_invF(5) * GradNB_invF[1] 
                - Dvelo_invF(8) * GradNB_invF[2] )
              - GradNA_invF[2] * tau_m * GradNB_invF_dot_Res_Mom
              + tau_m * detF * ( GradNA_invF_dot_GradP_invF * GradNB_invF[2]
                - GradNA_invF_dot_GradNB_invF * GradP_invF[2] )
              + tau_m * rho * detF * GradNA_invF_dot_part_Mom * GradNB_invF[2]
              ) );

        Tangent[4*nLocBas*(4*A+1)+4*B] += gwts * GradNA_invF[0] * detF * NB *
          (alpha_m * tau_c * mbeta - dd_dv*(1.0 - tau_c*dmbeta*p_t));

        Tangent[4*nLocBas*(4*A+2)+4*B] += gwts * GradNA_invF[1] * detF * NB *
          (alpha_m * tau_c * mbeta - dd_dv*(1.0 - tau_c*dmbeta*p_t));

        Tangent[4*nLocBas*(4*A+3)+4*B] += gwts * GradNA_invF[2] * detF * NB *
          (alpha_m * tau_c * mbeta - dd_dv*(1.0 - tau_c*dmbeta*p_t));

        const double mass_entry = gwts * NA * rho * detF * NB * alpha_m;

        Tangent[4*nLocBas*(4*A+1)+4*B+1] += mass_entry 
          + gwts * dd_dv * tau_c * detF * GradNA_invF[0] * GradNB_invF[0];

        Tangent[4*nLocBas*(4*A+1)+4*B+2] += gwts * dd_dv * tau_c * detF * GradNA_invF[0] * GradNB_invF[1];

        Tangent[4*nLocBas*(4*A+1)+4*B+3] += gwts * dd_dv * tau_c * detF * GradNA_invF[0] * GradNB_invF[2];

        Tangent[4*nLocBas*(4*A+2)+4*B+1] += gwts * dd_dv * tau_c * detF * GradNA_invF[1] * GradNB_invF[0];

        Tangent[4*nLocBas*(4*A+2)+4*B+2] += mass_entry
          + gwts * dd_dv * tau_c * detF * GradNA_invF[1] * GradNB_invF[1];

        Tangent[4*nLocBas*(4*A+2)+4*B+3] += gwts * dd_dv * tau_c * detF * GradNA_invF[1] * GradNB_invF[2];

        Tangent[4*nLocBas*(4*A+3)+4*B+1] += gwts * dd_dv * tau_c * detF * GradNA_invF[2] * GradNB_invF[0];

        Tangent[4*nLocBas*(4*A+3)+4*B+2] += gwts * dd_dv * tau_c * detF * GradNA_invF[2] * GradNB_invF[1];

        Tangent[4*nLocBas*(4*A+3)+4*B+3] += mass_entry
          + gwts * dd_dv * tau_c * detF * GradNA_invF[2] * GradNB_invF[2];

        const double geo_stiff = gwts * ddvm * (
            NA_x * ( S_iso(0) * NB_x + S_iso(1) * NB_y + S_iso(2) * NB_z )
            + NA_y * ( S_iso(3) * NB_x + S_iso(4) * NB_y + S_iso(5) * NB_z )
            + NA_z * ( S_iso(6) * NB_x + S_iso(7) * NB_y + S_iso(8) * NB_z) );

        Tangent[4*nLocBas*(4*A+1)+4*B+1] += geo_stiff;

        Tangent[4*nLocBas*(4*A+2)+4*B+2] += geo_stiff;

        Tangent[4*nLocBas*(4*A+3)+4*B+3] += geo_stiff;

        for(ii=1; ii<4; ++ii)
        {
          for(jj=1; jj<4; ++jj)
          {
            Tangent[ 4*nLocBas*(4*A + ii) + 4*B + jj ] += gwts * ddvm * (
                NA_x * (AA_iso(ii,0,jj,0) * NB_x + AA_iso(ii,0,jj,1) * NB_y
                  + AA_iso(ii,0,jj,2) * NB_z)
                + NA_y * (AA_iso(ii,1,jj,0) * NB_x + AA_iso(ii,1,jj,1) * NB_y
                  + AA_iso(ii,1,jj,2) * NB_z)
                + NA_z * (AA_iso(ii,2,jj,0) * NB_x + AA_iso(ii,2,jj,1) * NB_y
                  + AA_iso(ii,2,jj,2) * NB_z)
                - GradNA_invF[ii] * detF * p * GradNB_invF[jj]
                + GradNA_invF[jj] * detF * p * GradNB_invF[ii]
                - GradNA_invF[jj] * GradNB_invF[ii] * tau_c * Res_Mas
                + GradNA_invF[ii] * GradNB_invF[jj] * tau_c * mbeta * p_t * detF
                + GradNA_invF[ii] * GradNB_invF[jj] * tau_c * detF * invFDV_t
                - GradNA_invF[ii] * tau_c * detF * (Dvelo_invF(jj) * GradNB_invF[0] 
                  + Dvelo_invF(jj+3) * GradNB_invF[1] 
                  + Dvelo_invF(jj+6) * GradNB_invF[2]) );
          }
        }

        /* 
           Sub_Tan[5][index] += gwts * ddvm * NA * rho * detF * (vx_t - fx) * GradNB_invF[0];

           Sub_Tan[6][index] += gwts * ddvm * NA * rho * detF * (vx_t - fx) * GradNB_invF[1];

           Sub_Tan[7][index] += gwts * ddvm * NA * rho * detF * (vx_t - fx) * GradNB_invF[2];

           Sub_Tan[9][index] += gwts * ddvm * NA * rho * detF * (vy_t - fy) * GradNB_invF[0];

           Sub_Tan[10][index] += gwts * ddvm * NA * rho * detF * (vy_t - fy) * GradNB_invF[1];

           Sub_Tan[11][index] += gwts * ddvm * NA * rho * detF * (vy_t - fy) * GradNB_invF[2];

           Sub_Tan[13][index] += gwts * ddvm * NA * rho * detF * (vz_t - fz) * GradNB_invF[0];

           Sub_Tan[14][index] += gwts * ddvm * NA * rho * detF * (vz_t - fz) * GradNB_invF[1];

           Sub_Tan[15][index] += gwts * ddvm * NA * rho * detF * (vz_t - fz) * GradNB_invF[2];
           */
      } // Finish Loop-B
    } // Finish Loop-A
  } // Finish Loop-qua
}


void PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha::Assem_Mass_Residual(
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  double R[4], dR_dx[4], dR_dy[4], dR_dz[4];
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  int ii, qua, A, B, ii7;
  double p, ux_x, uy_x, uz_x, ux_y, uy_y, uz_y, ux_z, uy_z, uz_z;
  double vx_x, vy_x, vz_x, vx_y, vy_y, vz_y, vx_z, vy_z, vz_z;
  double vx, vy, vz;
  double fx, fy, fz, gwts, coor_x, coor_y, coor_z;
  double NA, NB, NA_x, NA_y, NA_z;

  double curr = 0.0;

  double invFDV_t;

  Zero_Tangent_Residual();

  for(qua=0; qua<nqp; ++qua)
  {
    p = 0.0; vx = 0.0; vy = 0.0; vz = 0.0;
    ux_x = 0.0; uy_x = 0.0; uz_x = 0.0;
    ux_y = 0.0; uy_y = 0.0; uz_y = 0.0;
    ux_z = 0.0; uy_z = 0.0; uz_z = 0.0;

    vx_x = 0.0; vy_x = 0.0; vz_x = 0.0;
    vx_y = 0.0; vy_y = 0.0; vz_y = 0.0;
    vx_z = 0.0; vy_z = 0.0; vz_z = 0.0;

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    element->get_R_gradR(qua, R, dR_dx, dR_dy, dR_dz);

    for(ii=0; ii<nLocBas; ++ii)
    {
      ii7 = 7 * ii;
      p += disp[ii7+3] * R[ii];

      vx += disp[ii7+4] * R[ii];
      vy += disp[ii7+5] * R[ii];
      vz += disp[ii7+6] * R[ii];

      ux_x += disp[ii7+0] * dR_dx[ii];
      uy_x += disp[ii7+1] * dR_dx[ii];
      uz_x += disp[ii7+2] * dR_dx[ii];

      ux_y += disp[ii7+0] * dR_dy[ii];
      uy_y += disp[ii7+1] * dR_dy[ii];
      uz_y += disp[ii7+2] * dR_dy[ii];

      ux_z += disp[ii7+0] * dR_dz[ii];
      uy_z += disp[ii7+1] * dR_dz[ii];
      uz_z += disp[ii7+2] * dR_dz[ii];

      vx_x += disp[ii7+4] * dR_dx[ii];
      vy_x += disp[ii7+5] * dR_dx[ii];
      vz_x += disp[ii7+6] * dR_dx[ii];

      vx_y += disp[ii7+4] * dR_dy[ii];
      vy_y += disp[ii7+5] * dR_dy[ii];
      vz_y += disp[ii7+6] * dR_dy[ii];

      vx_z += disp[ii7+4] * dR_dz[ii];
      vy_z += disp[ii7+5] * dR_dz[ii];
      vz_z += disp[ii7+6] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    gwts = element->get_detJac(qua) * quad->get_qw(qua);
    get_f(coor_x, coor_y, coor_z, curr, fx, fy, fz);

    const Matrix_3x3 F( ux_x + 1.0, ux_y, ux_z, uy_x, uy_y + 1.0, uy_z, uz_x, uz_y, uz_z + 1.0 );

    Matrix_3x3 invF( F ); invF.inverse();

    const Matrix_3x3 DVelo( vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z );

    // invF_Ii DV_i,I = v_i,i = div v
    invFDV_t = invF.MatTContraction(DVelo);

    Matrix_3x3 P_iso, S_iso;
    matmodel->get_PK(F, P_iso, S_iso);
    
    double mbeta = matmodel->get_beta(p);

    // use 1.0 in case of fully incompressible. 
    if( mbeta < 1.0e-5 ) mbeta = 1.0;

    const double rho = matmodel->get_rho(p);
    const double detF = F.det();

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      double GradNA_invF[3];
      invF.VecMultT( NA_x, NA_y, NA_z, GradNA_invF[0], GradNA_invF[1], GradNA_invF[2] );

      Residual[4*A  ] += gwts * NA * detF * invFDV_t;
      
      Residual[4*A+1] += gwts * ( NA_x * P_iso(0) + NA_y * P_iso(1) 
          + NA_z * P_iso(2) - GradNA_invF[0] * detF * p 
          - NA * rho * detF * fx );

      Residual[4*A+2] += gwts * ( NA_x * P_iso(3) + NA_y * P_iso(4) 
          + NA_z * P_iso(5) - GradNA_invF[1] * detF * p 
          - NA * rho * detF * fy );

      Residual[4*A+3] += gwts * ( NA_x * P_iso(6) + NA_y * P_iso(7) 
          + NA_z * P_iso(8) - GradNA_invF[2] * detF * p 
          - NA * rho * detF * fz );

      for(B=0; B<nLocBas; ++B)
      {
        NB = R[B];

        Tangent[4*nLocBas*(4*A)   + 4*B]   += gwts * NA * detF * mbeta * NB;
        Tangent[4*nLocBas*(4*A+1) + 4*B+1] += gwts * NA * rho * detF * NB;
        Tangent[4*nLocBas*(4*A+2) + 4*B+2] += gwts * NA * rho * detF * NB;
        Tangent[4*nLocBas*(4*A+3) + 4*B+3] += gwts * NA * rho * detF * NB;
      } // Finish loop-B
    } // Finish loop-A
  } // Finish loop-qua
}


void PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha::Assem_Residual_EBC(
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

  double gx, gy, gz, surface_area;
  const double curr = time + alpha_f * dt;

  Zero_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    get_ebc_fun( ebc_id, coor_x, coor_y, coor_z, curr,
        n_out.x(), n_out.y(), n_out.z(), gx, gy, gz );

    for(int A=0; A<snLocBas; ++A)
    {
      Residual[4*A+1] -= surface_area * quad -> get_qw(qua) * R[A] * gx;
      Residual[4*A+2] -= surface_area * quad -> get_qw(qua) * R[A] * gy;
      Residual[4*A+3] -= surface_area * quad -> get_qw(qua) * R[A] * gz;
    }
  }
}

// EOF
