#include "PLocAssem_Tet4_VMS_Seg_Incompressible.hpp"

PLocAssem_Tet4_VMS_Seg_Incompressible::PLocAssem_Tet4_VMS_Seg_Incompressible(
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


PLocAssem_Tet4_VMS_Seg_Incompressible::~PLocAssem_Tet4_VMS_Seg_Incompressible()
{
  delete [] Tangent; delete [] Residual; Tangent = nullptr; Residual = nullptr;
  if(num_ebc_fun > 0) delete [] flist;
}


void PLocAssem_Tet4_VMS_Seg_Incompressible::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("  Three-dimensional Hyper-elastic solid model with VMS mixed formulation, U-kinematic relation, Segregated, FEM formulation: \n");
  SYS_T::commPrint("\t  Spatial: Finite element with VMS stabilization \n");
  SYS_T::commPrint("\t  Temporal: Generalized-alpha method \n");
  SYS_T::commPrint("\t  Solid density rho0 = %e g/cm3\n\n", rho0);
  matmodel->print_info();
  SYS_T::commPrint("----------------------------------------------------------- \n");
}


void PLocAssem_Tet4_VMS_Seg_Incompressible::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 0.0;
}


void PLocAssem_Tet4_VMS_Seg_Incompressible::Zero_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
}


void PLocAssem_Tet4_VMS_Seg_Incompressible::Assem_Estimate()
{
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 1.0;
}


void PLocAssem_Tet4_VMS_Seg_Incompressible::get_tau( 
    double &tau_m_qua, double &tau_c_qua,
    const double &dt, const double &Jin,
    const double &dx ) const
{
  const double mu = matmodel->get_elastic_mu();
  //const double ka = matmodel->get_elastic_kappa();
  //const double c_max = std::pow( rho0 / (ka + 4*mu/3.0), -0.5);

  const double c_max = std::pow( rho0 / mu, -0.5); // Fully incompressible case

  //const double dt_ka = dx / c_max;
  const double dt_ka = 5.0e-3 / c_max;

  tau_m_qua = 1.0e-2 * dt_ka * Jin / rho0;
  tau_c_qua = 0.000 * dx * c_max * rho0 / Jin;
}


void PLocAssem_Tet4_VMS_Seg_Incompressible::Assem_Residual(
    const double &time, const double &dt,
    const double * const &velo,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &qua_prestress,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double h_e = element->get_h( eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double curr = time + alpha_f * dt;

  Zero_Residual();

  for(int qua=0; qua < nqp; ++qua)
  {
    double p = 0.0, p_t = 0.0, p_x = 0.0, p_y = 0.0, p_z = 0.0;

    double ux_t = 0.0, uy_t = 0.0, uz_t = 0.0;
    double vx   = 0.0, vy   = 0.0, vz   = 0.0;
    double vx_t = 0.0, vy_t = 0.0, vz_t = 0.0;

    double ux_x = 0.0, uy_x = 0.0, uz_x = 0.0;
    double ux_y = 0.0, uy_y = 0.0, uz_y = 0.0;
    double ux_z = 0.0, uy_z = 0.0, uz_z = 0.0;

    double vx_x = 0.0, vy_x = 0.0, vz_x = 0.0;
    double vx_y = 0.0, vy_y = 0.0, vz_y = 0.0;
    double vx_z = 0.0, vy_z = 0.0, vz_z = 0.0;

    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    double R[4], dR_dx[4], dR_dy[4], dR_dz[4];

    element->get_R_gradR( qua, R, dR_dx, dR_dy, dR_dz );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      p   += disp[ii*7+3] * R[ii];
      p_x += disp[ii*7+3] * dR_dx[ii];
      p_y += disp[ii*7+3] * dR_dy[ii];
      p_z += disp[ii*7+3] * dR_dz[ii];

      ux_t += velo[ii*7  ] * R[ii];
      uy_t += velo[ii*7+1] * R[ii];
      uz_t += velo[ii*7+2] * R[ii];
      p_t  += velo[ii*7+3] * R[ii];
      vx_t += velo[ii*7+4] * R[ii];
      vy_t += velo[ii*7+5] * R[ii];
      vz_t += velo[ii*7+6] * R[ii];

      vx   += disp[ii*7+4] * R[ii];
      vy   += disp[ii*7+5] * R[ii];
      vz   += disp[ii*7+6] * R[ii];

      ux_x += disp[ii*7+0] * dR_dx[ii];
      uy_x += disp[ii*7+1] * dR_dx[ii];
      uz_x += disp[ii*7+2] * dR_dx[ii];

      ux_y += disp[ii*7+0] * dR_dy[ii];
      uy_y += disp[ii*7+1] * dR_dy[ii];
      uz_y += disp[ii*7+2] * dR_dy[ii];

      ux_z += disp[ii*7+0] * dR_dz[ii];
      uy_z += disp[ii*7+1] * dR_dz[ii];
      uz_z += disp[ii*7+2] * dR_dz[ii];

      vx_x += disp[ii*7+4] * dR_dx[ii];
      vy_x += disp[ii*7+5] * dR_dx[ii];
      vz_x += disp[ii*7+6] * dR_dx[ii];

      vx_y += disp[ii*7+4] * dR_dy[ii];
      vy_y += disp[ii*7+5] * dR_dy[ii];
      vz_y += disp[ii*7+6] * dR_dy[ii];

      vx_z += disp[ii*7+4] * dR_dz[ii];
      vy_z += disp[ii*7+5] * dR_dz[ii];
      vz_z += disp[ii*7+6] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const Vector_3 f_body = get_f(coor_x, coor_y, coor_z, curr);

    const Matrix_3x3 F( ux_x + 1.0, ux_y, ux_z, uy_x, uy_y + 1.0, uy_z, uz_x, uz_y, uz_z + 1.0 );

    const Matrix_3x3 invF = inverse(F);

    const Matrix_3x3 DVelo( vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z );

    const double invFDV_t = invF.MatTContraction(DVelo); // invF_Ii V_i,I

    Matrix_3x3 P_iso, S_iso;
    matmodel->get_PK(F, P_iso, S_iso);

    // ------------------------------------------------------------------------
    // 1st PK stress corrected by prestress
    const Matrix_3x3 prestress( qua_prestress[qua*6+0], qua_prestress[qua*6+5], qua_prestress[qua*6+4],
        qua_prestress[qua*6+5], qua_prestress[qua*6+1], qua_prestress[qua*6+3],
        qua_prestress[qua*6+4], qua_prestress[qua*6+3], qua_prestress[qua*6+2] );

    P_iso += prestress * cofactor( F );
    // ------------------------------------------------------------------------

    const double rho = matmodel->get_rho(p);

    const double detF = F.det();

    // Get stabilization parameters
    double tau_m, tau_c;
    get_tau(tau_m, tau_c, dt, detF, h_e);

    // Residual of momentum equation
    const double Res_Mom[3] { rho * detF * ( vx_t - f_body.x() ) + detF * ( invF(0) * p_x + invF(3) * p_y + invF(6) * p_z ),
      rho * detF * ( vy_t - f_body.y() ) + detF * ( invF(1) * p_x + invF(4) * p_y + invF(7) * p_z ),
      rho * detF * ( vz_t - f_body.z() ) + detF * ( invF(2) * p_x + invF(5) * p_y + invF(8) * p_z ) };

    // Residual of mass equation
    const double Res_Mas = detF * invFDV_t;

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];
     
      // NA_I invF_Ii 
      const Vector_3 gradNA = invF.VecMultT( Vector_3(NA_x, NA_y, NA_z) );

      // tau_m stabilization term
      const double gradNA_ResMom = tau_m * ( gradNA.x() * Res_Mom[0] + gradNA.y() * Res_Mom[1] + gradNA.z() * Res_Mom[2] );

      Residual[4*A  ] += gwts * ( NA * Res_Mas + gradNA_ResMom );

      Residual[4*A+1] += gwts * ( NA * rho * detF * (vx_t - f_body.x())
          + NA_x * P_iso(0) + NA_y * P_iso(1) + NA_z * P_iso(2)
          - gradNA.x() * (detF * p - tau_c * Res_Mas) );

      Residual[4*A+2] += gwts * ( NA * rho * detF * (vy_t - f_body.y())
          + NA_x * P_iso(3) + NA_y * P_iso(4) + NA_z * P_iso(5)
          - gradNA.y() * (detF * p - tau_c * Res_Mas) );

      Residual[4*A+3] += gwts * ( NA * rho * detF * (vz_t - f_body.z())
          + NA_x * P_iso(6) + NA_y * P_iso(7) + NA_z * P_iso(8)
          - gradNA.z() * (detF * p - tau_c * Res_Mas) );
    }
  }
}


void PLocAssem_Tet4_VMS_Seg_Incompressible::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &velo,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &qua_prestress,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double h_e = element->get_h( eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double curr = time + alpha_f * dt;

  const double dd_dv = alpha_f * gamma * dt;
  
  const double ddvm = dd_dv * dd_dv / alpha_m;

  Zero_Tangent_Residual();

  for(int qua=0; qua < nqp; ++qua)
  {
    double p = 0.0, p_t = 0.0, p_x = 0.0, p_y = 0.0, p_z = 0.0;

    double ux_t = 0.0, uy_t = 0.0, uz_t = 0.0;
    double vx   = 0.0, vy   = 0.0, vz   = 0.0;
    double vx_t = 0.0, vy_t = 0.0, vz_t = 0.0;

    double ux_x = 0.0, uy_x = 0.0, uz_x = 0.0;
    double ux_y = 0.0, uy_y = 0.0, uz_y = 0.0;
    double ux_z = 0.0, uy_z = 0.0, uz_z = 0.0;

    double vx_x = 0.0, vy_x = 0.0, vz_x = 0.0;
    double vx_y = 0.0, vy_y = 0.0, vz_y = 0.0;
    double vx_z = 0.0, vy_z = 0.0, vz_z = 0.0;

    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    double R[4], dR_dx[4], dR_dy[4], dR_dz[4];
    
    element->get_R_gradR(qua, R, dR_dx, dR_dy, dR_dz);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      p   += disp[ii*7+3] * R[ii];
      p_x += disp[ii*7+3] * dR_dx[ii];
      p_y += disp[ii*7+3] * dR_dy[ii];
      p_z += disp[ii*7+3] * dR_dz[ii];

      ux_t += velo[ii*7  ] * R[ii];
      uy_t += velo[ii*7+1] * R[ii];
      uz_t += velo[ii*7+2] * R[ii];
      p_t  += velo[ii*7+3] * R[ii];
      vx_t += velo[ii*7+4] * R[ii];
      vy_t += velo[ii*7+5] * R[ii];
      vz_t += velo[ii*7+6] * R[ii];

      vx   += disp[ii*7+4] * R[ii];
      vy   += disp[ii*7+5] * R[ii];
      vz   += disp[ii*7+6] * R[ii];

      ux_x += disp[ii*7+0] * dR_dx[ii];
      uy_x += disp[ii*7+1] * dR_dx[ii];
      uz_x += disp[ii*7+2] * dR_dx[ii];

      ux_y += disp[ii*7+0] * dR_dy[ii];
      uy_y += disp[ii*7+1] * dR_dy[ii];
      uz_y += disp[ii*7+2] * dR_dy[ii];

      ux_z += disp[ii*7+0] * dR_dz[ii];
      uy_z += disp[ii*7+1] * dR_dz[ii];
      uz_z += disp[ii*7+2] * dR_dz[ii];

      vx_x += disp[ii*7+4] * dR_dx[ii];
      vy_x += disp[ii*7+5] * dR_dx[ii];
      vz_x += disp[ii*7+6] * dR_dx[ii];

      vx_y += disp[ii*7+4] * dR_dy[ii];
      vy_y += disp[ii*7+5] * dR_dy[ii];
      vz_y += disp[ii*7+6] * dR_dy[ii];

      vx_z += disp[ii*7+4] * dR_dz[ii];
      vy_z += disp[ii*7+5] * dR_dz[ii];
      vz_z += disp[ii*7+6] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const Vector_3 f_body = get_f(coor_x, coor_y, coor_z, curr);

    const Matrix_3x3 F( ux_x + 1.0, ux_y, ux_z, uy_x, uy_y + 1.0, uy_z, uz_x, uz_y, uz_z + 1.0 );

    const Matrix_3x3 invF = inverse(F);

    const Matrix_3x3 DVelo( vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z );

    const Matrix_3x3 Dvelo_invF = DVelo * invF;  // v_i,I invF_Ij = v_i,j

    double GradP_invF[3];
    invF.VecMultT( p_x, p_y, p_z, GradP_invF[0], GradP_invF[1], GradP_invF[2] ); // p_I invF_ii = p,i

    const double invFDV_t = invF.MatTContraction(DVelo); // invF_Ii V_i,I

    Matrix_3x3 P_iso, S_iso;
    Tensor4_3D AA_iso;
    matmodel->get_PK_FFStiffness(F, P_iso, S_iso, AA_iso);

    // ------------------------------------------------------------------------
    // 1st PK stress corrected by prestress
    const Matrix_3x3 prestress( qua_prestress[qua*6+0], qua_prestress[qua*6+5], qua_prestress[qua*6+4],
        qua_prestress[qua*6+5], qua_prestress[qua*6+1], qua_prestress[qua*6+3],
        qua_prestress[qua*6+4], qua_prestress[qua*6+3], qua_prestress[qua*6+2] );

    P_iso += prestress * cofactor( F );
    // ------------------------------------------------------------------------
    
    const double rho = matmodel->get_rho(p);
    const double detF = F.det();

    // Get stabilization parameters
    double tau_m, tau_c;
    get_tau(tau_m, tau_c, dt, detF, h_e);

    // Residual of momentum equation
    const double Res_Mom[3] { rho * detF * ( vx_t - f_body.x() ) + detF * ( invF(0) * p_x + invF(3) * p_y + invF(6) * p_z ),
      rho * detF * ( vy_t - f_body.y() ) + detF * ( invF(1) * p_x + invF(4) * p_y + invF(7) * p_z ),
      rho * detF * ( vz_t - f_body.z() ) + detF * ( invF(2) * p_x + invF(5) * p_y + invF(8) * p_z ) };
    
    // Residual of mass equation
    const double Res_Mas = detF * invFDV_t;

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];
      
      // NA_I invF_Ii 
      double GradNA_invF[3];
      invF.VecMultT( NA_x, NA_y, NA_z, GradNA_invF[0], GradNA_invF[1], GradNA_invF[2] );

      // tau_m stabilization term
      const double GradNA_invF_ResMom = tau_m * ( GradNA_invF[0] * Res_Mom[0] + GradNA_invF[1] * Res_Mom[1] + GradNA_invF[2] * Res_Mom[2] );

      const double GradNA_invF_dot_GradP_invF = GradNA_invF[0] * GradP_invF[0] +
        GradNA_invF[1] * GradP_invF[1] + GradNA_invF[2] * GradP_invF[2];

      const double GradNA_invF_dot_part_Mom = GradNA_invF[0] * (vx_t - f_body.x()) 
        + GradNA_invF[1] * (vy_t - f_body.y()) + GradNA_invF[2] * (vz_t - f_body.z());

      Residual[4*A  ] += gwts * ( NA * Res_Mas + GradNA_invF_ResMom );

      Residual[4*A+1] += gwts * ( NA * rho * detF * (vx_t - f_body.x())
          + NA_x * P_iso(0) + NA_y * P_iso(1) + NA_z * P_iso(2) 
          - GradNA_invF[0] * detF * p );

      Residual[4*A+2] += gwts * ( NA * rho * detF * (vy_t - f_body.y())
          + NA_x * P_iso(3) + NA_y * P_iso(4) + NA_z * P_iso(5) 
          - GradNA_invF[1] * detF * p );

      Residual[4*A+3] += gwts * ( NA * rho * detF * (vz_t - f_body.z())
          + NA_x * P_iso(6) + NA_y * P_iso(7) + NA_z * P_iso(8) 
          - GradNA_invF[2] * detF * p );

      for(int B=0; B<nLocBas; ++B)
      {
        const double NB = R[B], NB_x = dR_dx[B], NB_y = dR_dy[B], NB_z = dR_dz[B];
        
        double GradNB_invF[3];
        invF.VecMultT( NB_x, NB_y, NB_z, GradNB_invF[0], GradNB_invF[1], GradNB_invF[2] );

        const double GradNB_invF_dot_Res_Mom = GradNB_invF[0] * Res_Mom[0]
          + GradNB_invF[1] * Res_Mom[1] + GradNB_invF[2] * Res_Mom[2];

        const double GradNA_invF_dot_GradNB_invF = GradNA_invF[0] * GradNB_invF[0] +
          GradNA_invF[1] * GradNB_invF[1] + GradNA_invF[2] * GradNB_invF[2];

        Tangent[ 16*nLocBas*A + 4*B ] += gwts * dd_dv * GradNA_invF_dot_GradNB_invF * tau_m * detF;

        Tangent[ 16*nLocBas*A + 4*B + 1 ] += gwts * ( 
            alpha_m * GradNA_invF[0] * tau_m * detF * rho * NB
            + dd_dv * NA * detF * GradNB_invF[0]
            + ddvm * ( NA * detF * ( invFDV_t*GradNB_invF[0] 
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
            + ddvm * ( NA * detF * ( invFDV_t*GradNB_invF[1] 
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
            + ddvm * ( NA * detF * ( invFDV_t*GradNB_invF[2] 
                - Dvelo_invF(2) * GradNB_invF[0]
                - Dvelo_invF(5) * GradNB_invF[1] 
                - Dvelo_invF(8) * GradNB_invF[2] )
              - GradNA_invF[2] * tau_m * GradNB_invF_dot_Res_Mom
              + tau_m * detF * ( GradNA_invF_dot_GradP_invF * GradNB_invF[2]
                - GradNA_invF_dot_GradNB_invF * GradP_invF[2] ) 
              + tau_m * rho * detF * GradNA_invF_dot_part_Mom * GradNB_invF[2]
              ) );

        Tangent[4*nLocBas*(4*A+1)+4*B] -= gwts * GradNA_invF[0] * detF * NB * dd_dv;

        Tangent[4*nLocBas*(4*A+2)+4*B] -= gwts * GradNA_invF[1] * detF * NB * dd_dv;

        Tangent[4*nLocBas*(4*A+3)+4*B] -= gwts * GradNA_invF[2] * detF * NB * dd_dv;

        const double mass_entry = gwts * NA * NB * rho * detF * alpha_m;

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

        for(int ii=0; ii<3; ++ii)
        {
          for(int jj=0; jj<3; ++jj)
          {
            Tangent[ 4*nLocBas*(4*A + ii + 1) + 4*B + jj + 1 ] += gwts * ddvm * (
                NA_x * (AA_iso(ii,0,jj,0) * NB_x + AA_iso(ii,0,jj,1) * NB_y
                  + AA_iso(ii,0,jj,2) * NB_z)
                + NA_y * (AA_iso(ii,1,jj,0) * NB_x + AA_iso(ii,1,jj,1) * NB_y
                  + AA_iso(ii,1,jj,2) * NB_z)
                + NA_z * (AA_iso(ii,2,jj,0) * NB_x + AA_iso(ii,2,jj,1) * NB_y
                  + AA_iso(ii,2,jj,2) * NB_z)
                - GradNA_invF[ii] * detF * p * GradNB_invF[jj]
                + GradNA_invF[jj] * detF * p * GradNB_invF[ii]
                - GradNA_invF[jj] * GradNB_invF[ii] * tau_c * Res_Mas
                + GradNA_invF[ii] * GradNB_invF[jj] * tau_c * detF * invFDV_t
                - GradNA_invF[ii] * tau_c * detF * (Dvelo_invF(jj) * GradNB_invF[0] 
                  + Dvelo_invF(jj+3) * GradNB_invF[1] 
                  + Dvelo_invF(jj+6) * GradNB_invF[2]) );
          }
        }

      } // Finish Loop-B
    } // Finish Loop-A
  } // Finish Loop-qua
}


void PLocAssem_Tet4_VMS_Seg_Incompressible::Assem_Mass_Residual(
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &qua_prestress,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double curr = 0.0;

  Zero_Tangent_Residual();

  for(int qua=0; qua<nqp; ++qua)
  {
    double p = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
    double ux_x = 0.0, uy_x = 0.0, uz_x = 0.0;
    double ux_y = 0.0, uy_y = 0.0, uz_y = 0.0;
    double ux_z = 0.0, uy_z = 0.0, uz_z = 0.0;

    double vx_x = 0.0, vy_x = 0.0, vz_x = 0.0;
    double vx_y = 0.0, vy_y = 0.0, vz_y = 0.0;
    double vx_z = 0.0, vy_z = 0.0, vz_z = 0.0;

    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    double R[4], dR_dx[4], dR_dy[4], dR_dz[4];

    element->get_R_gradR(qua, R, dR_dx, dR_dy, dR_dz);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      p += disp[ii*7+3] * R[ii];

      vx += disp[ii*7+4] * R[ii];
      vy += disp[ii*7+5] * R[ii];
      vz += disp[ii*7+6] * R[ii];

      ux_x += disp[ii*7+0] * dR_dx[ii];
      uy_x += disp[ii*7+1] * dR_dx[ii];
      uz_x += disp[ii*7+2] * dR_dx[ii];

      ux_y += disp[ii*7+0] * dR_dy[ii];
      uy_y += disp[ii*7+1] * dR_dy[ii];
      uz_y += disp[ii*7+2] * dR_dy[ii];

      ux_z += disp[ii*7+0] * dR_dz[ii];
      uy_z += disp[ii*7+1] * dR_dz[ii];
      uz_z += disp[ii*7+2] * dR_dz[ii];

      vx_x += disp[ii*7+4] * dR_dx[ii];
      vy_x += disp[ii*7+5] * dR_dx[ii];
      vz_x += disp[ii*7+6] * dR_dx[ii];

      vx_y += disp[ii*7+4] * dR_dy[ii];
      vy_y += disp[ii*7+5] * dR_dy[ii];
      vz_y += disp[ii*7+6] * dR_dy[ii];

      vx_z += disp[ii*7+4] * dR_dz[ii];
      vy_z += disp[ii*7+5] * dR_dz[ii];
      vz_z += disp[ii*7+6] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);
    
    const Vector_3 f_body = get_f(coor_x, coor_y, coor_z, curr);

    const Matrix_3x3 F( ux_x + 1.0, ux_y, ux_z, uy_x, uy_y + 1.0, uy_z, uz_x, uz_y, uz_z + 1.0 );

    const Matrix_3x3 invF = inverse(F);

    const Matrix_3x3 DVelo(  vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z );

    // invF_Ii DV_i,I = v_i,i = div v
    const double invFDV_t = invF.MatTContraction(DVelo);

    Matrix_3x3 P_iso, S_iso;
    matmodel->get_PK(F, P_iso, S_iso);

    // ------------------------------------------------------------------------
    // 1st PK stress corrected by prestress
    const Matrix_3x3 prestress( qua_prestress[qua*6+0], qua_prestress[qua*6+5], qua_prestress[qua*6+4],
        qua_prestress[qua*6+5], qua_prestress[qua*6+1], qua_prestress[qua*6+3],
        qua_prestress[qua*6+4], qua_prestress[qua*6+3], qua_prestress[qua*6+2] );

    P_iso += prestress * cofactor( F );
    // ------------------------------------------------------------------------
    
    double mbeta = matmodel->get_beta(p);

    // use 1.0 in case of fully incompressible. 
    if( mbeta <1.0e-5 ) mbeta = 1.0;

    const double rho = matmodel->get_rho(p);
    const double detF = F.det();

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      const Vector_3 gradNA = invF.VecMultT( Vector_3(NA_x, NA_y, NA_z) );

      Residual[4*A  ] += gwts * NA * detF * invFDV_t;

      Residual[4*A+1] += gwts * ( NA_x * P_iso(0) + NA_y * P_iso(1) 
          + NA_z * P_iso(2) - gradNA.x() * detF * p 
          - NA * rho * detF * f_body.x() );

      Residual[4*A+2] += gwts * ( NA_x * P_iso(3) + NA_y * P_iso(4) 
          + NA_z * P_iso(5) - gradNA.y() * detF * p 
          - NA * rho * detF * f_body.y() );

      Residual[4*A+3] += gwts * ( NA_x * P_iso(6) + NA_y * P_iso(7) 
          + NA_z * P_iso(8) - gradNA.z() * detF * p 
          - NA * rho * detF * f_body.z() );

      for(int B=0; B<nLocBas; ++B)
      {
        Tangent[4*nLocBas*(4*A)   + 4*B]   += gwts * NA * detF * mbeta * R[B];
        Tangent[4*nLocBas*(4*A+1) + 4*B+1] += gwts * NA * detF * rho * R[B];
        Tangent[4*nLocBas*(4*A+2) + 4*B+2] += gwts * NA * detF * rho * R[B];
        Tangent[4*nLocBas*(4*A+3) + 4*B+3] += gwts * NA * detF * rho * R[B];
      } // Finish loop-B
    } // Finish loop-A
  } // Finish loop-qua
}


void PLocAssem_Tet4_VMS_Seg_Incompressible::Assem_Residual_EBC(
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

  const double curr = time + alpha_f * dt;

  Zero_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    double surface_area;
    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    const Vector_3 gg = get_ebc_fun( ebc_id, coor_x, coor_y, coor_z, curr,
        n_out.x(), n_out.y(), n_out.z() );

    for(int A=0; A<snLocBas; ++A)
    {
      Residual[4*A+1] -= surface_area * quad -> get_qw(qua) * R[A] * gg.x();
      Residual[4*A+2] -= surface_area * quad -> get_qw(qua) * R[A] * gg.y();
      Residual[4*A+3] -= surface_area * quad -> get_qw(qua) * R[A] * gg.z();
    }
  }
}


std::vector<Matrix_3x3> PLocAssem_Tet4_VMS_Seg_Incompressible::get_Wall_CauchyStress(
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad ) const
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  std::vector<Matrix_3x3> stress( nqp );

  for( int qua = 0; qua < nqp; ++qua )
  {
    double dR_dx[4], dR_dy[4], dR_dz[4];

    element->get_gradR( qua, dR_dx, dR_dy, dR_dz );

    double ux_x = 0.0, uy_x = 0.0, uz_x = 0.0;
    double ux_y = 0.0, uy_y = 0.0, uz_y = 0.0;
    double ux_z = 0.0, uy_z = 0.0, uz_z = 0.0;

    for(int ii=0; ii<nLocBas; ++ii)
    {
      ux_x += disp[ii*7+0] * dR_dx[ii];
      uy_x += disp[ii*7+1] * dR_dx[ii];
      uz_x += disp[ii*7+2] * dR_dx[ii];

      ux_y += disp[ii*7+0] * dR_dy[ii];
      uy_y += disp[ii*7+1] * dR_dy[ii];
      uz_y += disp[ii*7+2] * dR_dy[ii];

      ux_z += disp[ii*7+0] * dR_dz[ii];
      uy_z += disp[ii*7+1] * dR_dz[ii];
      uz_z += disp[ii*7+2] * dR_dz[ii];
    }

    const Matrix_3x3 F( ux_x + 1.0, ux_y, ux_z, uy_x, uy_y + 1.0, uy_z, uz_x, uz_y, uz_z + 1.0 );

    stress[qua] = matmodel -> get_Cauchy_stress( F );
  }

  return stress;
}


void PLocAssem_Tet4_VMS_Seg_Incompressible::Assem_Residual_EBC(
    const double &time,
    const double * const &vec,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  //const double factor = time >= 1.0 ? 1.0 : time;
  const double factor = 0.0;

  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  Zero_Residual();

  for(int qua = 0; qua < quad -> get_num_quadPts(); ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    double surface_area;
    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    double pp = 0.0;
    for(int ii=0; ii<snLocBas; ++ii) pp += vec[ii*7+3] * R[ii];

    for(int A=0; A<snLocBas; ++A)
    {
      Residual[4*A+1] -= surface_area * quad -> get_qw(qua) * R[A] * (-1.0) * factor * pp * n_out.x();
      Residual[4*A+2] -= surface_area * quad -> get_qw(qua) * R[A] * (-1.0) * factor * pp * n_out.y();
      Residual[4*A+3] -= surface_area * quad -> get_qw(qua) * R[A] * (-1.0) * factor * pp * n_out.z();
    }
  }
}

// EOF
