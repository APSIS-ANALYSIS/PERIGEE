#include "PLocAssem_LPS0_Seg_GenAlpha.hpp"

PLocAssem_LPS0_Seg_GenAlpha::PLocAssem_LPS0_Seg_GenAlpha(
    IMaterialModel * const &in_matmodel,
    const TimeMethod_GenAlpha * const &tm_gAlpha,
    const int &in_nlocbas, const int &in_nqp,
    const int &in_snlocbas )
: alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), num_ebc_fun(1),
  nLocBas(in_nlocbas), dof_per_node(7), vec_size(nLocBas * 4),
  nqp(in_nqp), snLocBas(in_snlocbas)
{
  matmodel = in_matmodel;

  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  Zero_Tangent_Residual();

  R = new double [nLocBas];
  dR_dx = new double [nLocBas];
  dR_dy = new double [nLocBas];
  dR_dz = new double [nLocBas];
  
  Sub_Tan = new double * [16];
  for(int ii=0; ii<16; ++ii) Sub_Tan[ii] = new double [nLocBas * nLocBas];

  if( num_ebc_fun == 0 ) flist = NULL;
  else flist = new locassem_lps0_funs [num_ebc_fun];

  flist[0] = &PLocAssem_LPS0_Seg_GenAlpha::get_rig_h;

  print_info();
}


PLocAssem_LPS0_Seg_GenAlpha::~PLocAssem_LPS0_Seg_GenAlpha()
{
  delete [] Tangent; Tangent = NULL;
  delete [] Residual; Residual = NULL;
  delete [] R; R = NULL;
  delete [] dR_dx; dR_dx = NULL;
  delete [] dR_dy; dR_dy = NULL;
  delete [] dR_dz; dR_dz = NULL;
  
  for(int ii=0; ii<16; ++ii) delete [] Sub_Tan[ii];

  delete [] Sub_Tan;
  if(num_ebc_fun > 0) delete [] flist;
}


void PLocAssem_LPS0_Seg_GenAlpha::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Three-dimensional Hyper-elastic solid model with LPS0 mixed formulation, Segregated, FEM formulation: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Spatial: Finite element with LPS stabilization to zero-th order polynomial \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Temporal: Generalized-alpha method \n");
  matmodel->print_info();
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PLocAssem_LPS0_Seg_GenAlpha::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 0.0;
}


void PLocAssem_LPS0_Seg_GenAlpha::Zero_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
}


void PLocAssem_LPS0_Seg_GenAlpha::Assem_Estimate()
{
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 1.0;
}


void PLocAssem_LPS0_Seg_GenAlpha::Assem_Residual(
    const double &time, const double &dt,
    const double * const &velo,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  // build the basis funciton based on the current configuration
  element -> buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  // Clean the residual array
  Zero_Residual();
 
  // get element size 
  const double h_e = element->get_h( eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );
  
  // get stabilization parameter from h
  const double tau_a = get_tau(dt, h_e);

  int ii, qua, A, ii7;
  double p, p_t, p_x, p_y, p_z;
  double ux_x, uy_x, uz_x, ux_y, uy_y, uz_y, ux_z, uy_z, uz_z;
  double vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z, vx_t, vy_t, vz_t;
  double fx, fy, fz, gwts, coor_x, coor_y, coor_z;
  double NA, NA_x, NA_y, NA_z;

  const double curr = time + alpha_f * dt;

  // Aux variables
  double ResMas, div_v;
  double gradNA[3];
  double gradP[3];
  Matrix_3x3 DVelo;

  for(qua=0; qua < nqp; ++qua)
  {
    p = 0.0; p_t = 0.0; p_x = 0.0; p_y = 0.0; p_z = 0.0;
    
    ux_x = 0.0; uy_x = 0.0; uz_x = 0.0;
    ux_y = 0.0; uy_y = 0.0; uz_y = 0.0;
    ux_z = 0.0; uy_z = 0.0; uz_z = 0.0;
    
    vx_x = 0.0; vx_y = 0.0; vx_z = 0.0;
    vy_x = 0.0; vy_y = 0.0; vy_z = 0.0;
    vz_x = 0.0; vz_y = 0.0; vz_z = 0.0;
    
    vx_t = 0.0; vy_t = 0.0; vz_t = 0.0;
    
    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    element -> get_R_gradR( qua, R, dR_dx, dR_dy, dR_dz );
    
    for(ii=0; ii<nLocBas; ++ii)
    {
      ii7 = 7 * ii;
      p   += disp[ii7+3] * R[ii];
      p_x += disp[ii7+3] * dR_dx[ii];
      p_y += disp[ii7+3] * dR_dy[ii];
      p_z += disp[ii7+3] * dR_dz[ii];
      
      p_t  += velo[ii7+3] * R[ii];
      vx_t += velo[ii7+4] * R[ii];
      vy_t += velo[ii7+5] * R[ii];
      vz_t += velo[ii7+6] * R[ii];

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
      vx_y += disp[ii7+4] * dR_dy[ii];
      vx_z += disp[ii7+4] * dR_dz[ii];

      vy_x += disp[ii7+5] * dR_dx[ii];
      vy_y += disp[ii7+5] * dR_dy[ii];
      vy_z += disp[ii7+5] * dR_dz[ii];

      vz_x += disp[ii7+6] * dR_dx[ii];
      vz_y += disp[ii7+6] * dR_dy[ii];
      vz_z += disp[ii7+6] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }
    
    gwts = element -> get_detJac(qua) * quad -> get_qw(qua);

    get_f(coor_x, coor_y, coor_z, curr, fx, fy, fz);
    
    F(0) = ux_x + 1.0; F(1) = ux_y;       F(2) = ux_z;
    F(3) = uy_x;       F(4) = uy_y + 1.0; F(5) = uy_z;
    F(6) = uz_x;       F(7) = uz_y;       F(8) = uz_z + 1.0;

    detF = F.det();

    invF.copy(F);
    invF.inverse();
    
    matmodel->get_PK(F, P_iso, S_iso);

    rho = matmodel -> get_rho(p);
    mbeta = matmodel -> get_beta(p);

    DVelo(0) = vx_x; DVelo(1) = vx_y; DVelo(2) = vx_z;
    DVelo(3) = vy_x; DVelo(4) = vy_y; DVelo(5) = vy_z;
    DVelo(6) = vz_x; DVelo(7) = vz_y; DVelo(8) = vz_z;

    div_v = invF.MatTContraction(DVelo); // invF_Ii V_i,I
    
    ResMas = detF * ( mbeta * p_t + div_v );

    invF.VecMultT( p_x, p_y, p_z, gradP ); // gradp_i = p_I invF_Ii

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];
      
      invF.VecMultT( NA_x, NA_y, NA_z, gradNA );

      Residual[4*A  ] += gwts * ( NA * ResMas + tau_a * detF * 
          (gradNA[0]*gradP[0] + gradNA[1]*gradP[1] + gradNA[2]*gradP[2]) );

      Residual[4*A+1] += gwts * ( NA * detF * rho * (vx_t - fx) 
          + NA_x * P_iso(0) + NA_y * P_iso(1) + NA_z * P_iso(2)
          - detF * gradNA[0] * p );

      Residual[4*A+2] += gwts * ( NA * detF * rho * (vy_t - fy) 
          + NA_x * P_iso(3) + NA_y * P_iso(4) + NA_z * P_iso(5)
          - detF * gradNA[1] * p );

      Residual[4*A+3] += gwts * ( NA * detF * rho * (vz_t - fz) 
          + NA_x * P_iso(6) + NA_y * P_iso(7) + NA_z * P_iso(8)
          - detF * gradNA[2] * p );
    }
  } // Finish-loop-over-qua
}


void PLocAssem_LPS0_Seg_GenAlpha::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &velo,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  // build the basis funciton based on the current configuration
  element -> buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  // Clean the residual array
  Zero_Tangent_Residual();
  Zero_Sub_Tan();
 
  // get element size 
  const double h_e = element->get_h( eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );
  
  // get stabilization parameter from h
  const double tau_a = get_tau(dt, h_e);

  int ii, jj, qua, A, B, index, ii7;
  double p, p_t, p_x, p_y, p_z;
  double ux_x, uy_x, uz_x, ux_y, uy_y, uz_y, ux_z, uy_z, uz_z;
  double vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z; 
  double vx_t, vy_t, vz_t;
  double fx, fy, fz, gwts, coor_x, coor_y, coor_z;
  double NA, NA_x, NA_y, NA_z, NB, NB_x, NB_y, NB_z;
  
  double gradNA[3];
  double gradNB[3];
  double gradP[3];
  double gradvx[3];
  double gradvy[3];
  double gradvz[3];
  double gradNA_dot_gradNB; // NA_xi NB_xi
  double gradNA_dot_gradP;  // NA_xi p_xi
  double gradNB_dot_gradP;  // NB_xi p_xi
  double div_v, ResMas, Res_mom_x, Res_mom_y, Res_mom_z;
  Matrix_3x3 DVelo;

  const double curr = time + alpha_f * dt;
  const double dd_dv = alpha_f * gamma * dt;
  const double ddvm = dd_dv * dd_dv / alpha_m;

  for(qua=0; qua < nqp; ++qua)
  {
    p = 0.0; p_t = 0.0; p_x = 0.0; p_y = 0.0; p_z = 0.0;
    
    ux_x = 0.0; uy_x = 0.0; uz_x = 0.0;
    ux_y = 0.0; uy_y = 0.0; uz_y = 0.0;
    ux_z = 0.0; uy_z = 0.0; uz_z = 0.0;
    vx_x = 0.0; vx_y = 0.0; vx_z = 0.0;
    vy_x = 0.0; vy_y = 0.0; vy_z = 0.0;
    vz_x = 0.0; vz_y = 0.0; vz_z = 0.0;
    vx_t = 0.0; vy_t = 0.0; vz_t = 0.0;
    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    element -> get_R_gradR( qua, R, dR_dx, dR_dy, dR_dz );
    
    for(ii=0; ii<nLocBas; ++ii)
    {
      ii7 = 7 * ii;
      p   += disp[ii7+3] * R[ii];
      p_x += disp[ii7+3] * dR_dx[ii];
      p_y += disp[ii7+3] * dR_dy[ii];
      p_z += disp[ii7+3] * dR_dz[ii];
    
      p_t  += velo[ii7+3] * R[ii];
      vx_t += velo[ii7+4] * R[ii];
      vy_t += velo[ii7+5] * R[ii];
      vz_t += velo[ii7+6] * R[ii];

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
      vx_y += disp[ii7+4] * dR_dy[ii];
      vx_z += disp[ii7+4] * dR_dz[ii];
      
      vy_x += disp[ii7+5] * dR_dx[ii];
      vy_y += disp[ii7+5] * dR_dy[ii];
      vy_z += disp[ii7+5] * dR_dz[ii];
      
      vz_x += disp[ii7+6] * dR_dx[ii];
      vz_y += disp[ii7+6] * dR_dy[ii];
      vz_z += disp[ii7+6] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }
    gwts = element -> get_detJac(qua) * quad -> get_qw(qua);

    get_f(coor_x, coor_y, coor_z, curr, fx, fy, fz);
    
    F(0) = ux_x + 1.0; F(1) = ux_y;       F(2) = ux_z;
    F(3) = uy_x;       F(4) = uy_y + 1.0; F(5) = uy_z;
    F(6) = uz_x;       F(7) = uz_y;       F(8) = uz_z + 1.0;
    
    detF = F.det();

    invF.copy(F);
    invF.inverse();

    matmodel->get_PK_FFStiffness(F, P_iso, S_iso, AA_iso);
    rho = matmodel -> get_rho(p);
    mbeta = matmodel -> get_beta(p);
    drho = matmodel -> get_drho_dp(p);
    dmbeta = matmodel -> get_dbeta_dp(p);
    
    DVelo(0) = vx_x; DVelo(1) = vx_y; DVelo(2) = vx_z;
    DVelo(3) = vy_x; DVelo(4) = vy_y; DVelo(5) = vy_z;
    DVelo(6) = vz_x; DVelo(7) = vz_y; DVelo(8) = vz_z;

    div_v = invF.MatTContraction(DVelo);

    ResMas = detF * ( mbeta * p_t + div_v );
    
    invF.VecMultT( p_x, p_y, p_z, gradP );
    invF.VecMultT( vx_x, vx_y, vx_z, gradvx );
    invF.VecMultT( vy_x, vy_y, vy_z, gradvy );
    invF.VecMultT( vz_x, vz_y, vz_z, gradvz );

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      invF.VecMultT( NA_x, NA_y, NA_z, gradNA );

      gradNA_dot_gradP = gradNA[0]*gradP[0]+gradNA[1]*gradP[1]+gradNA[2]*gradP[2];

      Res_mom_x = NA * detF * rho * (vx_t - fx) - detF * gradNA[0] * p;
      Res_mom_y = NA * detF * rho * (vy_t - fy) - detF * gradNA[1] * p;
      Res_mom_z = NA * detF * rho * (vz_t - fz) - detF * gradNA[2] * p;

      Residual[4*A  ] += gwts * ( NA * ResMas + tau_a * detF * gradNA_dot_gradP );

      Residual[4*A+1] += gwts * ( Res_mom_x 
          + NA_x * P_iso(0) + NA_y * P_iso(1) + NA_z * P_iso(2) );

      Residual[4*A+2] += gwts * ( Res_mom_y 
          + NA_x * P_iso(3) + NA_y * P_iso(4) + NA_z * P_iso(5) );

      Residual[4*A+3] += gwts * ( Res_mom_z
          + NA_x * P_iso(6) + NA_y * P_iso(7) + NA_z * P_iso(8) );

      for(B=0; B<nLocBas; ++B)
      {
        index = A * nLocBas + B;
        NB = R[B]; NB_x = dR_dx[B]; NB_y = dR_dy[B]; NB_z = dR_dz[B];
        
        invF.VecMultT( NB_x, NB_y, NB_z, gradNB );

        gradNA_dot_gradNB = gradNA[0]*gradNB[0] 
          + gradNA[1]*gradNB[1] + gradNA[2]*gradNB[2];

        gradNB_dot_gradP = gradNB[0] * gradP[0]
          + gradNB[1]*gradP[1] + gradNB[2]*gradP[2];

        Sub_Tan[0][index] += gwts * detF * ( alpha_m * mbeta * NA * NB
            + dd_dv * ( dmbeta * p_t * NA * NB + tau_a * gradNA_dot_gradNB ) );

        Sub_Tan[1][index] += gwts * detF * ( dd_dv * NA * gradNB[0]
            + ddvm * ( NA*gradNB[0]*(mbeta*p_t + div_v)
              - NA*(gradvx[0]*gradNB[0]+gradvx[1]*gradNB[1]+gradvx[2]*gradNB[2])
              + tau_a * gradNA_dot_gradP * gradNB[0] 
              - tau_a * gradNA_dot_gradNB * gradP[0]
              - tau_a * gradNA[0] * gradNB_dot_gradP ) );

        Sub_Tan[2][index] += gwts * ( dd_dv *  NA * gradNB[1] 
            + ddvm * ( NA*gradNB[1]*(mbeta*p_t + div_v) 
              - NA*(gradvy[0]*gradNB[0]+gradvy[1]*gradNB[1]+gradvy[2]*gradNB[2])
              + tau_a * gradNA_dot_gradP * gradNB[1] 
              - tau_a * gradNA_dot_gradNB * gradP[1]
              - tau_a * gradNA[1] * gradNB_dot_gradP ) );

        Sub_Tan[3][index] += gwts * ( dd_dv * NA * NB_z
            + ddvm * ( NA*gradNB[2]*(mbeta*p_t + div_v) 
              - NA*(gradvz[0]*gradNB[0]+gradvz[1]*gradNB[1]+gradvz[2]*gradNB[2])
              + tau_a * gradNA_dot_gradP * gradNB[2]
              - tau_a * gradNA_dot_gradNB * gradP[2]
              - tau_a * gradNA[2] * gradNB_dot_gradP ) );

        Sub_Tan[4][index]  -= gwts * dd_dv * detF * ( gradNA[0] * NB
            - NA * NB * drho * (vx_t - fx) );

        Sub_Tan[8][index]  -= gwts * dd_dv * detF * ( gradNA[1] * NB
            - NA * NB * drho * (vy_t - fy) );

        Sub_Tan[12][index] -= gwts * dd_dv * detF * ( gradNA[2] * NB
            - NA * NB * drho * (vz_t - fz) );

        const double mass_entry = gwts * alpha_m * NA * NB * detF * rho;
        const double geo_stiff = gwts * ddvm * (
            NA_x * ( S_iso(0) * NB_x + S_iso(1) * NB_y + S_iso(2) * NB_z )
            + NA_y * ( S_iso(3) * NB_x + S_iso(4) * NB_y + S_iso(5) * NB_z )
            + NA_z * ( S_iso(6) * NB_x + S_iso(7) * NB_y + S_iso(8) * NB_z) );

        Sub_Tan[5][index] += mass_entry + geo_stiff
          + gwts * ddvm * ( gradNB[0]*Res_mom_x + gradNA[0]*gradNB[0]*detF*p );

        Sub_Tan[6][index] += gwts * ddvm * (
            gradNB[1] * Res_mom_x + gradNA[1] * gradNB[0] * detF * p);

        Sub_Tan[7][index] += gwts * ddvm * (
            gradNB[2] * Res_mom_x + gradNA[2] * gradNB[0] * detF * p);

        Sub_Tan[9][index] += gwts * ddvm * (
            gradNB[0] * Res_mom_y + gradNA[0] * gradNB[1] * detF * p);

        Sub_Tan[10][index] += mass_entry + geo_stiff
          + gwts * ddvm * (gradNB[1] * Res_mom_y + gradNA[1]*gradNB[1]*detF*p);

        Sub_Tan[11][index] += gwts * ddvm * (
            gradNB[2] * Res_mom_y + gradNA[2] * gradNB[1] * detF * p);

        Sub_Tan[13][index] += gwts * ddvm * (
            gradNB[0] * Res_mom_z + gradNA[0] * gradNB[2] * detF * p);

        Sub_Tan[14][index] += gwts * ddvm * (
            gradNB[1] * Res_mom_z + gradNA[1] * gradNB[2] * detF * p);

        Sub_Tan[15][index] += mass_entry + geo_stiff
          + gwts * ddvm * (gradNB[2] * Res_mom_z + gradNA[2]*gradNB[2]*detF*p);

        // Material stiffness
        for(ii=0; ii<3; ++ii)
        {
          for(jj=0; jj<3; ++jj)
          {
            Sub_Tan[4*ii+jj+5][index] += gwts * ddvm * (
                NA_x * ( AA_iso(ii,0,jj,0)*NB_x + AA_iso(ii,0,jj,1)*NB_y
                  + AA_iso(ii,0,jj,2)*NB_z )
                + NA_y * ( AA_iso(ii,1,jj,0)*NB_x + AA_iso(ii,1,jj,1)*NB_y
                  + AA_iso(ii,1,jj,2) * NB_z )
                + NA_z * (AA_iso(ii,2,jj,0)*NB_x + AA_iso(ii,2,jj,1)*NB_y
                  + AA_iso(ii,2,jj,2) * NB_z ) );
          }
        }
      } // Finish loop-B
    } // Finish loop-A
  } // Finish-loop-over-qua

  for(ii=0; ii<4; ++ii)
  {
    for(jj=0; jj<4; ++jj)
    {
      for(A=0; A<nLocBas; ++A)
      {
        for(B=0; B<nLocBas; ++B)
        {
          Tangent[ 4*nLocBas*(4*A + ii) + 4*B + jj ]
            = Sub_Tan[ii*4+jj][A*nLocBas + B];
        }
      }
    }
  }
}


void PLocAssem_LPS0_Seg_GenAlpha::Assem_Mass_Residual(
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element -> buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  Zero_Tangent_Residual();
  Zero_Sub_Tan();

  int ii, jj, qua, A, B, index, ii7;
  double p, p_x, p_y, p_z;
  double ux_x, uy_x, uz_x, ux_y, uy_y, uz_y, ux_z, uy_z, uz_z;
  double vx_x, vy_x, vz_x, vx_y, vy_y, vz_y, vx_z, vy_z, vz_z;
  double fx, fy, fz, gwts, coor_x, coor_y, coor_z;
  double NA, NA_x, NA_y, NA_z, NB;

  const double curr = 0.0;

  // Aux variables
  double gradNA[3];
  double gradP[3];
  double div_v;
  Matrix_3x3 DVelo;

  for(qua=0; qua < nqp; ++qua)
  {
    p = 0.0; p_x = 0.0; p_y = 0.0; p_z = 0.0;

    ux_x = 0.0; uy_x = 0.0; uz_x = 0.0;
    ux_y = 0.0; uy_y = 0.0; uz_y = 0.0;
    ux_z = 0.0; uy_z = 0.0; uz_z = 0.0; 
    
    vx_x = 0.0; vy_x = 0.0; vz_x = 0.0;
    vx_y = 0.0; vy_y = 0.0; vz_y = 0.0;
    vx_z = 0.0; vy_z = 0.0; vz_z = 0.0;

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    element -> get_R_gradR( qua, R, dR_dx, dR_dy, dR_dz );

    for(ii=0; ii<nLocBas; ++ii)
    {
      ii7 = 7 * ii;
      p   += disp[ii7+3] * R[ii];
      p_x += disp[ii7+3] * dR_dx[ii];
      p_y += disp[ii7+3] * dR_dy[ii];
      p_z += disp[ii7+3] * dR_dz[ii];

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
    gwts = element -> get_detJac(qua) * quad -> get_qw(qua);

    get_f(coor_x, coor_y, coor_z, curr, fx, fy, fz);

    F(0) = ux_x + 1.0; F(1) = ux_y;       F(2) = ux_z;
    F(3) = uy_x;       F(4) = uy_y + 1.0; F(5) = uy_z;
    F(6) = uz_x;       F(7) = uz_y;       F(8) = uz_z + 1.0;

    detF = F.det();

    invF.copy(F);
    invF.inverse();

    DVelo(0) = vx_x; DVelo(1) = vx_y; DVelo(2) = vx_z;
    DVelo(3) = vy_x; DVelo(4) = vy_y; DVelo(5) = vy_z;
    DVelo(6) = vz_x; DVelo(7) = vz_y; DVelo(8) = vz_z;

    div_v = invF.MatTContraction(DVelo);

    matmodel->get_PK(F, P_iso, S_iso);
    rho = matmodel -> get_rho(p);
    mbeta = matmodel -> get_beta(p);

    if( std::abs(mbeta) < 1.0e-8 ) mbeta = 1.0;

    invF.VecMultT( p_x, p_y, p_z, gradP );

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      invF.VecMultT( NA_x, NA_y, NA_z, gradNA );

      Residual[4*A  ] += gwts * detF * NA * div_v;

      Residual[4*A+1] += gwts * ( NA_x * P_iso(0)
          + NA_y * P_iso(1) + NA_z * P_iso(2) - gradNA[0] * detF * p 
          - NA * rho * detF * fx );

      Residual[4*A+2] += gwts * ( NA_x * P_iso(3)
          + NA_y * P_iso(4) + NA_z * P_iso(5) - gradNA[1] * detF * p 
          - NA * rho * detF * fy );

      Residual[4*A+3] += gwts * ( NA_x * P_iso(6)
          + NA_y * P_iso(7) + NA_z * P_iso(8) - gradNA[2] * detF * p 
          - NA * rho * detF * fz );

      for(B=0; B<nLocBas; ++B)
      {
        index = A * nLocBas + B;
        NB = R[B];

        Sub_Tan[0][index]  += gwts * NA * mbeta * detF * NB;
        Sub_Tan[5][index]  += gwts * NA * rho * detF * NB;
        Sub_Tan[10][index] += gwts * NA * rho * detF * NB;
        Sub_Tan[15][index] += gwts * NA * rho * detF * NB;
      }
    }
  } // Finish-loop-over-qua

  for(ii=0; ii<4; ++ii)
  {
    for(jj=0; jj<4; ++jj)
    {
      for(A=0; A<nLocBas; ++A)
      {
        for(B=0; B<nLocBas; ++B)
        {
          Tangent[ 4*nLocBas*(4*A + ii) + 4*B + jj ]
            = Sub_Tan[ii*4+jj][A*nLocBas + B];
        }
      }
    }
  }
}


void PLocAssem_LPS0_Seg_GenAlpha::Assem_Residual_EBC(
    const int &ebc_id,
    const double &time, const double &dt,
    const double * const &velo,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element -> buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  int ii, qua, A;
  double gwts, coor_x, coor_y, coor_z, gx, gy, gz, nx, ny, nz, surface_area;
  const double curr = time + alpha_f * dt;

  Zero_Residual();

  for(qua=0; qua < face_nqp; ++qua)
  {
    element -> get_R(qua, R);
    element -> get_2d_normal_out(qua, nx, ny, nz, surface_area);

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
  } // finish loop over qua
}


// EOF
