#include "PLocAssem_Tet4_VMS_Seg_Debug.hpp"

PLocAssem_Tet4_VMS_Seg_Debug::PLocAssem_Tet4_VMS_Seg_Debug(
    IMaterialModel * const &in_matmodel,
    const TimeMethod_GenAlpha * const &tm_gAlpha,
    const int &in_nlocbas, const int &in_nqp,
    const int &in_snlocbas )
: rho0( in_matmodel->get_elastic_rho0() ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()),
  num_ebc_fun(6),  nLocBas(4), dof_per_node(7), vec_size(16),
  nqp(in_nqp), snLocBas(in_snlocbas)
{
  matmodel = in_matmodel;

  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size]; 
  
  Zero_Tangent_Residual();

  if( num_ebc_fun == 0 ) flist = NULL;
  else flist = new locassem_vms_seg_ela_fem_funs [num_ebc_fun];

  flist[0] = &PLocAssem_Tet4_VMS_Seg_Debug::get_bot_H;
  flist[1] = &PLocAssem_Tet4_VMS_Seg_Debug::get_top_H;
  flist[2] = &PLocAssem_Tet4_VMS_Seg_Debug::get_bac_H;
  flist[3] = &PLocAssem_Tet4_VMS_Seg_Debug::get_fro_H;
  flist[4] = &PLocAssem_Tet4_VMS_Seg_Debug::get_lef_H;
  flist[5] = &PLocAssem_Tet4_VMS_Seg_Debug::get_rig_H;

  print_info();
}


PLocAssem_Tet4_VMS_Seg_Debug::~PLocAssem_Tet4_VMS_Seg_Debug()
{
  delete [] Tangent; delete [] Residual; Tangent = NULL; Residual = NULL;
  if(num_ebc_fun > 0) delete [] flist;
}


void PLocAssem_Tet4_VMS_Seg_Debug::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Three-dimensional Segregated Debug code: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Spatial: Finite element with VMS stabilization \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Temporal: Generalized-alpha method \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Solid density rho0 = %e g/cm3\n\n", rho0);
  matmodel->print_info();
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PLocAssem_Tet4_VMS_Seg_Debug::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 0.0;
}


void PLocAssem_Tet4_VMS_Seg_Debug::Zero_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
}


void PLocAssem_Tet4_VMS_Seg_Debug::Assem_Estimate()
{
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 1.0;
}


void PLocAssem_Tet4_VMS_Seg_Debug::get_tau( 
    double &tau_m_qua, double &tau_c_qua,
    const double &dt, const double &Jin,
    const double &dx ) const
{
  tau_m_qua = 0.5 * dx * dx;
  tau_c_qua = 0.0;
}


void PLocAssem_Tet4_VMS_Seg_Debug::Assem_Residual(
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

  const double h_e = element->get_h( eleCtrlPts_x,
      eleCtrlPts_y, eleCtrlPts_z );

  int ii, qua, A, ii7;
  double p, p_t, p_x, p_y, p_z;
  double ux_t, uy_t, uz_t, vx, vy, vz;
  double vx_t, vy_t, vz_t;
  double ux_x, uy_x, uz_x, ux_y, uy_y, uz_y, ux_z, uy_z, uz_z;
  double vx_x, vy_x, vz_x, vx_y, vy_y, vz_y, vx_z, vy_z, vz_z;
  double fc, fx, fy, fz, gwts, coor_x, coor_y, coor_z;
  double NA, NA_x, NA_y, NA_z;

  double Res_Mom[3];
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
    get_fc(coor_x, coor_y, coor_z, curr, fc);

    // Stokes test
    get_tau(tau_m, tau_c, dt, detF, h_e);
    Res_Mom[0]= vx_t + p_x - fx;
    Res_Mom[1]= vy_t + p_y - fy;
    Res_Mom[2]= vz_t + p_z - fz;

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];
     
      Residual[4*A+0] += gwts * ( NA * (p_t + vx_x + vy_y + vz_z - fc)
          + tau_m * (NA_x * Res_Mom[0]+NA_y * Res_Mom[1]
            + NA_z * Res_Mom[2]) );

      Residual[4*A+1] += gwts * (NA * vx_t - NA_x * p
          + NA_x * vx_x + NA_y * vx_y + NA_z * vx_z - NA * fx );

      Residual[4*A+2] += gwts * (NA * vy_t - NA_y * p 
          + NA_x * vy_x + NA_y * vy_y + NA_z * vy_z - NA * fy );
      
      Residual[4*A+3] += gwts * (NA * vz_t - NA_z * p 
          + NA_x * vz_x + NA_y * vz_y + NA_z * vz_z - NA * fz );
    }
  }
}


void PLocAssem_Tet4_VMS_Seg_Debug::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &velo,
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y,
      eleCtrlPts_z );

  const double h_e = element->get_h( eleCtrlPts_x,
      eleCtrlPts_y, eleCtrlPts_z );

  int ii, jj, qua, A, B, index, ii7;
  double p, p_t, p_x, p_y, p_z;
  double ux_t, uy_t, uz_t, vx, vy, vz;
  double vx_t, vy_t, vz_t;
  double ux_x, uy_x, uz_x, ux_y, uy_y, uz_y, ux_z, uy_z, uz_z;
  double vx_x, vy_x, vz_x, vx_y, vy_y, vz_y, vx_z, vy_z, vz_z;
  double fc, fx, fy, fz, gwts, coor_x, coor_y, coor_z;
  double NA, NA_x, NA_y, NA_z;
  double NB, NB_x, NB_y, NB_z;
  double Res_Mom[3];
  //double Res_Mas;
  const double curr = time + alpha_f * dt;

  const double dd_dv = alpha_f * gamma * dt;
  
  Zero_Tangent_Residual();

  Zero_Sub_Tan();

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
    get_fc(coor_x, coor_y, coor_z, curr, fc);

    get_tau(tau_m, tau_c, dt, detF, h_e);

    // Residual of momentum equation
    Res_Mom[0]= vx_t + p_x - fx;
    Res_Mom[1]= vy_t + p_y - fy;
    Res_Mom[2]= vz_t + p_z - fz;

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      Residual[4*A+0] += gwts * ( NA * (p_t + vx_x + vy_y + vz_z - fc)
          + tau_m * (NA_x * Res_Mom[0] + NA_y * Res_Mom[1]
            + NA_z * Res_Mom[2]) );

      Residual[4*A+1] += gwts * (NA * vx_t - NA_x * p
          + NA_x * vx_x + NA_y * vx_y + NA_z * vx_z - NA * fx );

      Residual[4*A+2] += gwts * (NA * vy_t - NA_y * p
          + NA_x * vy_x + NA_y * vy_y + NA_z * vy_z - NA * fy );

      Residual[4*A+3] += gwts * (NA * vz_t - NA_z * p
          + NA_x * vz_x + NA_y * vz_y + NA_z * vz_z - NA * fz );

      for(B=0; B<nLocBas; ++B)
      {
        index = A * nLocBas + B;
        NB = R[B]; NB_x = dR_dx[B]; NB_y = dR_dy[B]; NB_z = dR_dz[B];

        Sub_Tan[0][index] += gwts * ( alpha_m * NA * NB 
            + dd_dv * tau_m * ( NA_x * NB_x 
            + NA_y * NB_y + NA_z * NB_z ) );

        Sub_Tan[1][index] += gwts * ( alpha_m * NA_x * tau_m * NB
            + dd_dv * NA * NB_x );

        Sub_Tan[2][index] += gwts * ( alpha_m * NA_y * tau_m * NB
            + dd_dv * NA * NB_y );

        Sub_Tan[3][index] += gwts * ( alpha_m * NA_z * tau_m * NB
            + dd_dv * NA * NB_z );

        Sub_Tan[4][index] -= gwts * dd_dv * NA_x * NB;
        Sub_Tan[8][index] -= gwts * dd_dv * NA_y * NB;
        Sub_Tan[12][index]-= gwts * dd_dv * NA_z * NB;

        const double mass_entry = gwts * NA * NB * alpha_m;

        Sub_Tan[5][index] += mass_entry 
          + gwts * dd_dv * (NA_x * NB_x + NA_y * NB_y + NA_z * NB_z);

        Sub_Tan[10][index] += mass_entry
          + gwts * dd_dv * (NA_x * NB_x + NA_y * NB_y + NA_z * NB_z);

        Sub_Tan[15][index] += mass_entry
          + gwts * dd_dv * (NA_x * NB_x + NA_y * NB_y + NA_z * NB_z);

      } // Finish Loop-B
    } // Finish Loop-A
  } // Finish Loop-qua

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


void PLocAssem_Tet4_VMS_Seg_Debug::Assem_Mass_Residual(
    const double * const &disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y,
      eleCtrlPts_z );

  int ii, jj, qua, A, B, index, ii7;
  double p, ux_x, uy_x, uz_x, ux_y, uy_y, uz_y, ux_z, uy_z, uz_z;
  double vx_x, vy_x, vz_x, vx_y, vy_y, vz_y, vx_z, vy_z, vz_z;
  double vx, vy, vz;
  double fc, fx, fy, fz, gwts, coor_x, coor_y, coor_z;
  double NA, NB, NA_x, NA_y, NA_z;

  double curr = 0.0;

  Zero_Tangent_Residual();

  Zero_Sub_Tan();

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
    get_fc(coor_x, coor_y, coor_z, curr, fc);

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      Residual[4*A+0] += gwts * NA * (vx_x + vy_y + vz_z - fc);

      Residual[4*A+1] += gwts * ( -NA_x * p
          + NA_x * vx_x + NA_y * vx_y + NA_z * vx_z - NA * fx );

      Residual[4*A+2] += gwts * ( -NA_y * p
          + NA_x * vy_x + NA_y * vy_y + NA_z * vy_z - NA * fy );

      Residual[4*A+3] += gwts * ( -NA_z * p
          + NA_x * vz_x + NA_y * vz_y + NA_z * vz_z - NA * fz );

      for(B=0; B<nLocBas; ++B)
      {
        index = A * nLocBas + B;
        NB = R[B];
        Sub_Tan[0][index]  += gwts * NA * NB;
        Sub_Tan[5][index]  += gwts * NA * NB;
        Sub_Tan[10][index] += gwts * NA * NB;
        Sub_Tan[15][index] += gwts * NA * NB;
      } // Finish loop-B
    } // Finish loop-A
  } // Finish loop-qua

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


void PLocAssem_Tet4_VMS_Seg_Debug::Assem_Residual_EBC(
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
