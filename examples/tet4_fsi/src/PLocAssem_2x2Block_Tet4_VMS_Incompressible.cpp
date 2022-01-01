#include "PLocAssem_2x2Block_Tet4_VMS_Incompressible.hpp"

PLocAssem_2x2Block_Tet4_VMS_Incompressible::PLocAssem_2x2Block_Tet4_VMS_Incompressible(
    IMaterialModel * const &in_matmodel,
    const TimeMethod_GenAlpha * const &tm_gAlpha, const int &in_nqp )
: rho0( in_matmodel->get_elastic_rho0() ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()),
  nLocBas(4), vec_size_0( nLocBas * 3 ), vec_size_1( nLocBas ), 
  nqp(in_nqp), snLocBas(3), matmodel( in_matmodel )
{
  Tangent00 = new PetscScalar[vec_size_0 * vec_size_0];
  Tangent01 = new PetscScalar[vec_size_0 * vec_size_1];
  Tangent10 = new PetscScalar[vec_size_1 * vec_size_0];
  Tangent11 = new PetscScalar[vec_size_1 * vec_size_1];

  Residual0 = new PetscScalar[vec_size_0];
  Residual1 = new PetscScalar[vec_size_1];

  Zero_Tangent_Residual();

  print_info();
}

PLocAssem_2x2Block_Tet4_VMS_Incompressible::~PLocAssem_2x2Block_Tet4_VMS_Incompressible()
{
  delete [] Tangent00; Tangent00 = nullptr;
  delete [] Tangent01; Tangent01 = nullptr;
  delete [] Tangent10; Tangent10 = nullptr;
  delete [] Tangent11; Tangent11 = nullptr;

  delete [] Residual0; Residual0 = nullptr;
  delete [] Residual1; Residual1 = nullptr;
}

void PLocAssem_2x2Block_Tet4_VMS_Incompressible::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint("  Three-dimensional Hyper-elastic solid model:\n");
  SYS_T::commPrint("  Spatial: Finite element with VMS stabilization \n");
  SYS_T::commPrint("  Temporal: Generalized-alpha method \n");
  SYS_T::commPrint("  Solid density rho0 = %e g/cm3\n", rho0);
  matmodel->print_info();
  SYS_T::print_sep_line();
}

void PLocAssem_2x2Block_Tet4_VMS_Incompressible::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size_0; ++ii) Residual0[ii] = 0.0;
  for(int ii=0; ii<vec_size_1; ++ii) Residual1[ii] = 0.0;

  for(int ii=0; ii<vec_size_0 * vec_size_0; ++ii) Tangent00[ii] = 0.0;
  for(int ii=0; ii<vec_size_0 * vec_size_1; ++ii) Tangent01[ii] = 0.0;
  for(int ii=0; ii<vec_size_1 * vec_size_0; ++ii) Tangent10[ii] = 0.0;
  for(int ii=0; ii<vec_size_1 * vec_size_1; ++ii) Tangent11[ii] = 0.0;
}

void PLocAssem_2x2Block_Tet4_VMS_Incompressible::Zero_Residual()
{
  for(int ii=0; ii<vec_size_0; ++ii) Residual0[ii] = 0.0;
  for(int ii=0; ii<vec_size_1; ++ii) Residual1[ii] = 0.0;
}

void PLocAssem_2x2Block_Tet4_VMS_Incompressible::get_tau(
    double &tau_m_qua, double &tau_c_qua,
    const double &dt, const double &Jin, const double &dx ) const
{
  const double mu = matmodel->get_elastic_mu();
  //const double ka = matmodel->get_elastic_kappa();
  //const double c_max = std::pow( rho0 / (ka + 4*mu/3.0), -0.5);

  const double c_max = std::pow( rho0 / mu, -0.5); // Fully incompressible case

  const double dt_ka = dx / c_max;

  tau_m_qua = 1.0e-2 * dt_ka * Jin / rho0;
  tau_c_qua = 0.000 * dx * c_max * rho0 / Jin;
}

void PLocAssem_2x2Block_Tet4_VMS_Incompressible::Assem_Estimate()
{
  for(int ii=0; ii<vec_size_0 * vec_size_0; ++ii) Tangent00[ii] = 1.0;
  for(int ii=0; ii<vec_size_0 * vec_size_1; ++ii) Tangent01[ii] = 1.0;
  for(int ii=0; ii<vec_size_1 * vec_size_0; ++ii) Tangent10[ii] = 1.0;
  for(int ii=0; ii<vec_size_1 * vec_size_1; ++ii) Tangent11[ii] = 1.0;
}

void PLocAssem_2x2Block_Tet4_VMS_Incompressible::Assem_Residual(
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
      p   += pres[ii] * R[ii];
      p_x += pres[ii] * dR_dx[ii];
      p_y += pres[ii] * dR_dy[ii];
      p_z += pres[ii] * dR_dz[ii];

      ux_t += dot_disp[ii*3  ] * R[ii];
      uy_t += dot_disp[ii*3+1] * R[ii];
      uz_t += dot_disp[ii*3+2] * R[ii];
      p_t  += dot_pres[ii]     * R[ii];
      vx_t += dot_velo[ii*3  ] * R[ii];
      vy_t += dot_velo[ii*3+1] * R[ii];
      vz_t += dot_velo[ii*3+2] * R[ii];

      vx   += velo[ii*3  ] * R[ii];
      vy   += velo[ii*3+1] * R[ii];
      vz   += velo[ii*3+2] * R[ii];

      ux_x += disp[ii*3  ] * dR_dx[ii];
      uy_x += disp[ii*3+1] * dR_dx[ii];
      uz_x += disp[ii*3+2] * dR_dx[ii];

      ux_y += disp[ii*3  ] * dR_dy[ii];
      uy_y += disp[ii*3+1] * dR_dy[ii];
      uz_y += disp[ii*3+2] * dR_dy[ii];

      ux_z += disp[ii*3  ] * dR_dz[ii];
      uy_z += disp[ii*3+1] * dR_dz[ii];
      uz_z += disp[ii*3+2] * dR_dz[ii];

      vx_x += velo[ii*3  ] * dR_dx[ii];
      vy_x += velo[ii*3+1] * dR_dx[ii];
      vz_x += velo[ii*3+2] * dR_dx[ii];

      vx_y += velo[ii*3  ] * dR_dy[ii];
      vy_y += velo[ii*3+1] * dR_dy[ii];
      vz_y += velo[ii*3+2] * dR_dy[ii];

      vx_z += velo[ii*3  ] * dR_dz[ii];
      vy_z += velo[ii*3+1] * dR_dz[ii];
      vz_z += velo[ii*3+2] * dR_dz[ii];

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

      Residual1[A] += gwts * ( NA * Res_Mas + gradNA_ResMom );

      Residual0[3*A  ] += gwts * ( NA * rho * detF * (vx_t - f_body.x())
          + NA_x * P_iso(0) + NA_y * P_iso(1) + NA_z * P_iso(2)
          - gradNA.x() * (detF * p - tau_c * Res_Mas) );

      Residual0[3*A+1] += gwts * ( NA * rho * detF * (vy_t - f_body.y())
          + NA_x * P_iso(3) + NA_y * P_iso(4) + NA_z * P_iso(5)
          - gradNA.y() * (detF * p - tau_c * Res_Mas) );

      Residual0[3*A+2] += gwts * ( NA * rho * detF * (vz_t - f_body.z())
          + NA_x * P_iso(6) + NA_y * P_iso(7) + NA_z * P_iso(8)
          - gradNA.z() * (detF * p - tau_c * Res_Mas) );
    }
  }
}




// EOF
