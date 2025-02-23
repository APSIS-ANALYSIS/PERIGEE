#include "PLocAssem_2x2Block_VMS_Incompressible.hpp"

PLocAssem_2x2Block_VMS_Incompressible::PLocAssem_2x2Block_VMS_Incompressible(
    const FEType &in_type, const int &in_nqp_v, const int &in_nqp_s,
    const TimeMethod_GenAlpha * const &tm_gAlpha,
    std::unique_ptr<MaterialModel_Mixed_Elasticity> in_matmodel )
: elemType(in_type), nqpv(in_nqp_v), nqps(in_nqp_s),
  elementv( ElementFactory::createVolElement(elemType, nqpv) ),
  elements( ElementFactory::createSurElement(elemType, nqps) ),
  quadv( QuadPtsFactory::createVolQuadrature(elemType, nqpv) ),
  quads( QuadPtsFactory::createSurQuadrature(elemType, nqps) ),
: rho0( in_matmodel->get_rho_0() ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()),
  nLocBas( elementv->get_nLocBas() ), snLocBas( elments->get_nLocBas() ), 
  vec_size_0( nLocBas * 3 ), vec_size_1( nLocBas ),
  sur_size_0( snLocBas * 3 ),
  matmodel( std::move(in_matmodel) )
{
  Tangent00 = new PetscScalar[vec_size_0 * vec_size_0];
  Tangent01 = new PetscScalar[vec_size_0 * vec_size_1];
  Tangent10 = new PetscScalar[vec_size_1 * vec_size_0];
  Tangent11 = new PetscScalar[vec_size_1 * vec_size_1];

  Residual0 = new PetscScalar[vec_size_0];
  Residual1 = new PetscScalar[vec_size_1];

  sur_Residual0 = new PetscScalar[sur_size_0];

  Zero_Tangent_Residual();
  Zero_sur_Residual();

  print_info();
}

PLocAssem_2x2Block_VMS_Incompressible::~PLocAssem_2x2Block_VMS_Incompressible()
{
  delete [] Tangent00; Tangent00 = nullptr;
  delete [] Tangent01; Tangent01 = nullptr;
  delete [] Tangent10; Tangent10 = nullptr;
  delete [] Tangent11; Tangent11 = nullptr;

  delete [] Residual0; Residual0 = nullptr;
  delete [] Residual1; Residual1 = nullptr;

  delete [] sur_Residual0; sur_Residual0 = nullptr;
}

void PLocAssem_2x2Block_VMS_Incompressible::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint("  Three-dimensional Hyper-elastic solid model:\n");
  elementv->print_info();
  SYS_T::commPrint("  Spatial: Finite element with VMS stabilization \n");
  SYS_T::commPrint("  Temporal: Generalized-alpha method \n");
  SYS_T::commPrint("  Solid density rho0 = %e g/cm3\n", rho0);
  matmodel->print_info();
  SYS_T::print_sep_line();
}

void PLocAssem_2x2Block_VMS_Incompressible::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size_0; ++ii) Residual0[ii] = 0.0;
  for(int ii=0; ii<vec_size_1; ++ii) Residual1[ii] = 0.0;

  for(int ii=0; ii<vec_size_0 * vec_size_0; ++ii) Tangent00[ii] = 0.0;
  for(int ii=0; ii<vec_size_0 * vec_size_1; ++ii) Tangent01[ii] = 0.0;
  for(int ii=0; ii<vec_size_1 * vec_size_0; ++ii) Tangent10[ii] = 0.0;
  for(int ii=0; ii<vec_size_1 * vec_size_1; ++ii) Tangent11[ii] = 0.0;
}

void PLocAssem_2x2Block_VMS_Incompressible::Zero_Residual()
{
  for(int ii=0; ii<vec_size_0; ++ii) Residual0[ii] = 0.0;
  for(int ii=0; ii<vec_size_1; ++ii) Residual1[ii] = 0.0;
}

std::array<double, 2> PLocAssem_2x2Block_VMS_Incompressible::get_tau(
    const double &dt, const double &Jin, const double &dx ) const
{
  const double mu = matmodel->get_elastic_mu();

  const double c_max = std::pow( rho0 / mu, -0.5); // Fully incompressible case

  const double dt_ka = dx / c_max;

  return {{1.0e-2 * dt_ka * Jin / rho0, 0.000 * dx * c_max * rho0 / Jin}};
}

void PLocAssem_2x2Block_VMS_Incompressible::Assem_Estimate()
{
  for(int ii=0; ii<vec_size_0 * vec_size_0; ++ii) Tangent00[ii] = 1.0;
  for(int ii=0; ii<vec_size_0 * vec_size_1; ++ii) Tangent01[ii] = 1.0;
  for(int ii=0; ii<vec_size_1 * vec_size_0; ++ii) Tangent10[ii] = 1.0;
  for(int ii=0; ii<vec_size_1 * vec_size_1; ++ii) Tangent11[ii] = 1.0;
}

void PLocAssem_2x2Block_VMS_Incompressible::Assem_Residual(
    const double &time, const double &dt,
    const double * const &dot_disp,
    const double * const &dot_velo,
    const double * const &dot_pres,
    const double * const &disp,
    const double * const &velo,
    const double * const &pres,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &qua_prestress )
{
  elementv->buildBasis( quadv.get(), eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double h_e = elementv->get_h( eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double curr = time + alpha_f * dt;

  Zero_Residual();

  for(int qua=0; qua < nqpv; ++qua)
  {
    double p = 0.0, p_x = 0.0, p_y = 0.0, p_z = 0.0;

    double vx_t = 0.0, vy_t = 0.0, vz_t = 0.0;

    double ux_x = 0.0, uy_x = 0.0, uz_x = 0.0;
    double ux_y = 0.0, uy_y = 0.0, uz_y = 0.0;
    double ux_z = 0.0, uy_z = 0.0, uz_z = 0.0;

    double vx_x = 0.0, vy_x = 0.0, vz_x = 0.0;
    double vx_y = 0.0, vy_y = 0.0, vz_y = 0.0;
    double vx_z = 0.0, vy_z = 0.0, vz_z = 0.0;

    Vector_3 coor(0.0, 0.0, 0.0);

    std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

    elementv->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      p   += pres[ii] * R[ii];
      p_x += pres[ii] * dR_dx[ii];
      p_y += pres[ii] * dR_dy[ii];
      p_z += pres[ii] * dR_dz[ii];

      vx_t += dot_velo[ii*3  ] * R[ii];
      vy_t += dot_velo[ii*3+1] * R[ii];
      vz_t += dot_velo[ii*3+2] * R[ii];

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

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    const double gwts = elementv->get_detJac(qua) * quadv->get_qw(qua);

    const Vector_3 f_body = get_f(coor, curr);

    const Tensor2_3D F( ux_x + 1.0, ux_y, ux_z, uy_x, uy_y + 1.0, uy_z, uz_x, uz_y, uz_z + 1.0 );

    const Tensor2_3D invF = Ten2::inverse(F);

    const Tensor2_3D DVelo( vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z );

    const double invFDV_t = invF.MatTContraction(DVelo); // invF_Ii V_i,I

    Tensor2_3D P_iso = matmodel->get_PK_1st( F );

    // ------------------------------------------------------------------------
    // 1st PK stress corrected by prestress
    const Tensor2_3D prestress( qua_prestress[qua*6+0], qua_prestress[qua*6+5], qua_prestress[qua*6+4],
        qua_prestress[qua*6+5], qua_prestress[qua*6+1], qua_prestress[qua*6+3],
        qua_prestress[qua*6+4], qua_prestress[qua*6+3], qua_prestress[qua*6+2] );

    P_iso += prestress * Ten2::cofactor( F );
    // ------------------------------------------------------------------------

    const double rho = matmodel->get_rho(p);

    const double detF = F.det();

    // Get stabilization parameters
    const std::array<double, 2> tau = get_tau(dt, detF, h_e);
    const double tau_m = tau[0];
    const double tau_c = tau[1];

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

void PLocAssem_2x2Block_VMS_Incompressible::Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_disp,
        const double * const &dot_velo,
        const double * const &dot_pres,
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &qua_prestress )
{
  elementv->buildBasis( quadv.get(), eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double h_e = elementv->get_h( eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double curr = time + alpha_f * dt;

  const double dd_dv = alpha_f * gamma * dt;
  
  const double ddvm = dd_dv * dd_dv / alpha_m;

  Zero_Tangent_Residual();

  for(int qua=0; qua < nqpv; ++qua)
  {
    double p = 0.0, p_x = 0.0, p_y = 0.0, p_z = 0.0;

    double vx_t = 0.0, vy_t = 0.0, vz_t = 0.0;

    double ux_x = 0.0, uy_x = 0.0, uz_x = 0.0;
    double ux_y = 0.0, uy_y = 0.0, uz_y = 0.0;
    double ux_z = 0.0, uy_z = 0.0, uz_z = 0.0;

    double vx_x = 0.0, vy_x = 0.0, vz_x = 0.0;
    double vx_y = 0.0, vy_y = 0.0, vz_y = 0.0;
    double vx_z = 0.0, vy_z = 0.0, vz_z = 0.0;

    Vector_3 coor(0.0, 0.0, 0.0);

    std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

    elementv->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      p   += pres[ii] * R[ii];
      p_x += pres[ii] * dR_dx[ii];
      p_y += pres[ii] * dR_dy[ii];
      p_z += pres[ii] * dR_dz[ii];

      vx_t += dot_velo[ii*3  ] * R[ii];
      vy_t += dot_velo[ii*3+1] * R[ii];
      vz_t += dot_velo[ii*3+2] * R[ii];

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

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    const double gwts = elementv->get_detJac(qua) * quadv->get_qw(qua);

    const Vector_3 f_body = get_f(coor, curr);

    const Tensor2_3D F( ux_x + 1.0, ux_y, ux_z, uy_x, uy_y + 1.0, uy_z, uz_x, uz_y, uz_z + 1.0 );

    const Tensor2_3D invF = Ten2::inverse(F);

    const Tensor2_3D DVelo( vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z );

    const Tensor2_3D Dvelo_invF = DVelo * invF;  // v_i,I invF_Ij = v_i,j

    double GradP_invF[3];
    invF.VecMultT( p_x, p_y, p_z, GradP_invF[0], GradP_invF[1], GradP_invF[2] ); // p_I invF_ii = p,i

    const double invFDV_t = invF.MatTContraction(DVelo); // invF_Ii V_i,I

    Tensor2_3D P_iso;
    SymmTensor2_3D S_iso;
    const Tensor4_3D AA_iso = matmodel->get_PK_FFStiffness(F, P_iso, S_iso);

    // ------------------------------------------------------------------------
    // 1st PK stress corrected by prestress
    const Tensor2_3D prestress( qua_prestress[qua*6+0], qua_prestress[qua*6+5], qua_prestress[qua*6+4],
        qua_prestress[qua*6+5], qua_prestress[qua*6+1], qua_prestress[qua*6+3],
        qua_prestress[qua*6+4], qua_prestress[qua*6+3], qua_prestress[qua*6+2] );

    P_iso += prestress * Ten2::cofactor( F );
    // ------------------------------------------------------------------------
    
    const double rho = matmodel->get_rho(p);
    const double detF = F.det();

    // Get stabilization parameters
    const std::array<double, 2> tau = get_tau(dt, detF, h_e);
    const double tau_m = tau[0];
    const double tau_c = tau[1];

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

      Residual1[A] += gwts * ( NA * Res_Mas + GradNA_invF_ResMom );

      Residual0[3*A  ] += gwts * ( NA * rho * detF * (vx_t - f_body.x())
          + NA_x * P_iso(0) + NA_y * P_iso(1) + NA_z * P_iso(2) 
          - GradNA_invF[0] * detF * p );

      Residual0[3*A+1] += gwts * ( NA * rho * detF * (vy_t - f_body.y())
          + NA_x * P_iso(3) + NA_y * P_iso(4) + NA_z * P_iso(5) 
          - GradNA_invF[1] * detF * p );

      Residual0[3*A+2] += gwts * ( NA * rho * detF * (vz_t - f_body.z())
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

        Tangent11[ nLocBas*A + B ] += gwts * dd_dv * GradNA_invF_dot_GradNB_invF * tau_m * detF;

        Tangent10[ 3*nLocBas*A + 3*B ] += gwts * ( 
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

        Tangent10[ 3*nLocBas*A + 3*B + 1 ] += gwts * (
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

        Tangent10[ 3*nLocBas*A + 3*B + 2 ] += gwts * (
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

        Tangent01[nLocBas*(3*A+0)+B] -= gwts * GradNA_invF[0] * detF * NB * dd_dv;

        Tangent01[nLocBas*(3*A+1)+B] -= gwts * GradNA_invF[1] * detF * NB * dd_dv;

        Tangent01[nLocBas*(3*A+2)+B] -= gwts * GradNA_invF[2] * detF * NB * dd_dv;

        const double mass_entry = gwts * NA * NB * rho * detF * alpha_m;

        Tangent00[3*nLocBas*(3*A+0)+3*B+0] += mass_entry
          + gwts * dd_dv * tau_c * detF * GradNA_invF[0] * GradNB_invF[0];

        Tangent00[3*nLocBas*(3*A+0)+3*B+1] += gwts * dd_dv * tau_c * detF * GradNA_invF[0] * GradNB_invF[1];

        Tangent00[3*nLocBas*(3*A+0)+3*B+2] += gwts * dd_dv * tau_c * detF * GradNA_invF[0] * GradNB_invF[2];

        Tangent00[3*nLocBas*(3*A+1)+3*B+0] += gwts * dd_dv * tau_c * detF * GradNA_invF[1] * GradNB_invF[0];

        Tangent00[3*nLocBas*(3*A+1)+3*B+1] += mass_entry
          + gwts * dd_dv * tau_c * detF * GradNA_invF[1] * GradNB_invF[1];

        Tangent00[3*nLocBas*(3*A+1)+3*B+2] += gwts * dd_dv * tau_c * detF * GradNA_invF[1] * GradNB_invF[2];

        Tangent00[3*nLocBas*(3*A+2)+3*B+0] += gwts * dd_dv * tau_c * detF * GradNA_invF[2] * GradNB_invF[0];

        Tangent00[3*nLocBas*(3*A+2)+3*B+1] += gwts * dd_dv * tau_c * detF * GradNA_invF[2] * GradNB_invF[1];

        Tangent00[3*nLocBas*(3*A+2)+3*B+2] += mass_entry
          + gwts * dd_dv * tau_c * detF * GradNA_invF[2] * GradNB_invF[2];

        const double geo_stiff = gwts * ddvm * (
            NA_x * ( S_iso.xx() * NB_x + S_iso.xy() * NB_y + S_iso.xz() * NB_z )
            + NA_y * ( S_iso.yx() * NB_x + S_iso.yy() * NB_y + S_iso.yz() * NB_z )
            + NA_z * ( S_iso.zx() * NB_x + S_iso.zy() * NB_y + S_iso.zz() * NB_z) );

        Tangent00[3*nLocBas*(3*A+0)+3*B+0] += geo_stiff;

        Tangent00[3*nLocBas*(3*A+1)+3*B+1] += geo_stiff;

        Tangent00[3*nLocBas*(3*A+2)+3*B+2] += geo_stiff;

        for(int ii=0; ii<3; ++ii)
        {
          for(int jj=0; jj<3; ++jj)
          {
            Tangent00[ 3*nLocBas*(3*A + ii) + 3*B + jj ] += gwts * ddvm * (
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

void PLocAssem_2x2Block_VMS_Incompressible::Assem_Mass_Residual(
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &qua_prestress )
{
  elementv->buildBasis( quadv.get(), eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double curr = 0.0;

  Zero_Tangent_Residual();

  for(int qua=0; qua<nqpv; ++qua)
  {
    double p = 0.0;
    double ux_x = 0.0, uy_x = 0.0, uz_x = 0.0;
    double ux_y = 0.0, uy_y = 0.0, uz_y = 0.0;
    double ux_z = 0.0, uy_z = 0.0, uz_z = 0.0;

    double vx_x = 0.0, vy_x = 0.0, vz_x = 0.0;
    double vx_y = 0.0, vy_y = 0.0, vz_y = 0.0;
    double vx_z = 0.0, vy_z = 0.0, vz_z = 0.0;

    Vector_3 coor(0.0, 0.0, 0.0);

    std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

    elementv->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      p   += pres[ii] * R[ii];

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

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    const double gwts = elementv->get_detJac(qua) * quadv->get_qw(qua);
    
    const Vector_3 f_body = get_f(coor, curr);

    const Tensor2_3D F( ux_x + 1.0, ux_y, ux_z, uy_x, uy_y + 1.0, uy_z, uz_x, uz_y, uz_z + 1.0 );

    const Tensor2_3D invF = Ten2::inverse(F);

    const Tensor2_3D DVelo(  vx_x, vx_y, vx_z, vy_x, vy_y, vy_z, vz_x, vz_y, vz_z );

    // invF_Ii DV_i,I = v_i,i = div v
    const double invFDV_t = invF.MatTContraction(DVelo);

    Tensor2_3D P_iso = matmodel->get_PK_1st( F );

    // ------------------------------------------------------------------------
    // 1st PK stress corrected by prestress
    const Tensor2_3D prestress( qua_prestress[qua*6+0], qua_prestress[qua*6+5], qua_prestress[qua*6+4],
        qua_prestress[qua*6+5], qua_prestress[qua*6+1], qua_prestress[qua*6+3],
        qua_prestress[qua*6+4], qua_prestress[qua*6+3], qua_prestress[qua*6+2] );

    P_iso += prestress * Ten2::cofactor( F );
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

      Residual1[A] += gwts * NA * detF * invFDV_t;

      Residual0[3*A  ] += gwts * ( NA_x * P_iso(0) + NA_y * P_iso(1) 
          + NA_z * P_iso(2) - gradNA.x() * detF * p 
          - NA * rho * detF * f_body.x() );

      Residual0[3*A+1] += gwts * ( NA_x * P_iso(3) + NA_y * P_iso(4) 
          + NA_z * P_iso(5) - gradNA.y() * detF * p 
          - NA * rho * detF * f_body.y() );

      Residual0[3*A+2] += gwts * ( NA_x * P_iso(6) + NA_y * P_iso(7) 
          + NA_z * P_iso(8) - gradNA.z() * detF * p 
          - NA * rho * detF * f_body.z() );

      for(int B=0; B<nLocBas; ++B)
      {
        Tangent11[nLocBas*A + B]   += gwts * NA * detF * mbeta * R[B];
        Tangent00[3*nLocBas*(3*A+0) + 3*B+0] += gwts * NA * detF * rho * R[B];
        Tangent00[3*nLocBas*(3*A+1) + 3*B+1] += gwts * NA * detF * rho * R[B];
        Tangent00[3*nLocBas*(3*A+2) + 3*B+2] += gwts * NA * detF * rho * R[B];
      } // Finish loop-B
    } // Finish loop-A
  } // Finish loop-qua
}

void PLocAssem_2x2Block_VMS_Incompressible::Assem_Residual_EBC(
    const int &ebc_id,
    const double &time, const double &dt,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z )
{
  elements->buildBasis( quads.get(), eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double curr = time + alpha_f * dt;

  Zero_sur_Residual();

  for(int qua = 0; qua < nqps; ++qua)
  {
    const std::vector<double> R = elements->get_R(qua);

    double surface_area;
    const Vector_3 n_out = elements->get_2d_normal_out(qua, surface_area);

    Vector_3 coor(0.0, 0.0, 0.0);
    for(int ii=0; ii<snLocBas; ++ii)
    {
      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    const Vector_3 gg = get_ebc_fun( ebc_id, coor, curr, n_out );

    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual0[3*A  ] -= surface_area * quads -> get_qw(qua) * R[A] * gg.x();
      sur_Residual0[3*A+1] -= surface_area * quads -> get_qw(qua) * R[A] * gg.y();
      sur_Residual0[3*A+2] -= surface_area * quads -> get_qw(qua) * R[A] * gg.z();
    }
  }
}

void PLocAssem_2x2Block_VMS_Incompressible::Assem_Residual_Interior_Wall_EBC(
    const double &time,
    const double * const &pres,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z )
{
  const double factor = 1.0; //time >= 1.0 ? 1.0 : time;

  elements->buildBasis( quads.get(), eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  Zero_sur_Residual();

  for(int qua = 0; qua < nqps; ++qua)
  {
    const std::vector<double> R = elements->get_R(qua);

    double surface_area;
    const Vector_3 n_out = elements->get_2d_normal_out(qua, surface_area);

    double pp = 0.0;
    for(int ii=0; ii<snLocBas; ++ii) pp += pres[ii] * R[ii];

    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual0[3*A  ] -= surface_area * quads -> get_qw(qua) * R[A] * factor * pp * n_out.x();
      sur_Residual0[3*A+1] -= surface_area * quads -> get_qw(qua) * R[A] * factor * pp * n_out.y();
      sur_Residual0[3*A+2] -= surface_area * quads -> get_qw(qua) * R[A] * factor * pp * n_out.z();
    }
  }
}

std::vector<SymmTensor2_3D> PLocAssem_2x2Block_VMS_Incompressible::get_Wall_CauchyStress(
    const double * const &disp,
    const double * const &pres,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z ) const
{
  elementv->buildBasis( quadv.get(), eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  std::vector<SymmTensor2_3D> stress( nqp );

  for( int qua = 0; qua < nqpv; ++qua )
  {
    std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

    elementv->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );

    double pp = 0.0;
    double ux_x = 0.0, uy_x = 0.0, uz_x = 0.0;
    double ux_y = 0.0, uy_y = 0.0, uz_y = 0.0;
    double ux_z = 0.0, uy_z = 0.0, uz_z = 0.0;

    for(int ii=0; ii<nLocBas; ++ii)
    {
      pp   += pres[ii] * R[ii];

      ux_x += disp[ii*3+0] * dR_dx[ii];
      uy_x += disp[ii*3+1] * dR_dx[ii];
      uz_x += disp[ii*3+2] * dR_dx[ii];

      ux_y += disp[ii*3+0] * dR_dy[ii];
      uy_y += disp[ii*3+1] * dR_dy[ii];
      uz_y += disp[ii*3+2] * dR_dy[ii];

      ux_z += disp[ii*3+0] * dR_dz[ii];
      uy_z += disp[ii*3+1] * dR_dz[ii];
      uz_z += disp[ii*3+2] * dR_dz[ii];
    }

    const Tensor2_3D F( ux_x + 1.0, ux_y, ux_z, uy_x, uy_y + 1.0, uy_z, uz_x, uz_y, uz_z + 1.0 );

    stress[qua] = matmodel -> get_Cauchy_stress( F );

    stress[qua].xx() -= pp;
    stress[qua].yy() -= pp;
    stress[qua].zz() -= pp;
  }

  return stress;
}

// EOF
