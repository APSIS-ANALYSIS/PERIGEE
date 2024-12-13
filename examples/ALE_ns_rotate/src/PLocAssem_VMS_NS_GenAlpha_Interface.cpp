#include "PLocAssem_VMS_NS_GenAlpha_Interface.hpp"

PLocAssem_VMS_NS_GenAlpha_Interface::PLocAssem_VMS_NS_GenAlpha_Interface(
  const TimeMethod_GenAlpha * const &tm_gAlpha,
  const int &in_nlocbas, const int &in_nqp,
  const int &in_snlocbas, const double &in_rho, 
  const double &in_vis_mu, const double &in_beta,
  const int &elemtype, const double &angular,
  const Vector_3 &point_xyz, const Vector_3 &angular_direc,
  const double &in_ct, const double &in_ctauc,
  const double &in_C_bI)
  : PLocAssem_VMS_NS_GenAlpha_WeakBC(tm_gAlpha, in_nlocbas, in_nqp, in_snlocbas,
  in_rho, in_vis_mu, in_beta, elemtype, angular, point_xyz, angular_direc, in_ct, in_ctauc, in_C_bI)
{
  Tangent_ss = new PetscScalar[vec_size * vec_size];
  Tangent_sr = new PetscScalar[vec_size * vec_size];
  Tangent_rs = new PetscScalar[vec_size * vec_size];
  Tangent_rr = new PetscScalar[vec_size * vec_size];

  Residual_s = new PetscScalar[vec_size];
  Residual_r = new PetscScalar[vec_size];

  Zero_Residual_s();
  Zero_Residual_r();
  Zero_Tangent_ss();
  Zero_Tangent_sr();
  Zero_Tangent_rs();
  Zero_Tangent_rr();
}

PLocAssem_VMS_NS_GenAlpha_Interface::~PLocAssem_VMS_NS_GenAlpha_Interface()
{
  delete [] Tangent_ss; Tangent_ss = nullptr;
  delete [] Tangent_sr; Tangent_sr = nullptr;
  delete [] Tangent_rs; Tangent_rs = nullptr;
  delete [] Tangent_rr; Tangent_rr = nullptr;

  delete [] Residual_s; Residual_s = nullptr;
  delete [] Residual_r; Residual_r = nullptr;
}

void PLocAssem_VMS_NS_GenAlpha_Interface::print_info() const
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
  SYS_T::commPrint("  Weakly enforced Dirichlet BC was applied on wall.\n");
  SYS_T::commPrint("  Nitsche's method was applied on interfaces.\n");
  SYS_T::commPrint("  Note: \n");
  SYS_T::commPrint("  1. Consistent tangent matrix used. \n");
  SYS_T::commPrint("  2. Nonlinear quadratic term is in advective form. \n");
  SYS_T::commPrint("  3. Pressure is evaluated at n+alpha_f rather than n+1. \n");
  SYS_T::commPrint("----------------------------------------------------------- \n");
}

void PLocAssem_VMS_NS_GenAlpha_Interface::Assem_Residual_itf_fixed(
  const int &qua, const double &fixed_qw, const double &dt,
  const FEAElement * const &fixed_elementv, const FEAElement * const &rotated_elementv,
  const double * const &fixed_local_sol, const double * const &rotated_local_sol,
  const double * const &rotated_local_mvelo)
{
  Zero_Residual_s();
  double ps {0.0};
  double us {0.0}, us_x {0.0}, us_y {0.0}, us_z {0.0};
  double vs {0.0}, vs_x {0.0}, vs_y {0.0}, vs_z {0.0};
  double ws {0.0}, ws_x {0.0}, ws_y {0.0}, ws_z {0.0};

  double pr {0.0};
  double ur {0.0}, ur_x {0.0}, ur_y {0.0}, ur_z {0.0};
  double vr {0.0}, vr_x {0.0}, vr_y {0.0}, vr_z {0.0};
  double wr {0.0}, wr_x {0.0}, wr_y {0.0}, wr_z {0.0};
  double mur {0.0}, mvr {0.0}, mwr{0.0};

  std::vector<double> Ns(nLocBas, 0.0), dNs_dx(nLocBas, 0.0), dNs_dy(nLocBas, 0.0), dNs_dz(nLocBas, 0.0);
  std::vector<double> Nr(nLocBas, 0.0), dNr_dx(nLocBas, 0.0), dNr_dy(nLocBas, 0.0), dNr_dz(nLocBas, 0.0);

  fixed_elementv -> get_R_gradR( qua, &Ns[0], &dNs_dx[0], &dNs_dy[0], &dNs_dz[0] );
  rotated_elementv -> get_R_gradR( 0, &Nr[0], &dNr_dx[0], &dNr_dy[0], &dNr_dz[0] );

  double fixed_J {0.0};
  const Vector_3 normal_s = fixed_elementv -> get_2d_normal_out(qua, fixed_J);
  const Vector_3 normal_r = -1 * normal_s;

  // Calculate h_b and tau_I
  const auto s_dxi_dx = fixed_elementv -> get_invJacobian(qua);
  const double h_s = get_h_b(s_dxi_dx, normal_s);

  const auto r_dxi_dx = rotated_elementv -> get_invJacobian(0);
  const double h_r = get_h_b(r_dxi_dx, normal_r);

  const double tau_I = 0.5 * vis_mu * C_bI * (1.0 / h_s + 1.0 / h_r);

  for(int ii{0}; ii<nLocBas; ++ii)
  {
    const int ii4{4 * ii};
    const int ii3{3 * ii};

    ps += fixed_local_sol[ii4 + 0] * Ns[ii];
    us += fixed_local_sol[ii4 + 1] * Ns[ii];
    vs += fixed_local_sol[ii4 + 2] * Ns[ii];
    ws += fixed_local_sol[ii4 + 3] * Ns[ii];

    us_x += fixed_local_sol[ii4 + 1] * dNs_dx[ii];
    us_y += fixed_local_sol[ii4 + 1] * dNs_dy[ii];
    us_z += fixed_local_sol[ii4 + 1] * dNs_dz[ii];

    vs_x += fixed_local_sol[ii4 + 2] * dNs_dx[ii];
    vs_y += fixed_local_sol[ii4 + 2] * dNs_dy[ii];
    vs_z += fixed_local_sol[ii4 + 2] * dNs_dz[ii];

    ws_x += fixed_local_sol[ii4 + 3] * dNs_dx[ii];
    ws_y += fixed_local_sol[ii4 + 3] * dNs_dy[ii];
    ws_z += fixed_local_sol[ii4 + 3] * dNs_dz[ii];

    pr += rotated_local_sol[ii4 + 0] * Nr[ii];
    ur += rotated_local_sol[ii4 + 1] * Nr[ii];
    vr += rotated_local_sol[ii4 + 2] * Nr[ii];
    wr += rotated_local_sol[ii4 + 3] * Nr[ii];

    mur += rotated_local_mvelo[ii3 + 0] * Nr[ii];
    mvr += rotated_local_mvelo[ii3 + 1] * Nr[ii];
    mwr += rotated_local_mvelo[ii3 + 2] * Nr[ii];

    ur_x += rotated_local_sol[ii4 + 1] * dNr_dx[ii];
    ur_y += rotated_local_sol[ii4 + 1] * dNr_dy[ii];
    ur_z += rotated_local_sol[ii4 + 1] * dNr_dz[ii];

    vr_x += rotated_local_sol[ii4 + 2] * dNr_dx[ii];
    vr_y += rotated_local_sol[ii4 + 2] * dNr_dy[ii];
    vr_z += rotated_local_sol[ii4 + 2] * dNr_dz[ii];

    wr_x += rotated_local_sol[ii4 + 3] * dNr_dx[ii];
    wr_y += rotated_local_sol[ii4 + 3] * dNr_dy[ii];
    wr_z += rotated_local_sol[ii4 + 3] * dNr_dz[ii];
  }

  // Mesh velocity in the quadrature point
  const Vector_3 velo_mesh = Vector_3(mur, mvr, mwr);

  const Vector_3 velo_jump(us - ur, vs - vr, ws - wr);
  const double nsx {normal_s.x()}, nsy {normal_s.y()}, nsz {normal_s.z()};
  const double nrx {normal_r.x()}, nry {normal_r.y()}, nrz {normal_r.z()};

  double inflow_s = (us * nsx + vs * nsy + ws * nsz < 0.0 ? us * nsx + vs * nsy + ws * nsz : 0.0);

  const double gwts = fixed_J * fixed_qw;

  for(int A{0}; A<nLocBas; ++A)
  {
    const double NAs {Ns[A]}, NAs_x {dNs_dx[A]}, NAs_y {dNs_dy[A]}, NAs_z {dNs_dz[A]};

    const int A4 = 4 * A;

    // Mass -- s
    Residual_s[A4] += -0.5 * gwts * NAs * velo_jump.dot_product(normal_s);

    // x-dir -- s
    Residual_s[A4+1] += gwts * (0.5 * ( NAs * (ps * nsx - vis_mu * (2 * us_x * nsx + (us_y + vs_x) * nsy + (us_z + ws_x) * nsz)
      - pr * nrx + vis_mu * (2 * ur_x * nrx + (ur_y + vr_x) * nry + (ur_z + wr_x) * nrz))
      - (2 * NAs_x * nsx + NAs_y * nsy + NAs_z * nsz) * vis_mu * velo_jump.x()
      - NAs_y * nsx * vis_mu * velo_jump.y()
      - NAs_z * nsx * vis_mu * velo_jump.z() )
      + NAs * (tau_I - rho0 * inflow_s) * velo_jump.x());

    // y-dir -- s
    Residual_s[A4+2] += gwts * (0.5 * ( NAs * (ps * nsy - vis_mu * ((us_y + vs_x) * nsx + 2 * vs_y * nsy + (vs_z + ws_y) * nsz)
      - pr * nry + vis_mu * ((ur_y + vr_x) * nrx + 2 * vr_y * nry + (vr_z + wr_y) * nrz))
      - NAs_x * nsy * vis_mu * velo_jump.x()
      - (NAs_x * nsx + 2 * NAs_y * nsy + NAs_z * nsz) * vis_mu * velo_jump.y()
      - NAs_z * nsy * vis_mu * velo_jump.z() )
      + NAs * (tau_I - rho0 * inflow_s) * velo_jump.y());

    // z-dir -- s
    Residual_s[A4+3] += gwts * (0.5 * ( NAs * (ps * nsz - vis_mu * ((us_z + ws_x) * nsx + (vs_z + ws_y) * nsy + 2 * ws_z * nsz)
      - pr * nrz + vis_mu * ((ur_z + wr_x) * nrx + (vr_z + wr_y) * nry + 2 * wr_z * nrz))
      - NAs_x * nsz * vis_mu * velo_jump.x()
      - NAs_y * nsz * vis_mu * velo_jump.y()
      - (NAs_x * nsx + NAs_y * nsy + 2 * NAs_z * nsz) * vis_mu * velo_jump.z() )
      + NAs * (tau_I - rho0 * inflow_s) * velo_jump.z());
  
  }
}

void PLocAssem_VMS_NS_GenAlpha_Interface::Assem_Residual_itf_rotated(
  const int &qua, const double &rotated_qw, const double &dt,
  const FEAElement * const &rotated_elementv, const FEAElement * const &fixed_elementv,
  const double * const &rotated_local_sol, const double * const &rotated_local_mvelo,
  const double * const &fixed_local_sol)
{
  Zero_Residual_r();
  double ps {0.0};
  double us {0.0}, us_x {0.0}, us_y {0.0}, us_z {0.0};
  double vs {0.0}, vs_x {0.0}, vs_y {0.0}, vs_z {0.0};
  double ws {0.0}, ws_x {0.0}, ws_y {0.0}, ws_z {0.0};

  double pr {0.0};
  double ur {0.0}, ur_x {0.0}, ur_y {0.0}, ur_z {0.0};
  double vr {0.0}, vr_x {0.0}, vr_y {0.0}, vr_z {0.0};
  double wr {0.0}, wr_x {0.0}, wr_y {0.0}, wr_z {0.0};
  double mur {0.0}, mvr {0.0}, mwr{0.0};
  
  std::vector<double> Ns(nLocBas, 0.0), dNs_dx(nLocBas, 0.0), dNs_dy(nLocBas, 0.0), dNs_dz(nLocBas, 0.0);
  std::vector<double> Nr(nLocBas, 0.0), dNr_dx(nLocBas, 0.0), dNr_dy(nLocBas, 0.0), dNr_dz(nLocBas, 0.0);

  fixed_elementv -> get_R_gradR( 0, &Ns[0], &dNs_dx[0], &dNs_dy[0], &dNs_dz[0] );
  rotated_elementv -> get_R_gradR( qua, &Nr[0], &dNr_dx[0], &dNr_dy[0], &dNr_dz[0] );

  double  rotated_J {0.0};
  const Vector_3 normal_r = rotated_elementv -> get_2d_normal_out(qua, rotated_J);
  const Vector_3 normal_s = -1 * normal_r;

  // Calculate h_b and tau_I
  const auto s_dxi_dx = fixed_elementv -> get_invJacobian(0);
  const double h_s = get_h_b(s_dxi_dx, normal_s);

  const auto r_dxi_dx = rotated_elementv -> get_invJacobian(qua);
  const double h_r = get_h_b(r_dxi_dx, normal_r);

  const double tau_I = 0.5 * vis_mu * C_bI * (1.0 / h_s + 1.0 / h_r);

  for(int ii{0}; ii<nLocBas; ++ii)
  {
    const int ii4{4 * ii};
    const int ii3{3 * ii};

    pr += rotated_local_sol[ii4 + 0] * Nr[ii];
    ur += rotated_local_sol[ii4 + 1] * Nr[ii];
    vr += rotated_local_sol[ii4 + 2] * Nr[ii];
    wr += rotated_local_sol[ii4 + 3] * Nr[ii];

    ur_x += rotated_local_sol[ii4 + 1] * dNr_dx[ii];
    ur_y += rotated_local_sol[ii4 + 1] * dNr_dy[ii];
    ur_z += rotated_local_sol[ii4 + 1] * dNr_dz[ii];

    vr_x += rotated_local_sol[ii4 + 2] * dNr_dx[ii];
    vr_y += rotated_local_sol[ii4 + 2] * dNr_dy[ii];
    vr_z += rotated_local_sol[ii4 + 2] * dNr_dz[ii];

    wr_x += rotated_local_sol[ii4 + 3] * dNr_dx[ii];
    wr_y += rotated_local_sol[ii4 + 3] * dNr_dy[ii];
    wr_z += rotated_local_sol[ii4 + 3] * dNr_dz[ii];

    mur += rotated_local_mvelo[ii3 + 0] * Nr[ii];
    mvr += rotated_local_mvelo[ii3 + 1] * Nr[ii];
    mwr += rotated_local_mvelo[ii3 + 2] * Nr[ii];

    ps += fixed_local_sol[ii4 + 0] * Ns[ii];
    us += fixed_local_sol[ii4 + 1] * Ns[ii];
    vs += fixed_local_sol[ii4 + 2] * Ns[ii];
    ws += fixed_local_sol[ii4 + 3] * Ns[ii];

    us_x += fixed_local_sol[ii4 + 1] * dNs_dx[ii];
    us_y += fixed_local_sol[ii4 + 1] * dNs_dy[ii];
    us_z += fixed_local_sol[ii4 + 1] * dNs_dz[ii];

    vs_x += fixed_local_sol[ii4 + 2] * dNs_dx[ii];
    vs_y += fixed_local_sol[ii4 + 2] * dNs_dy[ii];
    vs_z += fixed_local_sol[ii4 + 2] * dNs_dz[ii];

    ws_x += fixed_local_sol[ii4 + 3] * dNs_dx[ii];
    ws_y += fixed_local_sol[ii4 + 3] * dNs_dy[ii];
    ws_z += fixed_local_sol[ii4 + 3] * dNs_dz[ii];
  }

  // Mesh velocity in the quadrature point
  const Vector_3 velo_mesh = Vector_3(mur, mvr, mwr);

  const Vector_3 velo_jump(us - ur, vs - vr, ws - wr);
  const double nsx {normal_s.x()}, nsy {normal_s.y()}, nsz {normal_s.z()};
  const double nrx {normal_r.x()}, nry {normal_r.y()}, nrz {normal_r.z()};

  double inflow_r = ((ur - velo_mesh.x()) * nrx + (vr - velo_mesh.y()) * nry + (wr - velo_mesh.z()) * nrz ?
                     (ur - velo_mesh.x()) * nrx + (vr - velo_mesh.y()) * nry + (wr - velo_mesh.z()) * nrz : 0.0);

  const double gwts = rotated_J * rotated_qw;

  for(int A{0}; A<nLocBas; ++A)
  {
    const double NAr {Nr[A]}, NAr_x {dNr_dx[A]}, NAr_y {dNr_dy[A]}, NAr_z {dNr_dz[A]};

    const int A4 = 4 * A;

    // Mass -- r
    Residual_r[A4] += 0.5 * gwts * NAr * velo_jump.dot_product(normal_r);

    // x-dir -- r
    Residual_r[A4+1] += gwts * (0.5 * ( NAr * (-1.0 * ps * nsx + vis_mu * (2 * us_x * nsx + (us_y + vs_x) * nsy + (us_z + ws_x) * nsz)
      + pr * nrx - vis_mu * (2 * ur_x * nrx + (ur_y + vr_x) * nry + (ur_z + wr_x) * nrz))
      + (2 * NAr_x * nrx + NAr_y * nry + NAr_z * nrz) * vis_mu * velo_jump.x()
      + NAr_y * nrx * vis_mu * velo_jump.y()
      + NAr_z * nrx * vis_mu * velo_jump.z() )
      + NAr * (rho0 * inflow_r - tau_I) * velo_jump.x());

    // y-dir -- r
    Residual_r[A4+2] += gwts * (0.5 * ( NAr * (-1.0 * ps * nsy + vis_mu * ((us_y + vs_x) * nsx + 2 * vs_y * nsy + (vs_z + ws_y) * nsz)
      + pr * nry - vis_mu * ((ur_y + vr_x) * nrx + 2 * vr_y * nry + (vr_z + wr_y) * nrz))
      + NAr_x * nry * vis_mu * velo_jump.x()
      + (NAr_x * nrx + 2 * NAr_y * nry + NAr_z * nrz) * vis_mu * velo_jump.y()
      + NAr_z * nry * vis_mu * velo_jump.z() )
      + NAr * (rho0 * inflow_r - tau_I) * velo_jump.y());

    // z-dir -- r
    Residual_r[A4+3] += gwts * (0.5 * ( NAr * (-1.0 * ps * nsz + vis_mu * ((us_z + ws_x) * nsx + (vs_z + ws_y) * nsy + 2 * ws_z * nsz)
      + pr * nrz - vis_mu * ((ur_z + wr_x) * nrx + (vr_z + wr_y) * nry + 2 * wr_z * nrz))
      + NAr_x * nrz * vis_mu * velo_jump.x()
      + NAr_y * nrz * vis_mu * velo_jump.y()
      + (NAr_x * nrx + NAr_y * nry + 2 * NAr_z * nrz) * vis_mu * velo_jump.z() )
      + NAr * (rho0 * inflow_r - tau_I) * velo_jump.z());
  }
}

void PLocAssem_VMS_NS_GenAlpha_Interface::Assem_Diag_Tangent_Residual_itf_fixed(
  const int &qua, const double &fixed_qw, const double &dt,
  const FEAElement * const &fixed_elementv, const FEAElement * const &rotated_elementv,
  const double * const &fixed_local_sol, const double * const &rotated_local_sol,
  const double * const &rotated_local_mvelo)
{
  Zero_Residual_s();
  Zero_Tangent_ss();
  double ps {0.0};
  double us {0.0}, us_x {0.0}, us_y {0.0}, us_z {0.0};
  double vs {0.0}, vs_x {0.0}, vs_y {0.0}, vs_z {0.0};
  double ws {0.0}, ws_x {0.0}, ws_y {0.0}, ws_z {0.0};

  double pr {0.0};
  double ur {0.0}, ur_x {0.0}, ur_y {0.0}, ur_z {0.0};
  double vr {0.0}, vr_x {0.0}, vr_y {0.0}, vr_z {0.0};
  double wr {0.0}, wr_x {0.0}, wr_y {0.0}, wr_z {0.0};
  double mur {0.0}, mvr {0.0}, mwr{0.0};

  std::vector<double> Ns(nLocBas, 0.0), dNs_dx(nLocBas, 0.0), dNs_dy(nLocBas, 0.0), dNs_dz(nLocBas, 0.0);
  std::vector<double> Nr(nLocBas, 0.0), dNr_dx(nLocBas, 0.0), dNr_dy(nLocBas, 0.0), dNr_dz(nLocBas, 0.0);

  fixed_elementv -> get_R_gradR( qua, &Ns[0], &dNs_dx[0], &dNs_dy[0], &dNs_dz[0] );
  rotated_elementv -> get_R_gradR( 0, &Nr[0], &dNr_dx[0], &dNr_dy[0], &dNr_dz[0] );

  double fixed_J {0.0};
  const Vector_3 normal_s = fixed_elementv -> get_2d_normal_out(qua, fixed_J);
  const Vector_3 normal_r = -1 * normal_s;

  // Calculate h_b and tau_I
  const auto s_dxi_dx = fixed_elementv -> get_invJacobian(qua);
  const double h_s = get_h_b(s_dxi_dx, normal_s);

  const auto r_dxi_dx = rotated_elementv -> get_invJacobian(0);
  const double h_r = get_h_b(r_dxi_dx, normal_r);

  const double tau_I = 0.5 * vis_mu * C_bI * (1.0 / h_s + 1.0 / h_r);

  for(int ii{0}; ii<nLocBas; ++ii)
  {
    const int ii4{4 * ii};
    const int ii3{3 * ii};

    ps += fixed_local_sol[ii4 + 0] * Ns[ii];
    us += fixed_local_sol[ii4 + 1] * Ns[ii];
    vs += fixed_local_sol[ii4 + 2] * Ns[ii];
    ws += fixed_local_sol[ii4 + 3] * Ns[ii];

    us_x += fixed_local_sol[ii4 + 1] * dNs_dx[ii];
    us_y += fixed_local_sol[ii4 + 1] * dNs_dy[ii];
    us_z += fixed_local_sol[ii4 + 1] * dNs_dz[ii];

    vs_x += fixed_local_sol[ii4 + 2] * dNs_dx[ii];
    vs_y += fixed_local_sol[ii4 + 2] * dNs_dy[ii];
    vs_z += fixed_local_sol[ii4 + 2] * dNs_dz[ii];

    ws_x += fixed_local_sol[ii4 + 3] * dNs_dx[ii];
    ws_y += fixed_local_sol[ii4 + 3] * dNs_dy[ii];
    ws_z += fixed_local_sol[ii4 + 3] * dNs_dz[ii];

    pr += rotated_local_sol[ii4 + 0] * Nr[ii];
    ur += rotated_local_sol[ii4 + 1] * Nr[ii];
    vr += rotated_local_sol[ii4 + 2] * Nr[ii];
    wr += rotated_local_sol[ii4 + 3] * Nr[ii];

    mur += rotated_local_mvelo[ii3 + 0] * Nr[ii];
    mvr += rotated_local_mvelo[ii3 + 1] * Nr[ii];
    mwr += rotated_local_mvelo[ii3 + 2] * Nr[ii];

    ur_x += rotated_local_sol[ii4 + 1] * dNr_dx[ii];
    ur_y += rotated_local_sol[ii4 + 1] * dNr_dy[ii];
    ur_z += rotated_local_sol[ii4 + 1] * dNr_dz[ii];

    vr_x += rotated_local_sol[ii4 + 2] * dNr_dx[ii];
    vr_y += rotated_local_sol[ii4 + 2] * dNr_dy[ii];
    vr_z += rotated_local_sol[ii4 + 2] * dNr_dz[ii];

    wr_x += rotated_local_sol[ii4 + 3] * dNr_dx[ii];
    wr_y += rotated_local_sol[ii4 + 3] * dNr_dy[ii];
    wr_z += rotated_local_sol[ii4 + 3] * dNr_dz[ii];
  }

  // Mesh velocity in the quadrature point
  const Vector_3 velo_mesh = Vector_3(mur, mvr, mwr);

  const Vector_3 velo_jump(us - ur, vs - vr, ws - wr);
  const double nsx {normal_s.x()}, nsy {normal_s.y()}, nsz {normal_s.z()};
  const double nrx {normal_r.x()}, nry {normal_r.y()}, nrz {normal_r.z()};

  double inflow_s = (us * nsx + vs * nsy + ws * nsz < 0.0 ? us * nsx + vs * nsy + ws * nsz : 0.0);
  double delta_s = (us * nsx + vs * nsy + ws * nsz < 0.0 ? 1.0 : 0.0);

  const double gwts = fixed_J * fixed_qw;
  const double dd_dv = alpha_f * gamma * dt;
  const double common_coef = gwts * dd_dv;

  for(int A{0}; A<nLocBas; ++A)
  {
    const double NAs {Ns[A]}, NAs_x {dNs_dx[A]}, NAs_y {dNs_dy[A]}, NAs_z {dNs_dz[A]};

    const int A4 = 4 * A;

    // Mass -- s
    Residual_s[A4] += -0.5 * gwts * NAs * velo_jump.dot_product(normal_s);

    // x-dir -- s
    Residual_s[A4+1] += gwts * (0.5 * ( NAs * (ps * nsx - vis_mu * (2 * us_x * nsx + (us_y + vs_x) * nsy + (us_z + ws_x) * nsz)
      - pr * nrx + vis_mu * (2 * ur_x * nrx + (ur_y + vr_x) * nry + (ur_z + wr_x) * nrz))
      - (2 * NAs_x * nsx + NAs_y * nsy + NAs_z * nsz) * vis_mu * velo_jump.x()
      - NAs_y * nsx * vis_mu * velo_jump.y()
      - NAs_z * nsx * vis_mu * velo_jump.z() )
      + NAs * (tau_I - rho0 * inflow_s) * velo_jump.x());

    // y-dir -- s
    Residual_s[A4+2] += gwts * (0.5 * ( NAs * (ps * nsy - vis_mu * ((us_y + vs_x) * nsx + 2 * vs_y * nsy + (vs_z + ws_y) * nsz)
      - pr * nry + vis_mu * ((ur_y + vr_x) * nrx + 2 * vr_y * nry + (vr_z + wr_y) * nrz))
      - NAs_x * nsy * vis_mu * velo_jump.x()
      - (NAs_x * nsx + 2 * NAs_y * nsy + NAs_z * nsz) * vis_mu * velo_jump.y()
      - NAs_z * nsy * vis_mu * velo_jump.z() )
      + NAs * (tau_I - rho0 * inflow_s) * velo_jump.y());

    // z-dir -- s
    Residual_s[A4+3] += gwts * (0.5 * ( NAs * (ps * nsz - vis_mu * ((us_z + ws_x) * nsx + (vs_z + ws_y) * nsy + 2 * ws_z * nsz)
      - pr * nrz + vis_mu * ((ur_z + wr_x) * nrx + (vr_z + wr_y) * nry + 2 * wr_z * nrz))
      - NAs_x * nsz * vis_mu * velo_jump.x()
      - NAs_y * nsz * vis_mu * velo_jump.y()
      - (NAs_x * nsx + NAs_y * nsy + 2 * NAs_z * nsz) * vis_mu * velo_jump.z() )
      + NAs * (tau_I - rho0 * inflow_s) * velo_jump.z());

    for(int B{0}; B<nLocBas; ++B)
    {
      const double NBs {Ns[B]}, NBs_x {dNs_dx[B]}, NBs_y {dNs_dy[B]}, NBs_z {dNs_dz[B]};

      const int B4 = 4 * B, nLocBas4 = 4 * nLocBas;

      const double NAsNBs = NAs * NBs; 

      // Ks1s1
      // Tangent_ss[nLocBas4 * A4 + B4] += 0;

      // Ks1s2
      Tangent_ss[nLocBas4 * A4 + B4 + 1] += -0.5 * common_coef * NAsNBs * nsx;

      // Ks1s3
      Tangent_ss[nLocBas4 * A4 + B4 + 2] += -0.5 * common_coef * NAsNBs * nsy;

      // Ks1s4
      Tangent_ss[nLocBas4 * A4 + B4 + 3] += -0.5 * common_coef * NAsNBs * nsz;

      // Ks2s1
      Tangent_ss[nLocBas4 * (A4 + 1) + B4] += 0.5 * common_coef * NAsNBs * nsx;

      // Ks2s2
      Tangent_ss[nLocBas4 * (A4 + 1) + B4 + 1] += common_coef * (
        -0.5 * vis_mu * (NAs * (2 * NBs_x * nsx + NBs_y * nsy + NBs_z * nsz)
                      + (2 * NAs_x * nsx + NAs_y * nsy + NAs_z * nsz) * NBs)
        + NAsNBs * (tau_I - rho0 * (delta_s * velo_jump.x() * nsx + inflow_s)) );

      // Ks2s3
      Tangent_ss[nLocBas4 * (A4 + 1) + B4 + 2] += common_coef * (
        -0.5 * vis_mu * (NAs * NBs_x * nsy + NAs_y * nsx * NBs) - NAsNBs * rho0 * delta_s * velo_jump.x() * nsy );

      // Ks2s4
      Tangent_ss[nLocBas4 * (A4 + 1) + B4 + 3] += common_coef * (
        -0.5 * vis_mu * (NAs * NBs_x * nsz + NAs_z * nsx * NBs) - NAsNBs * rho0 * delta_s * velo_jump.x() * nsz );

      // Ks3s1
      Tangent_ss[nLocBas4 * (A4 + 2) + B4] += 0.5 * common_coef * NAsNBs * nsy;

      // Ks3s2
      Tangent_ss[nLocBas4 * (A4 + 2) + B4 + 1] += common_coef * (
        -0.5 * vis_mu * (NAs * NBs_y * nsx + NAs_x * nsy * NBs) - NAsNBs * rho0 * delta_s * velo_jump.y() * nsx );

      // Ks3s3
      Tangent_ss[nLocBas4 * (A4 + 2) + B4 + 2] += common_coef * (
        -0.5 * vis_mu * (NAs * (NBs_x * nsx + 2 * NBs_y * nsy + NBs_z * nsz)
                      + (NAs_x * nsx + 2 * NAs_y * nsy + NAs_z * nsz) * NBs)
        + NAsNBs * (tau_I - rho0 * (delta_s * velo_jump.y() * nsy + inflow_s)) );

      // Ks3s4
      Tangent_ss[nLocBas4 * (A4 + 2) + B4 + 3] += common_coef * (
        -0.5 * vis_mu * (NAs * NBs_y * nsz + NAs_z * nsy * NBs) - NAsNBs * rho0 * delta_s * velo_jump.y() * nsz );

      // Ks4s1
      Tangent_ss[nLocBas4 * (A4 + 3) + B4] += 0.5 * common_coef * NAsNBs * nsz;

      // Ks4s2
      Tangent_ss[nLocBas4 * (A4 + 3) + B4 + 1] += common_coef * (
        -0.5 * vis_mu * (NAs * NBs_z * nsx + NAs_x * nsz * NBs) - NAsNBs * rho0 * delta_s * velo_jump.z() * nsx );

      // Ks4s3
      Tangent_ss[nLocBas4 * (A4 + 3) + B4 + 2] += common_coef * (
        -0.5 * vis_mu * (NAs * NBs_z * nsy + NAs_y * nsz * NBs) - NAsNBs * rho0 * delta_s * velo_jump.z() * nsy );

      // Ks4s4
      Tangent_ss[nLocBas4 * (A4 + 3) + B4 + 3] += common_coef * (
        -0.5 * vis_mu * (NAs * (NBs_x * nsx + NBs_y * nsy + 2 * NBs_z * nsz)
                      + (NAs_x * nsx + NAs_y * nsy + 2 * NAs_z * nsz) * NBs)
        + NAsNBs * (tau_I - rho0 * (delta_s * velo_jump.z() * nsz + inflow_s)) );
    }
  }
}

void PLocAssem_VMS_NS_GenAlpha_Interface::Assem_Diag_Tangent_Residual_itf_rotated(
  const int &qua, const double &rotated_qw, const double &dt,
  const FEAElement * const &rotated_elementv, const FEAElement * const &fixed_elementv,
  const double * const &rotated_local_sol, const double * const &rotated_local_mvelo,
  const double * const &fixed_local_sol)
{
  Zero_Residual_r();
  Zero_Tangent_rr();
  double ps {0.0};
  double us {0.0}, us_x {0.0}, us_y {0.0}, us_z {0.0};
  double vs {0.0}, vs_x {0.0}, vs_y {0.0}, vs_z {0.0};
  double ws {0.0}, ws_x {0.0}, ws_y {0.0}, ws_z {0.0};

  double pr {0.0};
  double ur {0.0}, ur_x {0.0}, ur_y {0.0}, ur_z {0.0};
  double vr {0.0}, vr_x {0.0}, vr_y {0.0}, vr_z {0.0};
  double wr {0.0}, wr_x {0.0}, wr_y {0.0}, wr_z {0.0};
  double mur {0.0}, mvr {0.0}, mwr{0.0};
  
  std::vector<double> Ns(nLocBas, 0.0), dNs_dx(nLocBas, 0.0), dNs_dy(nLocBas, 0.0), dNs_dz(nLocBas, 0.0);
  std::vector<double> Nr(nLocBas, 0.0), dNr_dx(nLocBas, 0.0), dNr_dy(nLocBas, 0.0), dNr_dz(nLocBas, 0.0);

  fixed_elementv -> get_R_gradR( 0, &Ns[0], &dNs_dx[0], &dNs_dy[0], &dNs_dz[0] );
  rotated_elementv -> get_R_gradR( qua, &Nr[0], &dNr_dx[0], &dNr_dy[0], &dNr_dz[0] );

  double  rotated_J {0.0};
  const Vector_3 normal_r = rotated_elementv -> get_2d_normal_out(qua, rotated_J);
  const Vector_3 normal_s = -1 * normal_r;

  // Calculate h_b and tau_I
  const auto s_dxi_dx = fixed_elementv -> get_invJacobian(0);
  const double h_s = get_h_b(s_dxi_dx, normal_s);

  const auto r_dxi_dx = rotated_elementv -> get_invJacobian(qua);
  const double h_r = get_h_b(r_dxi_dx, normal_r);

  const double tau_I = 0.5 * vis_mu * C_bI * (1.0 / h_s + 1.0 / h_r);

  for(int ii{0}; ii<nLocBas; ++ii)
  {
    const int ii4{4 * ii};
    const int ii3{3 * ii};

    pr += rotated_local_sol[ii4 + 0] * Nr[ii];
    ur += rotated_local_sol[ii4 + 1] * Nr[ii];
    vr += rotated_local_sol[ii4 + 2] * Nr[ii];
    wr += rotated_local_sol[ii4 + 3] * Nr[ii];

    ur_x += rotated_local_sol[ii4 + 1] * dNr_dx[ii];
    ur_y += rotated_local_sol[ii4 + 1] * dNr_dy[ii];
    ur_z += rotated_local_sol[ii4 + 1] * dNr_dz[ii];

    vr_x += rotated_local_sol[ii4 + 2] * dNr_dx[ii];
    vr_y += rotated_local_sol[ii4 + 2] * dNr_dy[ii];
    vr_z += rotated_local_sol[ii4 + 2] * dNr_dz[ii];

    wr_x += rotated_local_sol[ii4 + 3] * dNr_dx[ii];
    wr_y += rotated_local_sol[ii4 + 3] * dNr_dy[ii];
    wr_z += rotated_local_sol[ii4 + 3] * dNr_dz[ii];

    mur += rotated_local_mvelo[ii3 + 0] * Nr[ii];
    mvr += rotated_local_mvelo[ii3 + 1] * Nr[ii];
    mwr += rotated_local_mvelo[ii3 + 2] * Nr[ii];

    ps += fixed_local_sol[ii4 + 0] * Ns[ii];
    us += fixed_local_sol[ii4 + 1] * Ns[ii];
    vs += fixed_local_sol[ii4 + 2] * Ns[ii];
    ws += fixed_local_sol[ii4 + 3] * Ns[ii];

    us_x += fixed_local_sol[ii4 + 1] * dNs_dx[ii];
    us_y += fixed_local_sol[ii4 + 1] * dNs_dy[ii];
    us_z += fixed_local_sol[ii4 + 1] * dNs_dz[ii];

    vs_x += fixed_local_sol[ii4 + 2] * dNs_dx[ii];
    vs_y += fixed_local_sol[ii4 + 2] * dNs_dy[ii];
    vs_z += fixed_local_sol[ii4 + 2] * dNs_dz[ii];

    ws_x += fixed_local_sol[ii4 + 3] * dNs_dx[ii];
    ws_y += fixed_local_sol[ii4 + 3] * dNs_dy[ii];
    ws_z += fixed_local_sol[ii4 + 3] * dNs_dz[ii];
  }

  // Mesh velocity in the quadrature point
  const Vector_3 velo_mesh = Vector_3(mur, mvr, mwr);

  const Vector_3 velo_jump(us - ur, vs - vr, ws - wr);
  const double nsx {normal_s.x()}, nsy {normal_s.y()}, nsz {normal_s.z()};
  const double nrx {normal_r.x()}, nry {normal_r.y()}, nrz {normal_r.z()};

  double inflow_r = ((ur - velo_mesh.x()) * nrx + (vr - velo_mesh.y()) * nry + (wr - velo_mesh.z()) * nrz ?
                     (ur - velo_mesh.x()) * nrx + (vr - velo_mesh.y()) * nry + (wr - velo_mesh.z()) * nrz : 0.0);
  double delta_r = ((ur - velo_mesh.x()) * nrx + (vr - velo_mesh.y()) * nry + (wr - velo_mesh.z()) * nrz ?
                     1.0 : 0.0);

  const double gwts = rotated_J * rotated_qw;
  const double dd_dv = alpha_f * gamma * dt;
  const double common_coef = gwts * dd_dv;

  for(int A{0}; A<nLocBas; ++A)
  {
    const double NAr {Nr[A]}, NAr_x {dNr_dx[A]}, NAr_y {dNr_dy[A]}, NAr_z {dNr_dz[A]};

    const int A4 = 4 * A;

    // Mass -- r
    Residual_r[A4] += 0.5 * gwts * NAr * velo_jump.dot_product(normal_r);

    // x-dir -- r
    Residual_r[A4+1] += gwts * (0.5 * ( NAr * (-1.0 * ps * nsx + vis_mu * (2 * us_x * nsx + (us_y + vs_x) * nsy + (us_z + ws_x) * nsz)
      + pr * nrx - vis_mu * (2 * ur_x * nrx + (ur_y + vr_x) * nry + (ur_z + wr_x) * nrz))
      + (2 * NAr_x * nrx + NAr_y * nry + NAr_z * nrz) * vis_mu * velo_jump.x()
      + NAr_y * nrx * vis_mu * velo_jump.y()
      + NAr_z * nrx * vis_mu * velo_jump.z() )
      + NAr * (rho0 * inflow_r - tau_I) * velo_jump.x());

    // y-dir -- r
    Residual_r[A4+2] += gwts * (0.5 * ( NAr * (-1.0 * ps * nsy + vis_mu * ((us_y + vs_x) * nsx + 2 * vs_y * nsy + (vs_z + ws_y) * nsz)
      + pr * nry - vis_mu * ((ur_y + vr_x) * nrx + 2 * vr_y * nry + (vr_z + wr_y) * nrz))
      + NAr_x * nry * vis_mu * velo_jump.x()
      + (NAr_x * nrx + 2 * NAr_y * nry + NAr_z * nrz) * vis_mu * velo_jump.y()
      + NAr_z * nry * vis_mu * velo_jump.z() )
      + NAr * (rho0 * inflow_r - tau_I) * velo_jump.y());

    // z-dir -- r
    Residual_r[A4+3] += gwts * (0.5 * ( NAr * (-1.0 * ps * nsz + vis_mu * ((us_z + ws_x) * nsx + (vs_z + ws_y) * nsy + 2 * ws_z * nsz)
      + pr * nrz - vis_mu * ((ur_z + wr_x) * nrx + (vr_z + wr_y) * nry + 2 * wr_z * nrz))
      + NAr_x * nrz * vis_mu * velo_jump.x()
      + NAr_y * nrz * vis_mu * velo_jump.y()
      + (NAr_x * nrx + NAr_y * nry + 2 * NAr_z * nrz) * vis_mu * velo_jump.z() )
      + NAr * (rho0 * inflow_r - tau_I) * velo_jump.z());

    for(int B{0}; B<nLocBas; ++B)
    {
      const double NBr {Nr[B]}, NBr_x {dNr_dx[B]}, NBr_y {dNr_dy[B]}, NBr_z {dNr_dz[B]};

      const int B4 = 4 * B, nLocBas4 = 4 * nLocBas;

      const double NArNBr = NAr * NBr; 

      // Kr1r1
      // Tangent_rr[nLocBas4 * A4 + B4] += 0;

      // Kr1r2
      Tangent_rr[nLocBas4 * A4 + B4 + 1] += -0.5 * common_coef * NArNBr * nrx;

      // Kr1r3
      Tangent_rr[nLocBas4 * A4 + B4 + 2] += -0.5 * common_coef * NArNBr * nry;

      // Kr1r4
      Tangent_rr[nLocBas4 * A4 + B4 + 3] += -0.5 * common_coef * NArNBr * nrz;

      // Kr2r1
      Tangent_rr[nLocBas4 * (A4 + 1) + B4] += 0.5 * common_coef * NArNBr * nrx;

      // Kr2r2
      Tangent_rr[nLocBas4 * (A4 + 1) + B4 + 1] += common_coef * (
        -0.5 * vis_mu * (NAr * (2 * NBr_x * nrx + NBr_y * nry + NBr_z * nrz)
                      + (2 * NAr_x * nrx + NAr_y * nry + NAr_z * nrz) * NBr)
        + NArNBr * (tau_I + rho0 * (delta_r * velo_jump.x() * nrx - inflow_r)) );
      // velo_jump is Us - Ur, so here is 'tau_I + ...'

      // Kr2r3
      Tangent_rr[nLocBas4 * (A4 + 1) + B4 + 2] += common_coef * (
        -0.5 * vis_mu * (NAr * NBr_x * nry + NAr_y * nrx * NBr) + NArNBr * rho0 * delta_r * velo_jump.x() * nry );

      // Kr2r4
      Tangent_rr[nLocBas4 * (A4 + 1) + B4 + 3] += common_coef * (
        -0.5 * vis_mu * (NAr * NBr_x * nrz + NAr_z * nrx * NBr) + NArNBr * rho0 * delta_r * velo_jump.x() * nrz );

      // Kr3r1
      Tangent_rr[nLocBas4 * (A4 + 2) + B4] += 0.5 * common_coef * NArNBr * nry;

      // Kr3r2
      Tangent_rr[nLocBas4 * (A4 + 2) + B4 + 1] += common_coef * (
        -0.5 * vis_mu * (NAr * NBr_y * nrx + NAr_x * nry * NBr) + NArNBr * rho0 * delta_r * velo_jump.y() * nrx );

      // Kr3r3
      Tangent_rr[nLocBas4 * (A4 + 2) + B4 + 2] += common_coef * (
        -0.5 * vis_mu * (NAr * (NBr_x * nrx + 2 * NBr_y * nry + NBr_z * nrz)
                      + (NAr_x * nrx + 2 * NAr_y * nry + NAr_z * nrz) * NBr)
        + NArNBr * (tau_I + rho0 * (delta_r * velo_jump.y() * nry - inflow_r)) );

      // Kr3r4
      Tangent_rr[nLocBas4 * (A4 + 2) + B4 + 3] += common_coef * (
        -0.5 * vis_mu * (NAr * NBr_y * nrz + NAr_z * nry * NBr) + NArNBr * rho0 * delta_r * velo_jump.y() * nrz );

      // Kr4r1
      Tangent_rr[nLocBas4 * (A4 + 3) + B4] += 0.5 * common_coef * NArNBr * nrz;

      // Kr4r2
      Tangent_rr[nLocBas4 * (A4 + 3) + B4 + 1] += common_coef * (
        -0.5 * vis_mu * (NAr * NBr_z * nrx + NAr_x * nrz * NBr) + NArNBr * rho0 * delta_r * velo_jump.z() * nrx );

      // Kr4r3
      Tangent_rr[nLocBas4 * (A4 + 3) + B4 + 2] += common_coef * (
        -0.5 * vis_mu * (NAr * NBr_z * nry + NAr_y * nrz * NBr) + NArNBr * rho0 * delta_r * velo_jump.z() * nry );

      // Kr4r4
      Tangent_rr[nLocBas4 * (A4 + 3) + B4 + 3] += common_coef * (
        -0.5 * vis_mu * (NAr * (NBr_x * nrx + NBr_y * nry + 2 * NBr_z * nrz)
                      + (NAr_x * nrx + NAr_y * nry + 2 * NAr_z * nrz) * NBr)
        + NArNBr * (tau_I + rho0 * (delta_r * velo_jump.z() * nrz - inflow_r)) );
    }
  }
}

void PLocAssem_VMS_NS_GenAlpha_Interface::Assem_Tangent_itf_MF_fixed(
  const int &qua, const double &fixed_qw, const double &dt,
  const FEAElement * const &fixed_elementv, const FEAElement * const &rotated_elementv,
  const double * const &fixed_local_sol, const double * const &rotated_local_sol,
  const double * const &rotated_local_mvelo)
{
  Zero_Tangent_sr();
  double ps {0.0}, us {0.0}, vs {0.0}, ws {0.0};
  double pr {0.0}, ur {0.0}, vr {0.0}, wr {0.0};
  double mur {0.0}, mvr {0.0}, mwr{0.0};

  std::vector<double> Ns(nLocBas, 0.0), dNs_dx(nLocBas, 0.0), dNs_dy(nLocBas, 0.0), dNs_dz(nLocBas, 0.0);
  std::vector<double> Nr(nLocBas, 0.0), dNr_dx(nLocBas, 0.0), dNr_dy(nLocBas, 0.0), dNr_dz(nLocBas, 0.0);

  fixed_elementv -> get_R_gradR( qua, &Ns[0], &dNs_dx[0], &dNs_dy[0], &dNs_dz[0] );
  rotated_elementv -> get_R_gradR( 0, &Nr[0], &dNr_dx[0], &dNr_dy[0], &dNr_dz[0] );

  double fixed_J {0.0};
  const Vector_3 normal_s = fixed_elementv -> get_2d_normal_out(qua, fixed_J);
  const Vector_3 normal_r = -1 * normal_s;

  // Calculate h_b and tau_I
  const auto s_dxi_dx = fixed_elementv -> get_invJacobian(qua);
  const double h_s = get_h_b(s_dxi_dx, normal_s);

  const auto r_dxi_dx = rotated_elementv -> get_invJacobian(0);
  const double h_r = get_h_b(r_dxi_dx, normal_r);

  const double tau_I = 0.5 * vis_mu * C_bI * (1.0 / h_s + 1.0 / h_r);

  for(int ii{0}; ii<nLocBas; ++ii)
  {
    const int ii4{4 * ii};
    const int ii3{3 * ii};

    ps += fixed_local_sol[ii4 + 0] * Ns[ii];
    us += fixed_local_sol[ii4 + 1] * Ns[ii];
    vs += fixed_local_sol[ii4 + 2] * Ns[ii];
    ws += fixed_local_sol[ii4 + 3] * Ns[ii];

    pr += rotated_local_sol[ii4 + 0] * Nr[ii];
    ur += rotated_local_sol[ii4 + 1] * Nr[ii];
    vr += rotated_local_sol[ii4 + 2] * Nr[ii];
    wr += rotated_local_sol[ii4 + 3] * Nr[ii];

    mur += rotated_local_mvelo[ii3 + 0] * Nr[ii];
    mvr += rotated_local_mvelo[ii3 + 1] * Nr[ii];
    mwr += rotated_local_mvelo[ii3 + 2] * Nr[ii];
  }

  // Mesh velocity in the quadrature point
  const Vector_3 velo_mesh = Vector_3(mur, mvr, mwr);

  const Vector_3 velo_jump(us - ur, vs - vr, ws - wr);
  const double nsx {normal_s.x()}, nsy {normal_s.y()}, nsz {normal_s.z()};
  const double nrx {normal_r.x()}, nry {normal_r.y()}, nrz {normal_r.z()};

  double inflow_s = (us * nsx + vs * nsy + ws * nsz < 0.0 ? us * nsx + vs * nsy + ws * nsz : 0.0);

  const double gwts = fixed_J * fixed_qw;
  const double dd_dv = alpha_f * gamma * dt;
  const double common_coef = gwts * dd_dv;

  for(int A{0}; A<nLocBas; ++A)
  {
    const double NAs {Ns[A]}, NAs_x {dNs_dx[A]}, NAs_y {dNs_dy[A]}, NAs_z {dNs_dz[A]};

    const int A4 = 4 * A;

    for(int B{0}; B<nLocBas; ++B)
    {
      const double NBr {Nr[B]}, NBr_x {dNr_dx[B]}, NBr_y {dNr_dy[B]}, NBr_z {dNr_dz[B]};

      const int B4 = 4 * B, nLocBas4 = 4 * nLocBas;

      const double NAsNBr = NAs * NBr;

      // Ks1r1
      // Tangent_sr[nLocBas4 * A4 + B4] += 0;

      // Ks1r2
      Tangent_sr[nLocBas4 * A4 + B4 + 1] += 0.5 * common_coef * NAsNBr * nsx;

      // Ks1r3
      Tangent_sr[nLocBas4 * A4 + B4 + 2] += 0.5 * common_coef * NAsNBr * nsy;

      // Ks1r4
      Tangent_sr[nLocBas4 * A4 + B4 + 3] += 0.5 * common_coef * NAsNBr * nsz;

      // Ks2r1
      Tangent_sr[nLocBas4 * (A4 + 1) + B4] += -0.5 * common_coef * NAsNBr * nrx;

      // Ks2r2
      Tangent_sr[nLocBas4 * (A4 + 1) + B4 + 1] += common_coef * (
        0.5 * vis_mu * (NAs * (2 * NBr_x * nrx + NBr_y * nry + NBr_z * nrz)
                     + (2 * NAs_x * nsx + NAs_y * nsy + NAs_z * nsz) * NBr)
        + NAsNBr * (rho0 * inflow_s - tau_I) );

      // Ks2r3
      Tangent_sr[nLocBas4 * (A4 + 1) + B4 + 2] += 0.5 * common_coef * vis_mu * (
        NAs * NBr_x * nry + NAs_y * nsx * NBr );

      // Ks2r4
      Tangent_sr[nLocBas4 * (A4 + 1) + B4 + 3] += 0.5 * common_coef * vis_mu * (
        NAs * NBr_x * nrz + NAs_z * nsx * NBr );

      // Ks3r1
      Tangent_sr[nLocBas4 * (A4 + 2) + B4] += -0.5 * common_coef * NAsNBr * nry;

      // Ks3r2
      Tangent_sr[nLocBas4 * (A4 + 2) + B4 + 1] += 0.5 * common_coef * vis_mu * (
        NAs * NBr_y * nrx + NAs_x * nsy * NBr );

      // Ks3r3
      Tangent_sr[nLocBas4 * (A4 + 2) + B4 + 2] += common_coef * (
        0.5 * vis_mu * (NAs * (NBr_x * nrx + 2 * NBr_y * nry + NBr_z * nrz)
                     + (NAs_x * nsx + 2 * NAs_y * nsy + NAs_z * nsz) * NBr)
        + NAsNBr * (rho0 * inflow_s - tau_I) );

      // Ks3r4
      Tangent_sr[nLocBas4 * (A4 + 2) + B4 + 3] += 0.5 * common_coef * vis_mu * (
        NAs * NBr_y * nrz + NAs_z * nsy * NBr );

      // Ks4r1
      Tangent_sr[nLocBas4 * (A4 + 3) + B4] += -0.5 * common_coef * NAsNBr * nrz;

      // Ks4r2
      Tangent_sr[nLocBas4 * (A4 + 3) + B4 + 1] += 0.5 * common_coef * vis_mu * (
        NAs * NBr_z * nrx + NAs_x * nsz * NBr );

      // Ks4r3
      Tangent_sr[nLocBas4 * (A4 + 3) + B4 + 2] += 0.5 * common_coef * vis_mu * (
        NAs * NBr_z * nry + NAs_y * nsz * NBr );

      // Ks4r4
      Tangent_sr[nLocBas4 * (A4 + 3) + B4 + 3] += common_coef * (
        0.5 * vis_mu * (NAs * (NBr_x * nrx + NBr_y * nry + 2 * NBr_z * nrz)
                     + (NAs_x * nsx + NAs_y * nsy + 2 * NAs_z * nsz) * NBr)
        + NAsNBr * (rho0 * inflow_s - tau_I) );
    }
  }
}

void PLocAssem_VMS_NS_GenAlpha_Interface::Assem_Tangent_itf_MF_rotated(
  const int &qua, const double &rotated_qw, const double &dt,
  const FEAElement * const &rotated_elementv, const FEAElement * const &fixed_elementv,
  const double * const &rotated_local_sol, const double * const &fixed_local_sol,
  const double * const &rotated_local_mvelo)
{
  Zero_Tangent_rs();
  double ps {0.0}, us {0.0}, vs {0.0}, ws {0.0};
  double pr {0.0}, ur {0.0}, vr {0.0}, wr {0.0};
  double mur {0.0}, mvr {0.0}, mwr{0.0};

  std::vector<double> Ns(nLocBas, 0.0), dNs_dx(nLocBas, 0.0), dNs_dy(nLocBas, 0.0), dNs_dz(nLocBas, 0.0);
  std::vector<double> Nr(nLocBas, 0.0), dNr_dx(nLocBas, 0.0), dNr_dy(nLocBas, 0.0), dNr_dz(nLocBas, 0.0);

  fixed_elementv -> get_R_gradR( 0, &Ns[0], &dNs_dx[0], &dNs_dy[0], &dNs_dz[0] );
  rotated_elementv -> get_R_gradR( qua, &Nr[0], &dNr_dx[0], &dNr_dy[0], &dNr_dz[0] );

  double rotated_J {0.0};
  const Vector_3 normal_r = rotated_elementv -> get_2d_normal_out(qua, rotated_J);
  const Vector_3 normal_s = -1 * normal_r;

  // Calculate h_b and tau_I
  const auto s_dxi_dx = fixed_elementv -> get_invJacobian(0);
  const double h_s = get_h_b(s_dxi_dx, normal_s);

  const auto r_dxi_dx = rotated_elementv -> get_invJacobian(qua);
  const double h_r = get_h_b(r_dxi_dx, normal_r);

  const double tau_I = 0.5 * vis_mu * C_bI * (1.0 / h_s + 1.0 / h_r);

  for(int ii{0}; ii<nLocBas; ++ii)
  {
    const int ii4{4 * ii};
    const int ii3{3 * ii};

    ps += fixed_local_sol[ii4 + 0] * Ns[ii];
    us += fixed_local_sol[ii4 + 1] * Ns[ii];
    vs += fixed_local_sol[ii4 + 2] * Ns[ii];
    ws += fixed_local_sol[ii4 + 3] * Ns[ii];

    pr += rotated_local_sol[ii4 + 0] * Nr[ii];
    ur += rotated_local_sol[ii4 + 1] * Nr[ii];
    vr += rotated_local_sol[ii4 + 2] * Nr[ii];
    wr += rotated_local_sol[ii4 + 3] * Nr[ii];

    mur += rotated_local_mvelo[ii3 + 0] * Nr[ii];
    mvr += rotated_local_mvelo[ii3 + 1] * Nr[ii];
    mwr += rotated_local_mvelo[ii3 + 2] * Nr[ii];
  }

  // Mesh velocity in the quadrature point
  const Vector_3 velo_mesh = Vector_3(mur, mvr, mwr);

  const Vector_3 velo_jump(us - ur, vs - vr, ws - wr);
  const double nsx {normal_s.x()}, nsy {normal_s.y()}, nsz {normal_s.z()};
  const double nrx {normal_r.x()}, nry {normal_r.y()}, nrz {normal_r.z()};

  double inflow_r = ((ur - velo_mesh.x()) * nrx + (vr - velo_mesh.y()) * nry + (wr - velo_mesh.z()) * nrz ?
                     (ur - velo_mesh.x()) * nrx + (vr - velo_mesh.y()) * nry + (wr - velo_mesh.z()) * nrz : 0.0);

  const double gwts = rotated_J * rotated_qw;
  const double dd_dv = alpha_f * gamma * dt;
  const double common_coef = gwts * dd_dv;

  for(int A{0}; A<nLocBas; ++A)
  {
    const double NAr {Nr[A]}, NAr_x {dNr_dx[A]}, NAr_y {dNr_dy[A]}, NAr_z {dNr_dz[A]};

    const int A4 = 4 * A;

    for(int B{0}; B<nLocBas; ++B)
    {
      const double NBs {Ns[B]}, NBs_x {dNs_dx[B]}, NBs_y {dNs_dy[B]}, NBs_z {dNs_dz[B]};

      const int B4 = 4 * B, nLocBas4 = 4 * nLocBas;

      const double NArNBs = NAr * NBs;

      // Kr1s1
      // Tangent_rs[nLocBas4 * A4 + B4] += 0;

      // Kr1s2
      Tangent_rs[nLocBas4 * A4 + B4 + 1] += 0.5 * common_coef * NArNBs * nrx;

      // Kr1s3
      Tangent_rs[nLocBas4 * A4 + B4 + 2] += 0.5 * common_coef * NArNBs * nry;

      // Kr1s4
      Tangent_rs[nLocBas4 * A4 + B4 + 3] += 0.5 * common_coef * NArNBs * nrz;

      // Kr2s1
      Tangent_rs[nLocBas4 * (A4 + 1) + B4] += -0.5 * common_coef * NArNBs * nsx;

      // Kr2s2
      Tangent_rs[nLocBas4 * (A4 + 1) + B4 + 1] += common_coef * (
        0.5 * vis_mu * (NAr * (2 * NBs_x * nsx + NBs_y * nsy + NBs_z * nsz)
                     + (2 * NAr_x * nrx + NAr_y * nry + NAr_z * nrz) * NBs)
        + NArNBs * (rho0 * inflow_r - tau_I) );

      // Kr2s3
      Tangent_rs[nLocBas4 * (A4 + 1) + B4 + 2] += 0.5 * common_coef * vis_mu * (
        NAr * NBs_x * nsy + NAr_y * nrx * NBs );

      // Kr2s4
      Tangent_rs[nLocBas4 * (A4 + 1) + B4 + 3] += 0.5 * common_coef * vis_mu * (
        NAr * NBs_x * nsz + NAr_z * nrx * NBs );

      // Kr3s1
      Tangent_rs[nLocBas4 * (A4 + 2) + B4] += -0.5 * common_coef * NArNBs * nsy;

      // Kr3s2
      Tangent_rs[nLocBas4 * (A4 + 2) + B4 + 1] += 0.5 * common_coef * vis_mu * (
        NAr * NBs_y * nsx + NAr_x * nry * NBs );

      // Kr3s3
      Tangent_rs[nLocBas4 * (A4 + 2) + B4 + 2] += common_coef * (
        0.5 * vis_mu * (NAr * (NBs_x * nsx + 2 * NBs_y * nsy + NBs_z * nsz)
                     + (NAr_x * nrx + 2 * NAr_y * nry + NAr_z * nrz) * NBs)
        + NArNBs * (rho0 * inflow_r - tau_I) );
      
      // Kr3s4
      Tangent_rs[nLocBas4 * (A4 + 2) + B4 + 3] += 0.5 * common_coef * vis_mu * (
        NAr * NBs_y * nsz + NAr_z * nry * NBs );

      // Kr4s1
      Tangent_rs[nLocBas4 * (A4 + 3) + B4] += -0.5 * common_coef * NArNBs * nsz;

      // Kr4s2
      Tangent_rs[nLocBas4 * (A4 + 3) + B4 + 1] += 0.5 * common_coef * vis_mu * (
        NAr * NBs_z * nsx + NAr_x * nrz * NBs );

      // Kr4s3
      Tangent_rs[nLocBas4 * (A4 + 3) + B4 + 2] += 0.5 * common_coef * vis_mu * (
        NAr * NBs_z * nsy + NAr_y * nrz * NBs );

      // Kr4s4
      Tangent_rs[nLocBas4 * (A4 + 3) + B4 + 3] += common_coef * (
        0.5 * vis_mu * (NAr * (NBs_x * nsx + NBs_y * nsy + 2 * NBs_z * nsz)
                     + (NAr_x * nrx + NAr_y * nry + 2 * NAr_z * nrz) * NBs)
        + NArNBs * (rho0 * inflow_r - tau_I) );
    }
  }
}
