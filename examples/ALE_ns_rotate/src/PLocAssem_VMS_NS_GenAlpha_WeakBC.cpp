#include "PLocAssem_VMS_NS_GenAlpha_WeakBC.hpp"

PLocAssem_VMS_NS_GenAlpha_WeakBC::PLocAssem_VMS_NS_GenAlpha_WeakBC(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas,
        const double &in_rho, const double &in_vis_mu,
        const double &in_beta, const int &elemtype,
        const double &angular,
        const Vector_3 &point_xyz, const Vector_3 &angular_direc,
        const double &in_ct, const double &in_ctauc,
        const double &in_C_bI )
: PLocAssem_VMS_NS_GenAlpha(tm_gAlpha, in_nlocbas, in_nqp, in_snlocbas,
  in_rho, in_vis_mu, in_beta, elemtype, angular, point_xyz, angular_direc, in_ct, in_ctauc), C_bI(in_C_bI)
{ }

PLocAssem_VMS_NS_GenAlpha_WeakBC::~PLocAssem_VMS_NS_GenAlpha_WeakBC()
{ }

void PLocAssem_VMS_NS_GenAlpha_WeakBC::print_info() const
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
  SYS_T::commPrint("  Penalty para C_bI = %e \n", C_bI);
  SYS_T::commPrint("  Note: \n");
  SYS_T::commPrint("  1. Consistent tangent matrix used. \n");
  SYS_T::commPrint("  2. Nonlinear quadratic term is in advective form. \n");
  SYS_T::commPrint("  3. Pressure is evaluated at n+alpha_f rather than n+1. \n");
  SYS_T::commPrint("----------------------------------------------------------- \n");
}

void PLocAssem_VMS_NS_GenAlpha_WeakBC::Assem_Residual_Weak(
    const double &time, const double &dt,
    const double * const &sol,
    FEAElement * const &elementvs,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quads,
    const int &face_id)
{
  // Build the basis function of volume element
  elementvs->buildBasis( face_id, quads, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double curr {time + alpha_f * dt};

  const int face_nqp {quads -> get_num_quadPts()};

  Zero_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua {0}; qua < face_nqp; ++qua)
  {
    elementvs->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );

    double surface_area {0.0};

    // Calculate surface Jacobian and normal_qua
    const Vector_3 n_out = elementvs->get_2d_normal_out(qua, surface_area);

    Vector_3 coor(0.0, 0.0, 0.0);

    // Calculate u_qua and grad_u_qua
    double p {0.0}, u {0.0}, v {0.0}, w {0.0};
    double u_x {0.0}, v_x {0.0}, w_x {0.0};
    double u_y {0.0}, v_y {0.0}, w_y {0.0};
    double u_z {0.0}, v_z {0.0}, w_z {0.0};
    for(int ii {0}; ii < nLocBas; ++ii)
    {
      const int ii4 {ii * 4};
      p += sol[ii4 + 0] * R[ii];
      u += sol[ii4 + 1] * R[ii];
      v += sol[ii4 + 2] * R[ii];
      w += sol[ii4 + 3] * R[ii];

      u_x += sol[ii4 + 1] * dR_dx[ii];
      v_x += sol[ii4 + 2] * dR_dx[ii];
      w_x += sol[ii4 + 3] * dR_dx[ii];

      u_y += sol[ii4 + 1] * dR_dy[ii];
      v_y += sol[ii4 + 2] * dR_dy[ii];
      w_y += sol[ii4 + 3] * dR_dy[ii];

      u_z += sol[ii4 + 1] * dR_dz[ii];
      v_z += sol[ii4 + 2] * dR_dz[ii];
      w_z += sol[ii4 + 3] * dR_dz[ii];

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    const Vector_3 u_vec (u, v, w);
    const double u_dot_n = u_vec.dot_product(n_out);
    const double inflow_factor = ( u_dot_n < 0.0 ? u_dot_n : 0.0 );

    // Calculate the g_qua
    const Vector_3 u_minus_g = u_vec - get_g_weak(coor, curr);

    // Calculate h_b and tau_B
    const auto dxi_dx = elementvs->get_invJacobian(qua);
    const double h_b = get_h_b(dxi_dx, n_out);

    const Vector_3 u_tan = u_vec - u_dot_n * n_out;
    const double tau_B = get_tau_B(u_tan, h_b / C_bI, vis_mu / rho0);

    // Assembly
    const double gwts = surface_area * quads->get_qw(qua);

    const double penalty_wall = C_bI * vis_mu / h_b - rho0 * tau_B;

    for( int A {0}; A < nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      const int A4 = 4 * A;

      Residual[A4] += gwts * (-1.0) * NA * u_minus_g.dot_product(n_out);

      Residual[A4 + 1] += gwts * ( NA * (rho0 * u_dot_n * u + p * n_out.x()
                                  - 2 * vis_mu * u_x * n_out.x() 
                                  - vis_mu * (u_y + v_x) * n_out.y() 
                                  - vis_mu * (u_z + w_x) * n_out.z())
          - (2 * NA_x * n_out.x() + NA_y * n_out.y() + NA_z * n_out.z()) * vis_mu * u_minus_g.x()
          - NA_y * n_out.x() * vis_mu * u_minus_g.y()
          - NA_z * n_out.x() * vis_mu * u_minus_g.z()
          + NA * rho0 * (tau_B - inflow_factor) * u_minus_g.x()
          + NA * n_out.x() * penalty_wall * u_minus_g.dot_product(n_out) );

      Residual[A4 + 2] += gwts * ( NA * (rho0 * u_dot_n * v + p * n_out.y()
                                  - vis_mu * (v_x + u_y) * n_out.x()
                                  - 2 * vis_mu * v_y * n_out.y()
                                  - vis_mu * (v_z + w_y) * n_out.z())
          - NA_x * n_out.y() * vis_mu * u_minus_g.x()
          - (NA_x * n_out.x() + 2 * NA_y * n_out.y() + NA_z * n_out.z()) * vis_mu * u_minus_g.y()
          - NA_z * n_out.y() * vis_mu * u_minus_g.z()
          + NA * rho0 * (tau_B - inflow_factor) * u_minus_g.y()
          + NA * n_out.y() * penalty_wall * u_minus_g.dot_product(n_out) );

      Residual[A4 + 3] += gwts * ( NA * (rho0 * u_dot_n * w + p * n_out.z()
                                  - vis_mu * (w_x + u_z) * n_out.x()
                                  - vis_mu * (w_y + v_z) * n_out.y()
                                  - 2 * vis_mu * w_z * n_out.z())
          - NA_x * n_out.z() * vis_mu * u_minus_g.x()
          - NA_y * n_out.z() * vis_mu * u_minus_g.y()
          - (NA_x * n_out.x() + NA_y * n_out.y() + 2 * NA_z * n_out.z()) * vis_mu * u_minus_g.z()
          + NA * rho0 * (tau_B - inflow_factor) * u_minus_g.z()
          + NA * n_out.z() * penalty_wall * u_minus_g.dot_product(n_out) );
    } // A-loop
  } // qua-loop
}

void PLocAssem_VMS_NS_GenAlpha_WeakBC::Assem_Tangent_Residual_Weak(
    const double &time, const double &dt,
    const double * const &sol,
    FEAElement * const &elementvs,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quads,
    const int &face_id)
{
  elementvs->buildBasis( face_id, quads, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double curr {time + alpha_f * dt};

  const double dd_dv {alpha_f * gamma * dt};

  const int face_nqp {quads -> get_num_quadPts()};

  Zero_Tangent_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua {0}; qua < face_nqp; ++qua)
  {
    elementvs->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );

    double surface_area {0.0};

    // Calculate surface Jacobian and normal_qua
    const Vector_3 n_out = elementvs->get_2d_normal_out(qua, surface_area);

    Vector_3 coor(0.0, 0.0, 0.0);

    // Calculate u_qua and grad_u_qua
    double p {0.0}, u {0.0}, v {0.0}, w {0.0};
    double u_x {0.0}, v_x {0.0}, w_x {0.0};
    double u_y {0.0}, v_y {0.0}, w_y {0.0};
    double u_z {0.0}, v_z {0.0}, w_z {0.0};
    for(int ii {0}; ii < nLocBas; ++ii)
    {
      const int ii4 {ii * 4};
      p += sol[ii4 + 0] * R[ii];
      u += sol[ii4 + 1] * R[ii];
      v += sol[ii4 + 2] * R[ii];
      w += sol[ii4 + 3] * R[ii];

      u_x += sol[ii4 + 1] * dR_dx[ii];
      v_x += sol[ii4 + 2] * dR_dx[ii];
      w_x += sol[ii4 + 3] * dR_dx[ii];

      u_y += sol[ii4 + 1] * dR_dy[ii];
      v_y += sol[ii4 + 2] * dR_dy[ii];
      w_y += sol[ii4 + 3] * dR_dy[ii];

      u_z += sol[ii4 + 1] * dR_dz[ii];
      v_z += sol[ii4 + 2] * dR_dz[ii];
      w_z += sol[ii4 + 3] * dR_dz[ii];

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    const Vector_3 u_vec (u, v, w);
    const double u_dot_n = u_vec.dot_product(n_out);
    const double inflow_factor = (u_dot_n < 0.0 ? u_dot_n : 0.0);
    const double inflow_flag   = (u_dot_n < 0.0 ? 1.0 : 0.0);
    
    // Calculate the g_qua
    const Vector_3 u_minus_g = u_vec - get_g_weak(coor, curr);

    // Calculate h_b and tau_B
    const auto dxi_dx = elementvs->get_invJacobian(qua);
    const double h_b = get_h_b(dxi_dx, n_out);

    const Vector_3 u_tan = u_vec - u_dot_n * n_out;
    const double tau_B = get_tau_B(u_tan, h_b / C_bI, vis_mu / rho0);

    // Assembly
    const double gwts = surface_area * quads->get_qw(qua);

    const double common_coef = gwts * dd_dv;

    const double penalty_wall = C_bI * vis_mu / h_b - rho0 * tau_B;

    for( int A {0}; A < nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      const int A4 = 4 * A;

      Residual[A4] += gwts * (-1.0) * NA * u_minus_g.dot_product(n_out);

      Residual[A4 + 1] += gwts * ( NA * (rho0 * u_dot_n * u + p * n_out.x()
                                  - 2 * vis_mu * u_x * n_out.x() 
                                  - vis_mu * (u_y + v_x) * n_out.y() 
                                  - vis_mu * (u_z + w_x) * n_out.z())
          - (2 * NA_x * n_out.x() + NA_y * n_out.y() + NA_z * n_out.z()) * vis_mu * u_minus_g.x()
          - NA_y * n_out.x() * vis_mu * u_minus_g.y()
          - NA_z * n_out.x() * vis_mu * u_minus_g.z()
          + NA * rho0 * (tau_B - inflow_factor) * u_minus_g.x()
          + NA * n_out.x() * penalty_wall * u_minus_g.dot_product(n_out) );

      Residual[A4 + 2] += gwts * ( NA * (rho0 * u_dot_n * v + p * n_out.y()
                                  - vis_mu * (v_x + u_y) * n_out.x()
                                  - 2 * vis_mu * v_y * n_out.y()
                                  - vis_mu * (v_z + w_y) * n_out.z())
          - NA_x * n_out.y() * vis_mu * u_minus_g.x()
          - (NA_x * n_out.x() + 2 * NA_y * n_out.y() + NA_z * n_out.z()) * vis_mu * u_minus_g.y()
          - NA_z * n_out.y() * vis_mu * u_minus_g.z()
          + NA * rho0 * (tau_B - inflow_factor) * u_minus_g.y()
          + NA * n_out.y() * penalty_wall * u_minus_g.dot_product(n_out) );

      Residual[A4 + 3] += gwts * ( NA * (rho0 * u_dot_n * w + p * n_out.z()
                                  - vis_mu * (w_x + u_z) * n_out.x()
                                  - vis_mu * (w_y + v_z) * n_out.y()
                                  - 2 * vis_mu * w_z * n_out.z())
          - NA_x * n_out.z() * vis_mu * u_minus_g.x()
          - NA_y * n_out.z() * vis_mu * u_minus_g.y()
          - (NA_x * n_out.x() + NA_y * n_out.y() + 2 * NA_z * n_out.z()) * vis_mu * u_minus_g.z()
          + NA * rho0 * (tau_B - inflow_factor) * u_minus_g.z()
          + NA * n_out.z() * penalty_wall * u_minus_g.dot_product(n_out) );

      for (int B {0}; B < nLocBas; ++B )
      {
        const double NB = R[B], NB_x = dR_dx[B], NB_y = dR_dy[B], NB_z = dR_dz[B];

        const double NANB = NA * NB;

        const int B4 = 4 * B, nLocBas4 = 4 * nLocBas;

        // K11
        // Tangent[nLocBas4 * A4 + B4] += 0;

        // K12
        Tangent[nLocBas4 * A4 + B4 + 1] += common_coef * (-1.0) * NANB * n_out.x();

        // K13
        Tangent[nLocBas4 * A4 + B4 + 2] += common_coef * (-1.0) * NANB * n_out.y();

        // K14
        Tangent[nLocBas4 * A4 + B4 + 3] += common_coef * (-1.0) * NANB * n_out.z();

        // K21
        Tangent[nLocBas4 * (A4 + 1) + B4] += common_coef * NANB * n_out.x();

        // K22
        Tangent[nLocBas4 * (A4 + 1) + B4 + 1] += common_coef * ( NANB * rho0 * (n_out.x() * u + u_dot_n)
            - vis_mu * NA * (2 * n_out.x() * NB_x + n_out.y() * NB_y + n_out.z() * NB_z)
            - vis_mu * NB * (2 * n_out.x() * NA_x + n_out.y() * NA_y + n_out.z() * NA_z)
            - NANB * rho0 * (inflow_flag * u_minus_g.x() * n_out.x() + inflow_factor)
            + NANB * (rho0 * tau_B + n_out.x() * n_out.x() * penalty_wall) );

        // K23
        Tangent[nLocBas4 * (A4 + 1) + B4 + 2] += common_coef * ( NANB * rho0 * n_out.y() * u
            - vis_mu * NA * n_out.y() * NB_x
            - vis_mu * NB * n_out.x() * NA_y
            - NANB * rho0 * inflow_flag * u_minus_g.x() * n_out.y()
            + NANB * n_out.x() * n_out.y() * penalty_wall );

        // K24
        Tangent[nLocBas4 * (A4 + 1) + B4 + 3] += common_coef * ( NANB * rho0 * n_out.z() * u
            - vis_mu * NA * n_out.z() * NB_x
            - vis_mu * NB * n_out.x() * NA_z
            - NANB * rho0 * inflow_flag * u_minus_g.x() * n_out.z()
            + NANB * n_out.x() * n_out.z() * penalty_wall );

        // K31
        Tangent[nLocBas4 * (A4 + 2) + B4] += common_coef * NANB * n_out.y();

        // K32
        Tangent[nLocBas4 * (A4 + 2) + B4 + 1] += common_coef * ( NANB * rho0 * n_out.x() * v
            - vis_mu * NA * n_out.x() * NB_y
            - vis_mu * NB * n_out.y() * NA_x
            - NANB * rho0 * inflow_flag * u_minus_g.y() * n_out.x()
            + NANB * n_out.x() * n_out.y() * penalty_wall );

        // K33
        Tangent[nLocBas4 * (A4 + 2) + B4 + 2] += common_coef * ( NANB * rho0 * (n_out.y() * v + u_dot_n)
            - vis_mu * NA * (n_out.x() * NB_x + 2 * n_out.y() * NB_y + n_out.z() * NB_z)
            - vis_mu * NB * (n_out.x() * NA_x + 2 * n_out.y() * NA_y + n_out.z() * NA_z)
            - NANB * rho0 * (inflow_flag * u_minus_g.y() * n_out.y() + inflow_factor)
            + NANB * (rho0 * tau_B + n_out.y() * n_out.y() * penalty_wall) );
        
        // K34
        Tangent[nLocBas4 * (A4 + 2) + B4 + 3] += common_coef * ( NANB * rho0 * n_out.z() * v
            - vis_mu * NA * n_out.z() * NB_y
            - vis_mu * NB * n_out.y() * NA_z
            - NANB * rho0 * inflow_flag * u_minus_g.y() * n_out.z()
            + NANB * n_out.y() * n_out.z() * penalty_wall );

        // K41
        Tangent[nLocBas4 * (A4 + 3) + B4] += common_coef * NANB * n_out.z();

        // K42
        Tangent[nLocBas4 * (A4 + 3) + B4 + 1] += common_coef * ( NANB * rho0 * n_out.x() * w
            - vis_mu * NA * n_out.x() * NB_z
            - vis_mu * NB * n_out.z() * NA_x
            - NANB * rho0 * inflow_flag * u_minus_g.z() * n_out.x()
            + NANB * n_out.x() * n_out.z() * penalty_wall );

        // K43
        Tangent[nLocBas4 * (A4 + 3) + B4 + 2] += common_coef * ( NANB * rho0 * n_out.y() * w
            - vis_mu * NA * n_out.y() * NB_z
            - vis_mu * NB * n_out.z() * NA_y
            - NANB * rho0 * inflow_flag * u_minus_g.z() * n_out.y()
            + NANB * n_out.y() * n_out.z() * penalty_wall );

        // K44
        Tangent[nLocBas4 * (A4 + 3) + B4 + 3] += common_coef * ( NANB * rho0 * (n_out.z() * w + u_dot_n)
            - vis_mu * NA * (n_out.x() * NB_x + n_out.y() * NB_y + 2 * n_out.z() * NB_z)
            - vis_mu * NB * (n_out.x() * NA_x + n_out.y() * NA_y + 2 * n_out.z() * NA_z)
            - NANB * rho0 * (inflow_flag * u_minus_g.z() * n_out.z() + inflow_factor)
            + NANB * (rho0 * tau_B + n_out.z() * n_out.z() * penalty_wall) );
      } // B-loop
    } // A-loop
  } // qua-loop
}

void PLocAssem_VMS_NS_GenAlpha_WeakBC::Assem_Residual_Weak_Rotated(
    const double &time, const double &dt,
    const double * const &sol,
    FEAElement * const &elementvs,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quads,
    const int &face_id)
{
  const double curr {time + alpha_f * dt};

  std::vector<double> curPt_x(nLocBas, 0.0), curPt_y(nLocBas, 0.0), curPt_z(nLocBas, 0.0);

  //Update coordinates
  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, curr, &curPt_x[0], &curPt_y[0], &curPt_z[0], 0);

  // Build the basis function of volume element
  elementvs->buildBasis( face_id, quads, &curPt_x[0], &curPt_y[0], &curPt_z[0] );

  const int face_nqp {quads -> get_num_quadPts()};

  Zero_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua {0}; qua < face_nqp; ++qua)
  {
    elementvs->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );

    double surface_area {0.0};

    // Calculate surface Jacobian and normal_qua
    const Vector_3 n_out = elementvs->get_2d_normal_out(qua, surface_area);

    Vector_3 coor(0.0, 0.0, 0.0);

    // Calculate u_qua and grad_u_qua
    double p {0.0}, u {0.0}, v {0.0}, w {0.0};
    double u_x {0.0}, v_x {0.0}, w_x {0.0};
    double u_y {0.0}, v_y {0.0}, w_y {0.0};
    double u_z {0.0}, v_z {0.0}, w_z {0.0};
    for(int ii {0}; ii < nLocBas; ++ii)
    {
      const int ii4 {ii * 4};
      p += sol[ii4 + 0] * R[ii];
      u += sol[ii4 + 1] * R[ii];
      v += sol[ii4 + 2] * R[ii];
      w += sol[ii4 + 3] * R[ii];

      u_x += sol[ii4 + 1] * dR_dx[ii];
      v_x += sol[ii4 + 2] * dR_dx[ii];
      w_x += sol[ii4 + 3] * dR_dx[ii];

      u_y += sol[ii4 + 1] * dR_dy[ii];
      v_y += sol[ii4 + 2] * dR_dy[ii];
      w_y += sol[ii4 + 3] * dR_dy[ii];

      u_z += sol[ii4 + 1] * dR_dz[ii];
      v_z += sol[ii4 + 2] * dR_dz[ii];
      w_z += sol[ii4 + 3] * dR_dz[ii];

      coor.x() += curPt_x[ii] * R[ii];
      coor.y() += curPt_y[ii] * R[ii];
      coor.z() += curPt_z[ii] * R[ii];
    }

    const Vector_3 u_vec (u, v, w);

    // Mesh velocity in the quadrature point
    const Vector_3 radius_qua = get_radius(coor);
    const Vector_3 velo_mesh = Vec3::cross_product(angular_velo, radius_qua);

    // v - hat(v)
    const Vector_3 c_vec = u_vec - velo_mesh;

    // Consider the mesh velocity
    const double c_dot_n = c_vec.dot_product(n_out);
    const double inflow_factor = ( c_dot_n < 0.0 ? c_dot_n : 0.0 );

    // Calculate the g_qua
    const Vector_3 u_minus_g = u_vec - get_g_weak(coor, curr);

    // Calculate h_b and tau_B
    const auto dxi_dx = elementvs->get_invJacobian(qua);
    const double h_b = get_h_b(dxi_dx, n_out);

    const Vector_3 c_tan = c_vec - c_dot_n * n_out;
    const double tau_B = get_tau_B(c_tan, h_b / C_bI, vis_mu / rho0);

    // Assembly
    const double gwts = surface_area * quads->get_qw(qua);

    const double penalty_wall = C_bI * vis_mu / h_b - rho0 * tau_B;

    for( int A {0}; A < nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      const int A4 = 4 * A;

      Residual[A4] += gwts * (-1.0) * NA * u_minus_g.dot_product(n_out);

      Residual[A4 + 1] += gwts * ( NA * (rho0 * c_dot_n * u + p * n_out.x()
                                  - 2 * vis_mu * u_x * n_out.x() 
                                  - vis_mu * (u_y + v_x) * n_out.y() 
                                  - vis_mu * (u_z + w_x) * n_out.z())
          - (2 * NA_x * n_out.x() + NA_y * n_out.y() + NA_z * n_out.z()) * vis_mu * u_minus_g.x()
          - NA_y * n_out.x() * vis_mu * u_minus_g.y()
          - NA_z * n_out.x() * vis_mu * u_minus_g.z()
          + NA * rho0 * (tau_B - inflow_factor) * u_minus_g.x()
          + NA * n_out.x() * penalty_wall * u_minus_g.dot_product(n_out) );

      Residual[A4 + 2] += gwts * ( NA * (rho0 * c_dot_n * v + p * n_out.y()
                                  - vis_mu * (v_x + u_y) * n_out.x()
                                  - 2 * vis_mu * v_y * n_out.y()
                                  - vis_mu * (v_z + w_y) * n_out.z())
          - NA_x * n_out.y() * vis_mu * u_minus_g.x()
          - (NA_x * n_out.x() + 2 * NA_y * n_out.y() + NA_z * n_out.z()) * vis_mu * u_minus_g.y()
          - NA_z * n_out.y() * vis_mu * u_minus_g.z()
          + NA * rho0 * (tau_B - inflow_factor) * u_minus_g.y()
          + NA * n_out.y() * penalty_wall * u_minus_g.dot_product(n_out) );

      Residual[A4 + 3] += gwts * ( NA * (rho0 * c_dot_n * w + p * n_out.z()
                                  - vis_mu * (w_x + u_z) * n_out.x()
                                  - vis_mu * (w_y + v_z) * n_out.y()
                                  - 2 * vis_mu * w_z * n_out.z())
          - NA_x * n_out.z() * vis_mu * u_minus_g.x()
          - NA_y * n_out.z() * vis_mu * u_minus_g.y()
          - (NA_x * n_out.x() + NA_y * n_out.y() + 2 * NA_z * n_out.z()) * vis_mu * u_minus_g.z()
          + NA * rho0 * (tau_B - inflow_factor) * u_minus_g.z()
          + NA * n_out.z() * penalty_wall * u_minus_g.dot_product(n_out) );
    } // A-loop
  } // qua-loop
}

void PLocAssem_VMS_NS_GenAlpha_WeakBC::Assem_Tangent_Residual_Weak_Rotated(
    const double &time, const double &dt,
    const double * const &sol,
    FEAElement * const &elementvs,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quads,
    const int &face_id)
{
  const double curr {time + alpha_f * dt};

  std::vector<double> curPt_x(nLocBas, 0.0), curPt_y(nLocBas, 0.0), curPt_z(nLocBas, 0.0);

  //Update coordinates
  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, curr, &curPt_x[0], &curPt_y[0], &curPt_z[0], 0);

  // Build the basis function of volume element
  elementvs->buildBasis( face_id, quads, &curPt_x[0], &curPt_y[0], &curPt_z[0] );

  const double dd_dv {alpha_f * gamma * dt};

  const int face_nqp {quads -> get_num_quadPts()};

  Zero_Tangent_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua {0}; qua < face_nqp; ++qua)
  {
    elementvs->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );

    double surface_area {0.0};

    // Calculate surface Jacobian and normal_qua
    const Vector_3 n_out = elementvs->get_2d_normal_out(qua, surface_area);

    Vector_3 coor(0.0, 0.0, 0.0);

    // Calculate u_qua and grad_u_qua
    double p {0.0}, u {0.0}, v {0.0}, w {0.0};
    double u_x {0.0}, v_x {0.0}, w_x {0.0};
    double u_y {0.0}, v_y {0.0}, w_y {0.0};
    double u_z {0.0}, v_z {0.0}, w_z {0.0};
    for(int ii {0}; ii < nLocBas; ++ii)
    {
      const int ii4 {ii * 4};
      p += sol[ii4 + 0] * R[ii];
      u += sol[ii4 + 1] * R[ii];
      v += sol[ii4 + 2] * R[ii];
      w += sol[ii4 + 3] * R[ii];

      u_x += sol[ii4 + 1] * dR_dx[ii];
      v_x += sol[ii4 + 2] * dR_dx[ii];
      w_x += sol[ii4 + 3] * dR_dx[ii];

      u_y += sol[ii4 + 1] * dR_dy[ii];
      v_y += sol[ii4 + 2] * dR_dy[ii];
      w_y += sol[ii4 + 3] * dR_dy[ii];

      u_z += sol[ii4 + 1] * dR_dz[ii];
      v_z += sol[ii4 + 2] * dR_dz[ii];
      w_z += sol[ii4 + 3] * dR_dz[ii];

      coor.x() += curPt_x[ii] * R[ii];
      coor.y() += curPt_y[ii] * R[ii];
      coor.z() += curPt_z[ii] * R[ii];
    }

    const Vector_3 u_vec (u, v, w);

    // Mesh velocity in the quadrature point
    const Vector_3 radius_qua = get_radius(coor);
    const Vector_3 velo_mesh = Vec3::cross_product(angular_velo, radius_qua);

    // v - hat(v)
    const Vector_3 c_vec = u_vec - velo_mesh;

    // Consider the mesh velocity
    const double c_dot_n = c_vec.dot_product(n_out);
    const double inflow_factor = ( c_dot_n < 0.0 ? c_dot_n : 0.0 );
    const double inflow_flag   = (c_dot_n < 0.0 ? 1.0 : 0.0);

    // Calculate the g_qua
    const Vector_3 u_minus_g = u_vec - get_g_weak(coor, curr);

    // Calculate h_b and tau_B
    const auto dxi_dx = elementvs->get_invJacobian(qua);
    const double h_b = get_h_b(dxi_dx, n_out);

    const Vector_3 c_tan = c_vec - c_dot_n * n_out;
    const double tau_B = get_tau_B(c_tan, h_b / C_bI, vis_mu / rho0);

    // Assembly
    const double gwts = surface_area * quads->get_qw(qua);

    const double common_coef = gwts * dd_dv;

    const double penalty_wall = C_bI * vis_mu / h_b - rho0 * tau_B;

    for( int A {0}; A < nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      const int A4 = 4 * A;

      Residual[A4] += gwts * (-1.0) * NA * u_minus_g.dot_product(n_out);

      Residual[A4 + 1] += gwts * ( NA * (rho0 * c_dot_n * u + p * n_out.x()
                                  - 2 * vis_mu * u_x * n_out.x() 
                                  - vis_mu * (u_y + v_x) * n_out.y() 
                                  - vis_mu * (u_z + w_x) * n_out.z())
          - (2 * NA_x * n_out.x() + NA_y * n_out.y() + NA_z * n_out.z()) * vis_mu * u_minus_g.x()
          - NA_y * n_out.x() * vis_mu * u_minus_g.y()
          - NA_z * n_out.x() * vis_mu * u_minus_g.z()
          + NA * rho0 * (tau_B - inflow_factor) * u_minus_g.x()
          + NA * n_out.x() * penalty_wall * u_minus_g.dot_product(n_out) );

      Residual[A4 + 2] += gwts * ( NA * (rho0 * c_dot_n * v + p * n_out.y()
                                  - vis_mu * (v_x + u_y) * n_out.x()
                                  - 2 * vis_mu * v_y * n_out.y()
                                  - vis_mu * (v_z + w_y) * n_out.z())
          - NA_x * n_out.y() * vis_mu * u_minus_g.x()
          - (NA_x * n_out.x() + 2 * NA_y * n_out.y() + NA_z * n_out.z()) * vis_mu * u_minus_g.y()
          - NA_z * n_out.y() * vis_mu * u_minus_g.z()
          + NA * rho0 * (tau_B - inflow_factor) * u_minus_g.y()
          + NA * n_out.y() * penalty_wall * u_minus_g.dot_product(n_out) );

      Residual[A4 + 3] += gwts * ( NA * (rho0 * c_dot_n * w + p * n_out.z()
                                  - vis_mu * (w_x + u_z) * n_out.x()
                                  - vis_mu * (w_y + v_z) * n_out.y()
                                  - 2 * vis_mu * w_z * n_out.z())
          - NA_x * n_out.z() * vis_mu * u_minus_g.x()
          - NA_y * n_out.z() * vis_mu * u_minus_g.y()
          - (NA_x * n_out.x() + NA_y * n_out.y() + 2 * NA_z * n_out.z()) * vis_mu * u_minus_g.z()
          + NA * rho0 * (tau_B - inflow_factor) * u_minus_g.z()
          + NA * n_out.z() * penalty_wall * u_minus_g.dot_product(n_out) );

      for (int B {0}; B < nLocBas; ++B )
      {
        const double NB = R[B], NB_x = dR_dx[B], NB_y = dR_dy[B], NB_z = dR_dz[B];

        const double NANB = NA * NB;

        const int B4 = 4 * B, nLocBas4 = 4 * nLocBas;

        // K11
        // Tangent[nLocBas4 * A4 + B4] += 0;

        // K12
        Tangent[nLocBas4 * A4 + B4 + 1] += common_coef * (-1.0) * NANB * n_out.x();

        // K13
        Tangent[nLocBas4 * A4 + B4 + 2] += common_coef * (-1.0) * NANB * n_out.y();

        // K14
        Tangent[nLocBas4 * A4 + B4 + 3] += common_coef * (-1.0) * NANB * n_out.z();

        // K21
        Tangent[nLocBas4 * (A4 + 1) + B4] += common_coef * NANB * n_out.x();

        // K22
        Tangent[nLocBas4 * (A4 + 1) + B4 + 1] += common_coef * ( NANB * rho0 * (n_out.x() * u + c_dot_n)
            - vis_mu * NA * (2 * n_out.x() * NB_x + n_out.y() * NB_y + n_out.z() * NB_z)
            - vis_mu * NB * (2 * n_out.x() * NA_x + n_out.y() * NA_y + n_out.z() * NA_z)
            - NANB * rho0 * (inflow_flag * u_minus_g.x() * n_out.x() + inflow_factor)
            + NANB * (rho0 * tau_B + n_out.x() * n_out.x() * penalty_wall) );

        // K23
        Tangent[nLocBas4 * (A4 + 1) + B4 + 2] += common_coef * ( NANB * rho0 * n_out.y() * u
            - vis_mu * NA * n_out.y() * NB_x
            - vis_mu * NB * n_out.x() * NA_y
            - NANB * rho0 * inflow_flag * u_minus_g.x() * n_out.y()
            + NANB * n_out.x() * n_out.y() * penalty_wall );

        // K24
        Tangent[nLocBas4 * (A4 + 1) + B4 + 3] += common_coef * ( NANB * rho0 * n_out.z() * u
            - vis_mu * NA * n_out.z() * NB_x
            - vis_mu * NB * n_out.x() * NA_z
            - NANB * rho0 * inflow_flag * u_minus_g.x() * n_out.z()
            + NANB * n_out.x() * n_out.z() * penalty_wall );

        // K31
        Tangent[nLocBas4 * (A4 + 2) + B4] += common_coef * NANB * n_out.y();

        // K32
        Tangent[nLocBas4 * (A4 + 2) + B4 + 1] += common_coef * ( NANB * rho0 * n_out.x() * v
            - vis_mu * NA * n_out.x() * NB_y
            - vis_mu * NB * n_out.y() * NA_x
            - NANB * rho0 * inflow_flag * u_minus_g.y() * n_out.x()
            + NANB * n_out.x() * n_out.y() * penalty_wall );

        // K33
        Tangent[nLocBas4 * (A4 + 2) + B4 + 2] += common_coef * ( NANB * rho0 * (n_out.y() * v + c_dot_n)
            - vis_mu * NA * (n_out.x() * NB_x + 2 * n_out.y() * NB_y + n_out.z() * NB_z)
            - vis_mu * NB * (n_out.x() * NA_x + 2 * n_out.y() * NA_y + n_out.z() * NA_z)
            - NANB * rho0 * (inflow_flag * u_minus_g.y() * n_out.y() + inflow_factor)
            + NANB * (rho0 * tau_B + n_out.y() * n_out.y() * penalty_wall) );
        
        // K34
        Tangent[nLocBas4 * (A4 + 2) + B4 + 3] += common_coef * ( NANB * rho0 * n_out.z() * v
            - vis_mu * NA * n_out.z() * NB_y
            - vis_mu * NB * n_out.y() * NA_z
            - NANB * rho0 * inflow_flag * u_minus_g.y() * n_out.z()
            + NANB * n_out.y() * n_out.z() * penalty_wall );

        // K41
        Tangent[nLocBas4 * (A4 + 3) + B4] += common_coef * NANB * n_out.z();

        // K42
        Tangent[nLocBas4 * (A4 + 3) + B4 + 1] += common_coef * ( NANB * rho0 * n_out.x() * w
            - vis_mu * NA * n_out.x() * NB_z
            - vis_mu * NB * n_out.z() * NA_x
            - NANB * rho0 * inflow_flag * u_minus_g.z() * n_out.x()
            + NANB * n_out.x() * n_out.z() * penalty_wall );

        // K43
        Tangent[nLocBas4 * (A4 + 3) + B4 + 2] += common_coef * ( NANB * rho0 * n_out.y() * w
            - vis_mu * NA * n_out.y() * NB_z
            - vis_mu * NB * n_out.z() * NA_y
            - NANB * rho0 * inflow_flag * u_minus_g.z() * n_out.y()
            + NANB * n_out.y() * n_out.z() * penalty_wall );

        // K44
        Tangent[nLocBas4 * (A4 + 3) + B4 + 3] += common_coef * ( NANB * rho0 * (n_out.z() * w + c_dot_n)
            - vis_mu * NA * (n_out.x() * NB_x + n_out.y() * NB_y + 2 * n_out.z() * NB_z)
            - vis_mu * NB * (n_out.x() * NA_x + n_out.y() * NA_y + 2 * n_out.z() * NA_z)
            - NANB * rho0 * (inflow_flag * u_minus_g.z() * n_out.z() + inflow_factor)
            + NANB * (rho0 * tau_B + n_out.z() * n_out.z() * penalty_wall) );
      } // B-loop
    } // A-loop
  } // qua-loop
}

// EOF