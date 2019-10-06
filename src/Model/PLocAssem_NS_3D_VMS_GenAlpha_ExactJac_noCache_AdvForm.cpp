#include "PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm.hpp"

PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm::PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm(
    const class TimeMethod_GenAlpha * const &tm_gAlpha,
    const int &in_nlocbas, const int &in_nqp,
    const int &in_sdeg, const int &in_tdeg, const int &in_udeg,
    const int &in_nqpx, const int &in_nqpy, const int &in_nqpz,
    const double &phy_len_x, const double &phy_len_y,
    const double &phy_len_z, const double &max_hx,
    const double &max_hy, const double &max_hz )
: nu(1.0e0), alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), nLocBas(in_nlocbas), dof_per_node(4), 
  vec_size(nLocBas * dof_per_node), sdeg(in_sdeg), tdeg(in_tdeg), udeg(in_udeg),
  nqp(in_nqp), nqpx(in_nqpx), nqpy(in_nqpy), nqpz(in_nqpz), CI(36.0), CT(0.0)
{
  if( nLocBas != (sdeg+1) * (tdeg+1) * (udeg+1) )
  {
    SYS_T::commPrint("Warning : In PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm, nLocBas != (sdeg+1)*(tdeg+1)*(udeg+1). \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  if( nqp != (nqpx * nqpy * nqpz) )
  {
    SYS_T::commPrint("PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm : nqp != nqpx * nqpy * nqpz. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  // Allocate memory layout
  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  for(int ii=0; ii<vec_size; ++ii)
    Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size*vec_size; ++ii)
    Tangent[ii] = 0.0;

  R       = new double [nLocBas];
  dR_dx   = new double [nLocBas];
  dR_dy   = new double [nLocBas];
  dR_dz   = new double [nLocBas];
  d2R_dxx = new double [nLocBas];
  d2R_dyy = new double [nLocBas];
  d2R_dzz = new double [nLocBas];

  dxi_dx  = new double [9]; // 9 components: ds_dx, ds_dy, ..., du_dz

  // 16 := dof_per_node * dof_per_node
  Sub_Tan_m = new double * [16];
  Sub_Tan_f = new double * [16];
  Sub_Tan_p = new double * [16];

  for(int ii=0; ii<16; ++ii)
  {
    Sub_Tan_m[ii] = new double [nLocBas * nLocBas];
    Sub_Tan_f[ii] = new double [nLocBas * nLocBas];
    Sub_Tan_p[ii] = new double [nLocBas * nLocBas];
  }

  // Print the model info on screen
  print_info();
}


PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm::~PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm()
{
  delete [] Tangent;
  delete [] Residual;
  delete [] R; delete [] dR_dx; delete [] dR_dy; delete [] dR_dz;
  delete [] d2R_dxx; delete [] d2R_dyy; delete [] d2R_dzz;
  delete [] dxi_dx;
  for(int ii=0; ii<16; ++ii)
  {
    delete [] Sub_Tan_m[ii];
    delete [] Sub_Tan_f[ii];
    delete [] Sub_Tan_p[ii];
  }
  delete [] Sub_Tan_m; delete [] Sub_Tan_f; delete [] Sub_Tan_p;
}


void PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Three-dimensional Incompressible Navier-Stokes: \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Spatial: Variational Multiscale Methods \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Quadrature info is NOT cached in memory \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Temporal: Generalized-alpha Method \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Viscosity nu = %e \n", nu);
  PetscPrintf(PETSC_COMM_WORLD, "  Stabilization para CI = %e \n", CI);
  PetscPrintf(PETSC_COMM_WORLD, "  Stabilization para CT = %e \n", CT);
  PetscPrintf(PETSC_COMM_WORLD, "  Consistent tangent matrix used. \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Nonlinear quadratic term is in advective form. \n");
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size; ++ii)
    Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size*vec_size; ++ii)
    Tangent[ii] = 0.0;
}


void PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm::Zero_Residual()
{
  for(int ii=0; ii<vec_size; ++ii)
    Residual[ii] = 0.0;
}


void PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm::Assem_Estimate()
{
  for(int ii=0; ii<vec_size*vec_size; ++ii)
    Tangent[ii] = 1.0;
}


void PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm::get_tau(
    double &tau_m_qua, double &tau_c_qua, const double &dt,
    const double * const &dxi_dx, const double &u,
    const double &v, const double &w) const
{
  // dxi_dx has 9 components: 0 ds_dx, 1 ds_dy, 2 ds_dz, 
  //                          3 dt_dx, 4 dt_dy, 5 dt_dz, 
  //                          6 du_dx, 7 du_dy, 8 du_dz.
  const double G11 = dxi_dx[0] * dxi_dx[0] + dxi_dx[3] * dxi_dx[3] + dxi_dx[6] * dxi_dx[6];
  const double G12 = dxi_dx[0] * dxi_dx[1] + dxi_dx[3] * dxi_dx[4] + dxi_dx[6] * dxi_dx[7];
  const double G13 = dxi_dx[0] * dxi_dx[2] + dxi_dx[3] * dxi_dx[5] + dxi_dx[6] * dxi_dx[8];

  //const double G21 = G12;
  const double G22 = dxi_dx[1] * dxi_dx[1] + dxi_dx[4] * dxi_dx[4] + dxi_dx[7] * dxi_dx[7];
  const double G23 = dxi_dx[1] * dxi_dx[2] + dxi_dx[4] * dxi_dx[5] + dxi_dx[7] * dxi_dx[8];

  //const double G31 = G13;
  //const double G32 = G23;
  const double G33 = dxi_dx[2] * dxi_dx[2] + dxi_dx[5] * dxi_dx[5] + dxi_dx[8] * dxi_dx[8]; 

  const double GdG = G11 * G11 + 2.0 * G12 * G12 + 2.0 * G13 * G13
    + G22 * G22 + 2.0 * G23 * G23 + G33 * G33;

  const double uGu = G11 * u * u + 2.0 * G12 * u * v + 2.0 * G13 * u * w
    + G22 * v * v + 2.0 * G23 * v * w + G33 * w * w;

  const double g1 = dxi_dx[0] + dxi_dx[3] + dxi_dx[6];
  const double g2 = dxi_dx[1] + dxi_dx[4] + dxi_dx[7];
  const double g3 = dxi_dx[2] + dxi_dx[5] + dxi_dx[8];

  const double g_dot_g = g1 * g1 + g2 * g2 + g3 * g3;

  const double denom_m = CT / (dt*dt) + uGu + CI * nu * nu * GdG;

  tau_m_qua = 1.0 / sqrt(denom_m);

  const double denom_c = tau_m_qua * g_dot_g;

  tau_c_qua = 1.0 / denom_c;
}


void PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm::Assem_Residual(
    const double &time, const double &dt,
    const double * const &velo,
    const double * const &disp,
    const int &eindex,
    const double &hx, const double &hy, const double &hz,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const BernsteinBasis_Array * const &bu,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &eleCtrlPts_w,
    const double * const &ext_x,
    const double * const &ext_y,
    const double * const &ext_z,
    const class AInt_Weight * const &weight )
{
  // Generate the element basis function value at quad-points
  FEAElement * element = new FEAElement_NURBS_3D_der1_lap_vms( eindex, 
      hx, hy, hz, bs, bt, bu, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z,
      eleCtrlPts_w, ext_x, ext_y, ext_z ); 

  int ii, qua, A; // iterator
  double u, u_t, u_x, u_y, u_z, u_xx, u_yy, u_zz;
  double v, v_t, v_x, v_y, v_z, v_xx, v_yy, v_zz;
  double w, w_t, w_x, w_y, w_z, w_xx, w_yy, w_zz;
  double p, f1, f2, f3;
  //double mf1, mf2, mf3, mf4;
  double p_x, p_y, p_z;
  double gwts, coor_x, coor_y, coor_z;

  // Momentum residual and their derivatives
  double rx, ry, rz;

  // Stabilization parameter
  double tau_m, tau_c;

  // Frequently used quantities in assembly
  const double two_nu = 2.0 * nu;
  double NA, NA_x, NA_y, NA_z;
  double velo_dot_gradR, div_vel;
  double r_dot_gradR; // rx * NA_x + ry * NA_y + rz * NA_z
  double tau_m_2; // tau_m * tau_m

  const double curr = time + alpha_f * dt;

  Zero_Residual();

  for(qua=0; qua<nqp; ++qua)
  {
    u = 0.0; u_t = 0.0; u_x = 0.0; u_y = 0.0; u_z = 0.0;
    u_xx = 0.0; u_yy = 0.0; u_zz = 0.0;
    v = 0.0; v_t = 0.0; v_x = 0.0; v_y = 0.0; v_z = 0.0;
    v_xx = 0.0; v_yy = 0.0; v_zz = 0.0;
    w = 0.0; w_t = 0.0; w_x = 0.0; w_y = 0.0; w_z = 0.0;
    w_xx = 0.0; w_yy = 0.0; w_zz = 0.0;
    p = 0.0; coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    p_x = 0.0; p_y = 0.0; p_z = 0.0;

    element->get_3D_R_gradR_LaplacianR(qua, R, dR_dx, dR_dy, dR_dz,
        d2R_dxx, d2R_dyy, d2R_dzz);
    element->get_invJacobian(qua, dxi_dx); 

    for(ii=0; ii<nLocBas; ++ii)
    {
      u_t += velo[4*ii]   * R[ii];
      v_t += velo[4*ii+1] * R[ii];
      w_t += velo[4*ii+2] * R[ii];

      u += disp[4*ii]   * R[ii];
      v += disp[4*ii+1] * R[ii];
      w += disp[4*ii+2] * R[ii];
      p += disp[4*ii+3] * R[ii];

      u_x += disp[4*ii]   * dR_dx[ii];
      v_x += disp[4*ii+1] * dR_dx[ii];
      w_x += disp[4*ii+2] * dR_dx[ii];
      p_x += disp[4*ii+3] * dR_dx[ii];

      u_y += disp[4*ii]   * dR_dy[ii];
      v_y += disp[4*ii+1] * dR_dy[ii];
      w_y += disp[4*ii+2] * dR_dy[ii];
      p_y += disp[4*ii+3] * dR_dy[ii];

      u_z += disp[4*ii]   * dR_dz[ii];
      v_z += disp[4*ii+1] * dR_dz[ii];
      w_z += disp[4*ii+2] * dR_dz[ii];
      p_z += disp[4*ii+3] * dR_dz[ii];

      u_xx += disp[4*ii]   * d2R_dxx[ii];
      v_xx += disp[4*ii+1] * d2R_dxx[ii];
      w_xx += disp[4*ii+2] * d2R_dxx[ii];

      u_yy += disp[4*ii]   * d2R_dyy[ii];
      v_yy += disp[4*ii+1] * d2R_dyy[ii];
      w_yy += disp[4*ii+2] * d2R_dyy[ii];

      u_zz += disp[4*ii]   * d2R_dzz[ii];
      v_zz += disp[4*ii+1] * d2R_dzz[ii];
      w_zz += disp[4*ii+2] * d2R_dzz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    get_tau(tau_m, tau_c, dt, dxi_dx, u, v, w);

    tau_m_2 = tau_m * tau_m;

    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    get_f(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    //get_mf(coor_x, coor_y, coor_z, curr, mf1, mf2, mf3, mf4);

    rx = u_t + u_x * u + u_y * v + u_z * w + p_x - nu * (u_xx + u_yy + u_zz) - f1; 
    ry = v_t + v_x * u + v_y * v + v_z * w + p_y - nu * (v_xx + v_yy + v_zz) - f2;
    rz = w_t + w_x * u + w_y * v + w_z * w + p_z - nu * (w_xx + w_yy + w_zz) - f3;

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];
      velo_dot_gradR = NA_x * u + NA_y * v + NA_z * w;
      r_dot_gradR = NA_x * rx + NA_y * ry + NA_z * rz;
      div_vel = u_x + v_y + w_z;

      Residual[4*A] += gwts * (NA * u_t 
          + NA * (u*u_x + v*u_y + w * u_z) 
          - NA_x * p
          + NA_x * two_nu * u_x 
          + NA_y * nu * (u_y + v_x) 
          + NA_z * nu * (u_z + w_x)
          + velo_dot_gradR * tau_m * rx 
          + r_dot_gradR * u * tau_m
          + NA_x * tau_c * div_vel 
          - r_dot_gradR * tau_m_2 * rx
          - NA * f1 );

      Residual[4*A+1] += gwts * (NA * v_t 
          + NA * (u * v_x + v * v_y + w * v_z) 
          - NA_y * p
          + NA_x * nu * (u_y + v_x)
          + NA_y * two_nu * v_y
          + NA_z * nu * (v_z + w_y)
          + velo_dot_gradR * tau_m * ry
          + r_dot_gradR * v * tau_m
          + NA_y * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * ry
          - NA * f2 );

      Residual[4*A+2] += gwts * (NA * w_t 
          + NA * (u * w_x + v * w_y + w * w_z) 
          - NA_z * p
          + NA_x * nu * (u_z + w_x)
          + NA_y * nu * (w_y + v_z)
          + NA_z * two_nu * w_z
          + velo_dot_gradR * tau_m * rz
          + r_dot_gradR * w * tau_m 
          + NA_z * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * rz
          - NA * f3 );

      Residual[4*A+3] += gwts * ( NA * div_vel + tau_m * r_dot_gradR );
    }
  }

  delete element;
}


void PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm::Assem_Residual_TopFace(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    const int &eindex,
    const double &hx, const double &hy, const double &hz,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &eleCtrlPts_w,
    const double * const &ext_x,
    const double * const &ext_y,
    const double * const &ext_z
    )
{
  // Build the boundary quadrature info for this element : eindex
  IQuadPts * face_quad_z = new QuadPts_bc1();
  IQuadPts * face_quad_y = new QuadPts_Gauss(nqpy);
  IQuadPts * face_quad_x = new QuadPts_Gauss(nqpx);

  AInt_Weight * face_Int_w = new AInt_Weight(face_quad_x, face_quad_y, face_quad_z);

  BernsteinBasis_Array * ba_s = new BernsteinBasis_Array(sdeg, face_quad_x);
  BernsteinBasis_Array * ba_t = new BernsteinBasis_Array(tdeg, face_quad_y);
  BernsteinBasis_Array * ba_u = new BernsteinBasis_Array(udeg, face_quad_z);

  FEAElement * element = new FEAElement_NURBS_3D_der1_lap_vms_jac( eindex,
      hx, hy, hz, ba_s, ba_t, ba_u, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z,
      eleCtrlPts_w, ext_x, ext_y, ext_z );

  const int face_nqp = face_Int_w->get_num();
  
  int ii, qua, A; // iterator
  double gwts, coor_x, coor_y, coor_z;
  double gn_x, gn_y, gn_z; // cartesian component of Neumann bc
  double nx, ny, nz; // normal vector on the element's top surface in physical domain
  double surface_area;
  const double curr = time + alpha_f * dt; 
  
  Zero_Residual();
  
  for(qua = 0; qua<face_nqp; ++qua)
  {
    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    
    element->get_R(qua, R);
    element->get_3d_normal_top(qua, nx, ny, nz, surface_area);

    for(ii=0; ii<nLocBas; ++ii)
    {
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    gwts = surface_area * hx * hy * face_Int_w->get_weight(qua);
    
    get_top_gn(coor_x, coor_y, coor_z, curr, nx, ny, nz, gn_x, gn_y, gn_z);
 
    for(A=0; A<nLocBas; ++A)
    {
      Residual[4*A]   -= gwts * R[A] * gn_x;
      Residual[4*A+1] -= gwts * R[A] * gn_y;
      Residual[4*A+2] -= gwts * R[A] * gn_z; 
    }
  }

  // clean up the memory allocation
  delete element;
  delete ba_s; delete ba_t; delete ba_u; delete face_Int_w;
  delete face_quad_x; delete face_quad_y; delete face_quad_z;
}



void PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm::Assem_Residual_BotFace(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    const int &eindex,
    const double &hx, const double &hy, const double &hz,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &eleCtrlPts_w,
    const double * const &ext_x,
    const double * const &ext_y,
    const double * const &ext_z
    )
{
  // Build the boundary quadrature info for this element : eindex
  IQuadPts * face_quad_z = new QuadPts_bc0();
  IQuadPts * face_quad_y = new QuadPts_Gauss(nqpy);
  IQuadPts * face_quad_x = new QuadPts_Gauss(nqpx);

  AInt_Weight * face_Int_w = new AInt_Weight(face_quad_x, face_quad_y, face_quad_z);

  BernsteinBasis_Array * ba_s = new BernsteinBasis_Array(sdeg, face_quad_x);
  BernsteinBasis_Array * ba_t = new BernsteinBasis_Array(tdeg, face_quad_y);
  BernsteinBasis_Array * ba_u = new BernsteinBasis_Array(udeg, face_quad_z);

  FEAElement * element = new FEAElement_NURBS_3D_der1_lap_vms_jac( eindex,
      hx, hy, hz, ba_s, ba_t, ba_u, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z,
      eleCtrlPts_w, ext_x, ext_y, ext_z );

  const int face_nqp = face_Int_w->get_num();
  
  int ii, qua, A; // iterator
  double gwts, coor_x, coor_y, coor_z;
  double gn_x, gn_y, gn_z; // cartesian component of Neumann bc
  double nx, ny, nz; // normal vector on the element's top surface in physical domain
  double surface_area; // surface area element, given as the norm of r_s x r_t
  const double curr = time + alpha_f * dt; 
  
  Zero_Residual();
  
  for(qua = 0; qua<face_nqp; ++qua)
  {
    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    
    element->get_R(qua, R);
    element->get_3d_normal_bottom(qua, nx, ny, nz, surface_area);

    for(ii=0; ii<nLocBas; ++ii)
    {
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }
    
    gwts = surface_area * hx * hy * face_Int_w->get_weight(qua);
    
    get_bot_gn(coor_x, coor_y, coor_z, curr, nx, ny, nz, gn_x, gn_y, gn_z);
 
    for(A=0; A<nLocBas; ++A)
    {
      Residual[4*A]   -= gwts * R[A] * gn_x;
      Residual[4*A+1] -= gwts * R[A] * gn_y;
      Residual[4*A+2] -= gwts * R[A] * gn_z; 
    }
  }

  // clean up the memory allocation
  delete element;
  delete ba_s; delete ba_t; delete ba_u; delete face_Int_w;
  delete face_quad_x; delete face_quad_y; delete face_quad_z;
}



void PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &velo,
    const double * const &disp,
    const int &eindex,
    const double &hx, const double &hy, const double &hz,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const BernsteinBasis_Array * const &bu,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &eleCtrlPts_w,
    const double * const &ext_x,
    const double * const &ext_y,
    const double * const &ext_z,
    const class AInt_Weight * const &weight )
{
  // Generate the element basis function value at quad-points
  FEAElement * element = new FEAElement_NURBS_3D_der1_lap_vms( eindex, 
      hx, hy, hz, bs, bt, bu, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z,
      eleCtrlPts_w, ext_x, ext_y, ext_z ); 

  int ii, qua, A; // iterator
  int jj, B, index; // iterator for tangent assembly
  double u, u_t, u_x, u_y, u_z, u_xx, u_yy, u_zz;
  double v, v_t, v_x, v_y, v_z, v_xx, v_yy, v_zz;
  double w, w_t, w_x, w_y, w_z, w_xx, w_yy, w_zz;
  double p, f1, f2, f3;
  //double mf1, mf2, mf3, mf4;
  double p_x, p_y, p_z;
  double gwts, coor_x, coor_y, coor_z;

  // Momentum residual and their derivatives
  double rx, ry, rz;

  // Stabilization parameter
  double tau_m, tau_c;

  // Frequently used quantities in assembly
  const double two_nu = 2.0 * nu;
  double NA, NA_x, NA_y, NA_z;
  double velo_dot_gradR, div_vel;
  double r_dot_gradR; // rx * NA_x + ry * NA_y + rz * NA_z
  double tau_m_2; // tau_m * tau_m

  double NB, NB_x, NB_y, NB_z, NB_xx, NB_yy, NB_zz, NB_Lap;
  double NANB, NAxNB, NAyNB, NAzNB;
  double NANBx, NAxNBx, NAyNBx, NAzNBx;
  double NANBy, NAxNBy, NAyNBy, NAzNBy;
  double NANBz, NAxNBz, NAyNBz, NAzNBz;

  double drx_du_B, drx_dv_B, drx_dw_B, drx_dp_B;
  double dry_du_B, dry_dv_B, dry_dw_B, dry_dp_B;
  double drz_du_B, drz_dv_B, drz_dw_B, drz_dp_B;

  double velo_dot_gradNB;

  const double curr = time + alpha_f * dt;

  Zero_Tangent_Residual();

  // Initialize the Sub_Tan_x matrix 
  for(ii=0; ii<16; ++ii)
  {
    for(jj=0; jj<nLocBas*nLocBas; ++jj)
    {
      Sub_Tan_m[ii][jj] = 0.0;
      Sub_Tan_f[ii][jj] = 0.0;
      Sub_Tan_p[ii][jj] = 0.0;
    }
  }

  for(qua=0; qua<nqp; ++qua)
  {
    u = 0.0; u_t = 0.0; u_x = 0.0; u_y = 0.0; u_z = 0.0;
    u_xx = 0.0; u_yy = 0.0; u_zz = 0.0;
    v = 0.0; v_t = 0.0; v_x = 0.0; v_y = 0.0; v_z = 0.0;
    v_xx = 0.0; v_yy = 0.0; v_zz = 0.0;
    w = 0.0; w_t = 0.0; w_x = 0.0; w_y = 0.0; w_z = 0.0;
    w_xx = 0.0; w_yy = 0.0; w_zz = 0.0;
    p = 0.0; coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    p_x = 0.0; p_y = 0.0; p_z = 0.0;

    element->get_3D_R_gradR_LaplacianR(qua, R, dR_dx, dR_dy, dR_dz,
        d2R_dxx, d2R_dyy, d2R_dzz);
    element->get_invJacobian(qua, dxi_dx); 

    for(ii=0; ii<nLocBas; ++ii)
    {
      u_t += velo[4*ii]   * R[ii];
      v_t += velo[4*ii+1] * R[ii];
      w_t += velo[4*ii+2] * R[ii];

      u += disp[4*ii]   * R[ii];
      v += disp[4*ii+1] * R[ii];
      w += disp[4*ii+2] * R[ii];
      p += disp[4*ii+3] * R[ii];

      u_x += disp[4*ii]   * dR_dx[ii];
      v_x += disp[4*ii+1] * dR_dx[ii];
      w_x += disp[4*ii+2] * dR_dx[ii];
      p_x += disp[4*ii+3] * dR_dx[ii];

      u_y += disp[4*ii]   * dR_dy[ii];
      v_y += disp[4*ii+1] * dR_dy[ii];
      w_y += disp[4*ii+2] * dR_dy[ii];
      p_y += disp[4*ii+3] * dR_dy[ii];

      u_z += disp[4*ii]   * dR_dz[ii];
      v_z += disp[4*ii+1] * dR_dz[ii];
      w_z += disp[4*ii+2] * dR_dz[ii];
      p_z += disp[4*ii+3] * dR_dz[ii];

      u_xx += disp[4*ii]   * d2R_dxx[ii];
      v_xx += disp[4*ii+1] * d2R_dxx[ii];
      w_xx += disp[4*ii+2] * d2R_dxx[ii];

      u_yy += disp[4*ii]   * d2R_dyy[ii];
      v_yy += disp[4*ii+1] * d2R_dyy[ii];
      w_yy += disp[4*ii+2] * d2R_dyy[ii];

      u_zz += disp[4*ii]   * d2R_dzz[ii];
      v_zz += disp[4*ii+1] * d2R_dzz[ii];
      w_zz += disp[4*ii+2] * d2R_dzz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    get_tau(tau_m, tau_c, dt, dxi_dx, u, v, w);

    tau_m_2 = tau_m * tau_m;

    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    get_f(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    //get_mf(coor_x, coor_y, coor_z, curr, mf1, mf2, mf3, mf4);

    rx = u_t + u_x * u + u_y * v + u_z * w + p_x - nu * (u_xx + u_yy + u_zz) - f1; 
    ry = v_t + v_x * u + v_y * v + v_z * w + p_y - nu * (v_xx + v_yy + v_zz) - f2;
    rz = w_t + w_x * u + w_y * v + w_z * w + p_z - nu * (w_xx + w_yy + w_zz) - f3;

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];
      velo_dot_gradR = NA_x * u + NA_y * v + NA_z * w;
      r_dot_gradR = NA_x * rx + NA_y * ry + NA_z * rz;
      div_vel = u_x + v_y + w_z;

      Residual[4*A] += gwts * (NA * u_t 
          + NA * (u * u_x + v * u_y + w * u_z) 
          - NA_x * p
          + NA_x * two_nu * u_x 
          + NA_y * nu * (u_y + v_x) 
          + NA_z * nu * (u_z + w_x)
          + velo_dot_gradR * tau_m * rx 
          + r_dot_gradR * u * tau_m
          + NA_x * tau_c * div_vel 
          - r_dot_gradR * tau_m_2 * rx
          - NA * f1 );

      Residual[4*A+1] += gwts * (NA * v_t 
          + NA * (u * v_x + v * v_y + w * v_z) 
          - NA_y * p
          + NA_x * nu * (u_y + v_x)
          + NA_y * two_nu * v_y
          + NA_z * nu * (v_z + w_y)
          + velo_dot_gradR * tau_m * ry
          + r_dot_gradR * v * tau_m
          + NA_y * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * ry
          - NA * f2 );

      Residual[4*A+2] += gwts * (NA * w_t 
          + NA * (u * w_x + v * w_y + w * w_z) 
          - NA_z * p
          + NA_x * nu * (u_z + w_x)
          + NA_y * nu * (w_y + v_z)
          + NA_z * two_nu * w_z
          + velo_dot_gradR * tau_m * rz
          + r_dot_gradR * w * tau_m 
          + NA_z * tau_c * div_vel
          - r_dot_gradR * tau_m_2 * rz
          - NA * f3 );

      Residual[4*A+3] += gwts * ( NA * div_vel + tau_m * r_dot_gradR );


      for(B=0; B<nLocBas; ++B)
      {
        index = B + A * nLocBas;
        NB = R[B]; NB_x = dR_dx[B]; NB_y = dR_dy[B]; NB_z = dR_dz[B];
        NB_xx = d2R_dxx[B]; NB_yy = d2R_dyy[B]; NB_zz = d2R_dzz[B];
        NB_Lap = NB_xx + NB_yy + NB_zz;

        velo_dot_gradNB = u * NB_x + v * NB_y + w * NB_z;

        NANB  = NA * NB;   NANBx  = NA * NB_x;   NANBy  = NA * NB_y;   NANBz  = NA * NB_z;
        NAxNB = NA_x * NB; NAxNBx = NA_x * NB_x; NAxNBy = NA_x * NB_y; NAxNBz = NA_x * NB_z;
        NAyNB = NA_y * NB; NAyNBx = NA_y * NB_x; NAyNBy = NA_y * NB_y; NAyNBz = NA_y * NB_z;
        NAzNB = NA_z * NB; NAzNBx = NA_z * NB_x; NAzNBy = NA_z * NB_y; NAzNBz = NA_z * NB_z; 

        drx_du_B = u_x * NB + velo_dot_gradNB - nu * NB_Lap;  
        drx_dv_B = u_y * NB;
        drx_dw_B = u_z * NB;
        drx_dp_B = NB_x;

        dry_du_B = v_x * NB;
        dry_dv_B = v_y * NB + velo_dot_gradNB - nu * NB_Lap;
        dry_dw_B = v_z * NB;
        dry_dp_B = NB_y;

        drz_du_B = w_x * NB;
        drz_dv_B = w_y * NB;
        drz_dw_B = w_z * NB + velo_dot_gradNB - nu * NB_Lap;
        drz_dp_B = NB_z;

        // first row of tangent matrix  
        Sub_Tan_m[0][index] += gwts * ( NANB 
            + tau_m * (2.0 * u * NAxNB + v * NAyNB + w * NAzNB)
            - tau_m_2 * (2.0 * rx * NAxNB + ry * NAyNB + rz * NAzNB ) );

        Sub_Tan_f[0][index] += gwts * ( NA * velo_dot_gradNB + NANB * u_x 
            + nu * (2.0*NAxNBx + NAyNBy + NAzNBz) 
            + tau_m * (2.0 * rx * NAxNB + ry * NAyNB + rz * NAzNB)
            + (velo_dot_gradR + u * NA_x) * tau_m * drx_du_B
            + tau_m * u * v_x * NAyNB + tau_m * u * w_x * NAzNB + tau_c * NAxNBx
            - 2.0 * tau_m_2 * rx  * NA_x * drx_du_B
            - tau_m_2 * ry * NA_y * drx_du_B 
            - tau_m_2 * rz * NA_z * drx_du_B
            - tau_m_2 * rx * (v_x * NAyNB + w_x * NAzNB) );

        Sub_Tan_m[1][index] += gwts * ( tau_m * u - tau_m_2 * rx ) * NAyNB;

        Sub_Tan_f[1][index] += gwts * ( NANB * u_y + nu * NAyNBx + tau_m * rx * NAyNB
            + (velo_dot_gradR + u * NA_x) * tau_m * drx_dv_B 
            + tau_m * u * ( NA_y * dry_dv_B + NA_z * drz_dv_B ) + tau_c * NAxNBy
            - 2.0 * tau_m_2 * rx * NA_x * drx_dv_B
            - tau_m_2 * NA_y * (rx * dry_dv_B + ry * drx_dv_B)
            - tau_m_2 * NA_z * (rx * drz_dv_B + rz * drx_dv_B) );

        Sub_Tan_m[2][index] += gwts * ( tau_m * u - tau_m_2 * rx ) * NAzNB;

        Sub_Tan_f[2][index] += gwts * ( NANB * u_z + nu * NAzNBx + tau_m * rx * NAzNB
            + (velo_dot_gradR + u * NA_x) * tau_m * drx_dw_B
            + tau_m * u * ( NA_y * dry_dw_B + NA_z * drz_dw_B ) + tau_c * NAxNBz
            - 2.0 * tau_m_2 * rx * NA_x * drx_dw_B
            - tau_m_2 * NA_y * (rx * dry_dw_B + ry * drx_dw_B)
            - tau_m_2 * NA_z * (rx * drz_dw_B + rz * drx_dw_B) );

        Sub_Tan_p[3][index] += gwts * ((-1.0) * NAxNB 
            + (velo_dot_gradR + u * NA_x) * tau_m * NB_x 
            + tau_m * u * (NAyNBy + NAzNBz)
            - 2.0 * tau_m_2 * rx * NA_x * drx_dp_B
            - tau_m_2 * NA_y * (rx * dry_dp_B + ry * drx_dp_B)
            - tau_m_2 * NA_z * (rx * drz_dp_B + rz * drx_dp_B) );

        // second row of tangent matrix
        Sub_Tan_m[4][index] += gwts * (tau_m * v - tau_m_2 * ry) * NAxNB; 

        Sub_Tan_f[4][index] += gwts * ( NANB * v_x + nu * NAxNBy
            + tau_m * ry * NAxNB 
            + NA_x * v * tau_m * drx_du_B
            + (velo_dot_gradR + NA_y * v) * tau_m * dry_du_B
            + NA_z * v * tau_m * drz_du_B 
            + tau_c * NAyNBx
            - tau_m_2 * NA_x * (ry * drx_du_B + rx * dry_du_B)
            - 2.0 * tau_m_2 * ry * NA_y * dry_du_B
            - tau_m_2 * NA_z * (ry * drz_du_B + rz * dry_du_B) );


        Sub_Tan_m[5][index] += gwts * ( NANB + (velo_dot_gradR + NA_y * v) * tau_m * NB 
            - tau_m_2 * (rx * NAxNB + 2.0 * ry * NAyNB + rz * NAzNB) );

        Sub_Tan_f[5][index] += gwts * ( NA * velo_dot_gradNB + NANB * v_y
            + nu * (NAxNBx + 2.0 * NAyNBy + NAzNBz)
            + tau_m * (rx * NAxNB + 2.0 * ry * NAyNB + rz * NAzNB)
            + tau_m * v * NA_x * drx_dv_B
            + (velo_dot_gradR + v * NA_y) * tau_m * dry_dv_B
            + tau_m * v * NA_z * drz_dv_B
            + tau_c * NAyNBy
            - tau_m_2 * NA_x * (rx * dry_dv_B + ry * drx_dv_B)
            - 2.0 * tau_m_2 * ry * NA_y * dry_dv_B
            - tau_m_2 * NA_z * (ry * drz_dv_B + rz * dry_dv_B) );


        Sub_Tan_m[6][index] += gwts * (tau_m * v - tau_m_2 * ry) * NAzNB;

        Sub_Tan_f[6][index] += gwts * ( NANB * v_z + nu * NAzNBy
            + tau_m * ry * NAzNB
            + NA_x * v * tau_m * drx_dw_B
            + (velo_dot_gradR  + NA_y * v) * tau_m * dry_dw_B
            + NA_z * v * tau_m * drz_dw_B
            + tau_c * NAyNBz
            - tau_m_2 * NA_x * (rx * dry_dw_B + ry * drx_dw_B)
            - tau_m_2 * 2.0 * ry * NA_y * dry_dw_B
            - tau_m_2 * NA_z * (ry * drz_dw_B + rz * dry_dw_B) );

        Sub_Tan_p[7][index] += gwts * ( (-1.0) * NAyNB
            + tau_m * v * (NAxNBx + NAzNBz)
            + (velo_dot_gradR + NA_y * v) * tau_m * NB_y
            - tau_m_2 * NA_x * (rx * NB_y + ry * NB_x)
            - 2.0 * tau_m_2 * ry * NAyNBy
            - tau_m_2 * NA_z * (ry * NB_z + rz * NB_y) );

        // third row of tangent matrix
        Sub_Tan_m[8][index] += gwts * (tau_m * w - tau_m_2 * rz) * NAxNB;

        Sub_Tan_f[8][index] += gwts * ( NANB * w_x + nu * NAxNBz
            + tau_m * rz * NAxNB 
            + tau_m * w * NA_x * drx_du_B 
            + tau_m * w * NA_y * dry_du_B
            + (velo_dot_gradR + w * NA_z) * tau_m * drz_du_B
            + tau_c * NAzNBx
            - tau_m_2 * NA_x * (rx * drz_du_B + rz * drx_du_B)
            - tau_m_2 * NA_y * (ry * drz_du_B + rz * dry_du_B)
            - 2.0 * tau_m_2 * rz * NA_z * drz_du_B );

        Sub_Tan_m[9][index] += gwts * (tau_m * w - tau_m_2 * rz) * NAyNB;

        Sub_Tan_f[9][index] += gwts * ( NANB * w_y + nu * NAyNBz
            + tau_m * rz * NAyNB
            + tau_m * w * (NA_x * drx_dv_B + NA_y * dry_dv_B)
            + (velo_dot_gradR + w * NA_z) * tau_m * drz_dv_B
            + tau_c * NAzNBy
            - tau_m_2 * NA_x * (rx * drz_dv_B + rz * drx_dv_B)
            - tau_m_2 * NA_y * (ry * drz_dv_B + rz * dry_dv_B)
            - 2.0 * tau_m_2 * rz * NA_z * drz_dv_B );

        Sub_Tan_m[10][index] += gwts * ( NANB + (velo_dot_gradR + w * NA_z) * tau_m * NB
            - tau_m_2 * (rx * NAxNB + ry * NAyNB + 2.0 * rz * NAzNB) );


        Sub_Tan_f[10][index] += gwts * ( NA * velo_dot_gradNB + NANB * w_z 
            + nu * (NAxNBx + NAyNBy + 2.0 * NAzNBz)
            + tau_m * (rx * NAxNB + ry * NAyNB + 2.0 * rz * NAzNB)
            + tau_m * w * NA_x * drx_dw_B
            + tau_m * w * NA_y * dry_dw_B
            + tau_m * (velo_dot_gradR + w * NA_z ) * drz_dw_B
            + tau_c * NAzNBz
            - tau_m_2 * NA_x * (rx * drz_dw_B + rz * drx_dw_B)
            - tau_m_2 * NA_y * (ry * drz_dw_B + rz * dry_dw_B)
            - 2.0 * tau_m_2 * NA_z * rz * drz_dw_B );

        Sub_Tan_p[11][index] += gwts * ((-1.0) * NAzNB 
            + w * tau_m * NAxNBx
            + w * tau_m * NAyNBy
            + tau_m * (velo_dot_gradR + w * NA_z) * NB_z
            - tau_m_2 * NA_x * (rx * NB_z + rz * NB_x)
            - tau_m_2 * NA_y * (ry * NB_z + rz * NB_y)
            - 2.0 * tau_m_2 * rz * NAzNBz );

        // fourth row of tangent matrix
        Sub_Tan_m[12][index] += gwts * tau_m * NAxNB;

        Sub_Tan_f[12][index] += gwts * ( NANBx + tau_m * NA_x * drx_du_B
            + tau_m * NA_y * dry_du_B + tau_m * NA_z * drz_du_B );

        Sub_Tan_m[13][index] += gwts * tau_m * NAyNB;

        Sub_Tan_f[13][index] += gwts * ( NANBy + tau_m * NA_x * drx_dv_B
            + tau_m * NA_y * dry_dv_B + tau_m * NA_z * drz_dv_B );

        Sub_Tan_m[14][index] += gwts * tau_m * NAzNB;

        Sub_Tan_f[14][index] += gwts * ( NANBz + tau_m * NA_x * drx_dw_B
            + tau_m * NA_y * dry_dw_B + tau_m * NA_z * drz_dw_B );

        Sub_Tan_p[15][index] += gwts * tau_m * (NAxNBx + NAyNBy + NAzNBz);


      } // B loop
    } // A loop
  } // qua loop

  for(ii=0; ii<4; ++ii)
  {
    for(jj=0; jj<4; ++jj)
    {
      for(A=0; A<nLocBas; ++A)
      {
        for(B=0; B<nLocBas; ++B)
        {
          Tangent[ 4*nLocBas*(4*A+ii) + 4*B + jj  ] =
            alpha_m * Sub_Tan_m[ii*4+jj][A*nLocBas + B]
            + alpha_f * gamma * dt * Sub_Tan_f[ii*4+jj][A*nLocBas + B]
            + Sub_Tan_p[ii*4+jj][A*nLocBas + B];
        }
      }
    }
  }

  delete element;
}


void PLocAssem_NS_3D_VMS_GenAlpha_ExactJac_noCache_AdvForm::Assem_Mass_Residual(
    const double * const &disp,
    const int &eindex,
    const double &hx, const double &hy, const double &hz,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const BernsteinBasis_Array * const &bu,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &eleCtrlPts_w,
    const double * const &ext_x,
    const double * const &ext_y,
    const double * const &ext_z,
    const class AInt_Weight * const &weight )
{
  // Generate the element basis function value at quad-points
  FEAElement * element = new FEAElement_NURBS_3D_der1_lap_vms( eindex, 
      hx, hy, hz, bs, bt, bu, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z,
      eleCtrlPts_w, ext_x, ext_y, ext_z ); 

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

  for(ii=0; ii<16; ++ii)
  {
    for(jj=0; jj<nLocBas * nLocBas; ++jj)
      Sub_Tan_m[ii][jj] = 0.0;
  }


  for(qua=0; qua<nqp; ++qua)
  {
    u = 0.0; u_x = 0.0; u_y = 0.0; u_z = 0.0;
    v = 0.0; v_x = 0.0; v_y = 0.0; v_z = 0.0;
    w = 0.0; w_x = 0.0; w_y = 0.0; w_z = 0.0;
    p = 0.0; coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    element->get_3D_R_gradR_LaplacianR(qua, R, dR_dx, dR_dy, dR_dz,
        d2R_dxx, d2R_dyy, d2R_dzz);

    for(ii=0; ii<nLocBas; ++ii)
    {
      u += disp[4*ii]   * R[ii];
      v += disp[4*ii+1] * R[ii];
      w += disp[4*ii+2] * R[ii];
      p += disp[4*ii+3] * R[ii];

      u_x += disp[4*ii]   * dR_dx[ii];
      v_x += disp[4*ii+1] * dR_dx[ii];
      w_x += disp[4*ii+2] * dR_dx[ii];

      u_y += disp[4*ii]   * dR_dy[ii];
      v_y += disp[4*ii+1] * dR_dy[ii];
      w_y += disp[4*ii+2] * dR_dy[ii]; 

      u_z += disp[4*ii]   * dR_dz[ii];
      v_z += disp[4*ii+1] * dR_dz[ii];
      w_z += disp[4*ii+2] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii]; 
    }

    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    get_f(coor_x, coor_y, coor_z, curr, f1, f2, f3);

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      Residual[4*A] += gwts * ( NA*(u*u_x + v*u_y + w*u_z) 
          - NA_x * p
          + two_nu * NA_x * u_x 
          + nu * NA_y * (u_y + v_x) 
          + nu * NA_z * (u_z + w_x)
          - NA * f1 );

      Residual[4*A+1] += gwts * ( NA*(u*v_x + v*v_y + w*v_z) 
          - NA_y * p
          + nu * NA_x * (u_y + v_x)
          + two_nu * NA_y * v_y
          + nu * NA_z * (v_z + w_y)
          - NA * f2 );

      Residual[4*A+2] += gwts * ( NA*(u*w_x + v*w_y + w*w_z) 
          - NA_z * p 
          + nu * NA_x * (u_z + w_x)
          + nu * NA_y * (w_y + v_z)
          + two_nu * NA_z * w_z
          - NA * f3 );

      for(B=0; B<nLocBas; ++B)
      {
        index = A * nLocBas + B;

        Sub_Tan_m[0][index]  += gwts * NA * R[B];
        Sub_Tan_m[5][index]  += gwts * NA * R[B];
        Sub_Tan_m[10][index] += gwts * NA * R[B];
        Sub_Tan_m[15][index] += gwts * NA * R[B];
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
            Sub_Tan_m[ii*4+jj][A*nLocBas + B];
        }
      }
    }
  }

  delete element;
}

// EOF
