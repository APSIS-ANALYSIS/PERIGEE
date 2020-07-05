#include "PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha.hpp"

PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha(
    const TimeMethod_GenAlpha * const &tm_gAlpha,
    const int &in_nlocbas, const int &in_nqp,
    const int &in_snlocbas, const double &in_rho,
    const double &in_vis_mu, const double &in_beta,
    const double &in_ctauc, const int &elemtype )
: rho0( in_rho ), vis_mu( in_vis_mu ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), beta(in_beta), nqp(in_nqp),
  Ctauc( in_ctauc )
{
  if(elemtype == 501)
  {
    // 501 is linear element
    CI = 36.0; CT = 4.0;
    nLocBas = 4; snLocBas = 3;
  }
  else if(elemtype == 502)
  {
    // 502 is quadratic element
    CI = 60.0; CT = 4.0;
    nLocBas = 10; snLocBas = 6;
  }
  else SYS_T::print_fatal("Error: unknown elem type.\n");

  vec_size_v = nLocBas * 3;
  vec_size_p = nLocBas;
  sur_size_v = snLocBas * 3;

  R.resize(nLocBas);
  dR_dx.resize(nLocBas); dR_dy.resize(nLocBas); dR_dz.resize(nLocBas);
  d2R_dxx.resize(nLocBas); d2R_dyy.resize(nLocBas); d2R_dzz.resize(nLocBas);

  Tangent00 = new PetscScalar[vec_size_v * vec_size_v];
  Tangent01 = new PetscScalar[vec_size_v * vec_size_p];
  Tangent10 = new PetscScalar[vec_size_p * vec_size_v];
  Tangent11 = new PetscScalar[vec_size_p * vec_size_p];
  
  Residual0 = new PetscScalar[vec_size_v];
  Residual1 = new PetscScalar[vec_size_p];

  sur_Tangent00 = new PetscScalar[sur_size_v * sur_size_v];
  sur_Residual0 = new PetscScalar[sur_size_v];

  Zero_Tangent_Residual();

  Zero_sur_Tangent_Residual();

  print_info();
}


PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::~PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha()
{
  delete Tangent00; Tangent00 = nullptr;
  delete Tangent01; Tangent01 = nullptr;
  delete Tangent10; Tangent10 = nullptr;
  delete Tangent11; Tangent11 = nullptr;

  delete sur_Tangent00; sur_Tangent00 = nullptr;

  delete Residual0; Residual0 = nullptr;
  delete Residual1; Residual1 = nullptr;
  
  delete sur_Residual0; sur_Residual0 = nullptr;
}


void PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("  Three-dimensional Incompressible Navier-Stokes equations: \n");
  if(nLocBas == 4)
    SYS_T::commPrint("  FEM: 4-node Tetrahedral element \n");
  else if(nLocBas == 10)
    SYS_T::commPrint("  FEM: 10-node Tetrahedral element \n");
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
  SYS_T::commPrint("  Note: \n");
  SYS_T::commPrint("  1. Consistent tangent matrix used. \n");
  SYS_T::commPrint("  2. Nonlinear quadratic term is in advective form. \n");
  SYS_T::commPrint("  3. Pressure is evaluated at n+alpha_f rather than n+1. \n");
  SYS_T::commPrint("----------------------------------------------------------- \n");
}


void PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::get_metric( const double * const &f,
    double &G11, double &G12, double &G13,
    double &G22, double &G23, double &G33 ) const
{
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


void PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::get_tau( double &tau_m_qua, double &tau_c_qua,
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

  const double temp_nu = vis_mu / rho0;

  const double denom_m = CT / (dt*dt) + uGu + CI * temp_nu * temp_nu * GdG;

  tau_m_qua = 1.0 / ( rho0 * sqrt(denom_m) );

  const double denom_c = tau_m_qua * g_dot_g;

  tau_c_qua = Ctauc / denom_c;
}


void PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::get_DC( 
    double &dc_tau, const double * const &dxi_dx,
    const double &u, const double &v, const double &w ) const
{
  //double G11, G12, G13, G22, G23, G33;
  //get_metric( dxi_dx, G11, G12, G13, G22, G23, G33 );

  //dc_tau = G11 * u * u + 2.0 * G12 * u * v + 2.0 * G13 * u * w + G22 * v * v
  //  + 2.0 * G23 * v * w + G33 * w * w;

  //if(dc_tau > 1.0e-15) dc_tau = rho0 * std::pow(dc_tau, -0.5);
  //else dc_tau = 0.0;

  dc_tau = 0.0;
}


void PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::Assem_Residual(
    const double &time, const double &dt,
    const double * const &dot_velo,
    const double * const &dot_pres,
    const double * const &velo,
    const double * const &pres,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{}


void PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::Assem_Tangent_Residual( 
    const double &time, const double &dt,
    const double * const &dot_velo,
    const double * const &dot_pres,
    const double * const &velo,
    const double * const &pres,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{}


void PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::Assem_Mass_Residual(
    const double * const &velo,
    const double * const &pres,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{}


void PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::Assem_Residual_EBC(
    const int &ebc_id,
    const double &time, const double &dt,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{}


void PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::Assem_Residual_EBC_Resistance(
    const int &ebc_id,
    const double &val,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{}


void PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::Assem_Residual_BackFlowStab(
    const double * const &velo,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{}


void PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::Assem_Tangent_Residual_BackFlowStab(
    const double &dt,
    const double * const &velo,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{}


double PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::get_flowrate( const double * const &velo,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  double nx, ny, nz, surface_area;

  double flrate = 0.0;

  for(int qua =0; qua< face_nqp; ++qua)
  {
    element->get_R(qua, &R[0]);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);

    double u = 0.0, v = 0.0, w = 0.0;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      u += velo[ii*3]   * R[ii];
      v += velo[ii*3+1] * R[ii];
      w += velo[ii*3+2] * R[ii];
    }

    flrate += surface_area * quad->get_qw(qua) * ( u * nx + v * ny + w * nz );
  }

  return flrate;
}


void PLocAssem_2x2Block_Tet_VMS_NS_GenAlpha::get_pressure_area( 
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad,
    double &pres, double &area )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  double nx, ny, nz, surface_area;

  pres = 0.0;
  area = 0.0;

  for(int qua =0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, &R[0]);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);

    double pp = 0.0;
    for(int ii=0; ii<snLocBas; ++ii) pp += sol[ii] * R[ii];

    pres += surface_area * quad->get_qw(qua) * pp;
    area += surface_area * quad->get_qw(qua);
  }
}

// EOF
