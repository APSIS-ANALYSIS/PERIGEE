#include "PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha.hpp"

PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha::PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha(
    const TimeMethod_GenAlpha * const &tm_gAlpha, const int &in_nqp, 
    const double &in_rho, const double &in_vis_mu, const double &in_beta )
: rho0( in_rho ), vis_mu( in_vis_mu ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), beta(in_beta), CI(36.0), CT(4.0),
  nLocBas_v(4), nLocBas_p(4), snLocBas_v(3), snLocBas_p(3),
  vec_size_0( nLocBas_v * 3 ), vec_size_1( nLocBas_p ), 
  sur_size_0( snLocBas_v * 3 ), sur_size_1( snLocBas_p ), nqp( in_nqp )
{
  Tangent00 = new PetscScalar[vec_size_0 * vec_size_0];
  Tangent01 = new PetscScalar[vec_size_0 * vec_size_1];
  Tangent10 = new PetscScalar[vec_size_1 * vec_size_0];
  Tangent11 = new PetscScalar[vec_size_1 * vec_size_1];

  Residual0 = new PetscScalar[vec_size_0];
  Residual1 = new PetscScalar[vec_size_1];
  
  sur_Tangent00 = new PetscScalar[sur_size_0 * sur_size_0];
  sur_Tangent01 = new PetscScalar[sur_size_0 * sur_size_1];
  sur_Tangent10 = new PetscScalar[sur_size_1 * sur_size_0];
  sur_Tangent11 = new PetscScalar[sur_size_1 * sur_size_1];

  sur_Residual0 = new PetscScalar[sur_size_0];
  sur_Residual1 = new PetscScalar[sur_size_1];

  Zero_Tangent_Residual();
  Zero_sur_Tangent_Residual();

  // print info of this assembly routine
  print_info();
}


PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha::~PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha()
{}


void PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint("  Three-dimensional Incompressible Navier-Stokes equations: \n");
  SYS_T::commPrint("  FEM: 4-node Tetrahedral \n");
  SYS_T::commPrint("  Spatial: ALE-VMS \n");
  SYS_T::commPrint("  Temporal: Generalized-alpha Method \n");
  SYS_T::commPrint("  Density rho = %e \n", rho0);
  SYS_T::commPrint("  Dynamic Viscosity mu = %e \n", vis_mu);
  SYS_T::commPrint("  Kienmatic Viscosity nu = %e \n", vis_mu / rho0);
  SYS_T::commPrint("  Stabilization para CI = %e \n", CI);
  SYS_T::commPrint("  Stabilization para CT = %e \n", CT);
  SYS_T::commPrint("  Backflow Stab. para beta = %e \n", beta);
  SYS_T::commPrint("  Consistent tangent matrix used. \n");
  SYS_T::commPrint("  Nonlinear quadratic term is in advective form. \n");
  SYS_T::commPrint("  Pressure is evaluated at n+alpha_f rather than n+1. \n");
  SYS_T::commPrint("  Input solution vector should have 7 dofs including the mesh motion\n      information in the first three slots. \n");
  SYS_T::commPrint("  Density rho is in front of the du/dt and the dimension of the\n      equations is momentum time rate. \n");
  SYS_T::print_sep_line();
}

void PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha::get_metric(
    const double * const &f,
    double &G11, double &G12, double &G13,
    double &G22, double &G23, double &G33 ) const
{
  // PHASTA definition 
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

void PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha::get_tau(
    double &tau_m_qua, double &tau_c_qua,
    const double &dt, const double * const &dxi_dx,
    const double &u, const double &v, const double &w ) const
{
  // Use K matrix to correct the metric
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

  tau_c_qua = 1.0 / denom_c;
}


double PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha::get_DC(
    const double * const &dxi_dx,
    const double &u, const double &v, const double &w ) const
{
  //double G11, G12, G13, G22, G23, G33;
  //get_metric( dxi_dx, G11, G12, G13, G22, G23, G33 );

  //dc_tau = G11 * u * u + 2.0 * G12 * u * v + 2.0 * G13 * u * w + G22 * v * v
  //  + 2.0 * G23 * v * w + G33 * w * w;

  //if(dc_tau > 1.0e-15) dc_tau = rho0 * std::pow(dc_tau, -0.5);
  //else dc_tau = 0.0;
  // return dc_tau;

  return 0.0;
}



// EOF
