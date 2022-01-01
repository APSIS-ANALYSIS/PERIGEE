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


// EOF
