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

// EOF
