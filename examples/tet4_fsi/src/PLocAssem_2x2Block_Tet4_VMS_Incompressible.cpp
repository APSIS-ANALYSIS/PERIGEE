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


// EOF
