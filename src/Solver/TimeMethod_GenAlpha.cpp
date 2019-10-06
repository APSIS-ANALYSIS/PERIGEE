#include "TimeMethod_GenAlpha.hpp"

TimeMethod_GenAlpha::TimeMethod_GenAlpha(const double &input_spectral )
: rho_infty( input_spectral ), 
  alpha_m( 0.5 * (3.0 - input_spectral) / (1.0 + input_spectral) ),
  alpha_f( 1.0/(1.0 + input_spectral) ),
  gamma( 0.5 + alpha_m - alpha_f ),
  beta( 0.25 * (1.0 - alpha_f + alpha_m) * (1.0 - alpha_f + alpha_m) ), 
  is2nd(false), is_rho_set(true)
{
}


TimeMethod_GenAlpha::TimeMethod_GenAlpha( const double &in_alpha_m,
    const double &in_alpha_f, const double &in_gamma )
: rho_infty(0.0), alpha_m(in_alpha_m), 
  alpha_f(in_alpha_f), gamma(in_gamma), 
  beta( 0.25 * (1.0 - alpha_f + alpha_m) * (1.0 - alpha_f + alpha_m) ), 
  is2nd(false), is_rho_set(false)
{
}


TimeMethod_GenAlpha::~TimeMethod_GenAlpha()
{}


TimeMethod_GenAlpha::TimeMethod_GenAlpha( const double &input_spectral,
    const bool &input_is2ndorder )
: rho_infty(input_spectral), is2nd( input_is2ndorder ), is_rho_set(true)
{
  if(is2nd)
  {
    alpha_m = (2.0-rho_infty) / (1.0+rho_infty);
    alpha_f = 1.0 / (1.0+rho_infty);
  }
  else
  {
    alpha_m = 0.5 * (3.0 - rho_infty) / (1.0 + rho_infty);
    alpha_f = 1.0 / (1.0 + rho_infty);
  }
  
  gamma = 0.5 - alpha_f + alpha_m;
  beta = 0.25 * (1.0 - alpha_f + alpha_m) * (1.0 - alpha_f + alpha_m);
}



void TimeMethod_GenAlpha::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "Generalized-alpha method ");
  if(is2nd) PetscPrintf(PETSC_COMM_WORLD, "for 2nd-order system: \n");
  else PetscPrintf(PETSC_COMM_WORLD, "for 1st-order system: \n");
  if(is_rho_set) PetscPrintf(PETSC_COMM_WORLD, "  --- rho_inf: %e \n", rho_infty);
  else PetscPrintf(PETSC_COMM_WORLD, "  --- rho_inf: unset \n");
  PetscPrintf(PETSC_COMM_WORLD, "  --- Alpha_f: %e \n", alpha_f);
  PetscPrintf(PETSC_COMM_WORLD, "  --- Alpha_m: %e \n", alpha_m);
  PetscPrintf(PETSC_COMM_WORLD, "  --- Gamma  : %e \n", gamma);
  PetscPrintf(PETSC_COMM_WORLD, "  --- Beta   : %e \n", beta);
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}

// EOF
