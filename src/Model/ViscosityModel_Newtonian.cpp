#include "ViscosityModel_Newtonian.hpp"

ViscosityModel_Newtonian::ViscosityModel_Newtonian( const double &in_mu )
: mu( in_mu )
{}

void ViscosityModel_Newtonian::print_info() const
{
  SYS_T::commPrint("\t  ViscosityModel_Newtonian:: \n");
  SYS_T::commPrint("\t  Viscosity mu = %e \n", mu);
}
