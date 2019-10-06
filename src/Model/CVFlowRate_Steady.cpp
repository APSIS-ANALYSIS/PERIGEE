#include "CVFlowRate_Steady.hpp"

CVFlowRate_Steady::CVFlowRate_Steady(const double &in_flrate)
: flrate( in_flrate )
{}

CVFlowRate_Steady::~CVFlowRate_Steady()
{}

void CVFlowRate_Steady::print_info() const
{
  SYS_T::commPrint("---- CVFlowRate_steady: ");
  PetscPrintf(PETSC_COMM_WORLD, "steady flow rate is %e \n", flrate);
}

// EOF
