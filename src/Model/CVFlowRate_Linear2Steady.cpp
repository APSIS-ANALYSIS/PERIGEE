#include "CVFlowRate_Linear2Steady.hpp"

CVFlowRate_Linear2Steady::CVFlowRate_Linear2Steady(
    const double &in_time, const double &flrate)
: thred_time(in_time), target_flow_rate(flrate)
{
  // Calculate flow rate and record in txt file 
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if( rank == 0 )
  {
    std::ofstream ofile;
    ofile.open( "Inlet_flowrate.txt", std::ofstream::out | std::ofstream::trunc );
    for(double tt=0; tt <= thred_time * 2.0; tt += 0.001 )
      ofile<<tt<<'\t'<<get_flow_rate(tt)<<'\n';
    ofile.close();
  }
}

CVFlowRate_Linear2Steady::~CVFlowRate_Linear2Steady()
{}

double CVFlowRate_Linear2Steady::get_flow_rate(const double &time) const
{
  double out_rate = 0.0;
  if( time < thred_time && time >= 0.0 ) 
    out_rate = target_flow_rate * time / thred_time;
  else 
    out_rate = target_flow_rate;

  return out_rate;
}

void CVFlowRate_Linear2Steady::print_info() const
{
  SYS_T::commPrint("---- CVFlowRate_Linear2Steady:\n");
  SYS_T::commPrint("  Steady flow rate is %e \n", target_flow_rate);
  SYS_T::commPrint("  Time to steady state is %e \n", thred_time);
}

// EOF
