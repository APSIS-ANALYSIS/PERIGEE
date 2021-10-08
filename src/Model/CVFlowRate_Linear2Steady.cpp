#include "CVFlowRate_Linear2Steady.hpp"

CVFlowRate_Linear2Steady::CVFlowRate_Linear2Steady(
    const ALocal_Inflow_NodalBC * const &infnbc_part,
    const double &in_thred_time, const double &flrate )
: num_nbc(infnbc_part->get_num_nbc()), thred_time(in_thred_time), target_flow_rate(flrate)
{
  // Calculate flow rate and record in txt file 
  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    if( SYS_T::get_MPI_rank() == 0 )
    {
      std::ofstream ofile;
      ofile.open( gen_flowfile_name(nbc_id).c_str(), std::ofstream::out | std::ofstream::trunc );
      for( double tt = 0; tt <= thred_time * 2.0; tt += 0.001 )
        ofile<<tt<<'\t'<<get_flow_rate(nbc_id, tt)<<'\n';
      ofile.close();
    }
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

CVFlowRate_Linear2Steady::~CVFlowRate_Linear2Steady()
{}

double CVFlowRate_Linear2Steady::get_flow_rate( const int &nbc_id,
    const double &time ) const
{
  double out_rate = target_flow_rate;

  if( time < thred_time && time >= 0.0 ) 
    out_rate = target_flow_rate * time / thred_time;

  return out_rate;
}

void CVFlowRate_Linear2Steady::print_info() const
{
  SYS_T::commPrint("     CVFlowRate_Linear2Steady:\n");
  SYS_T::commPrint("     Steady flow rate is %e \n", target_flow_rate);
  SYS_T::commPrint("     Time to reach steady state is %e \n", thred_time);
}

// EOF
