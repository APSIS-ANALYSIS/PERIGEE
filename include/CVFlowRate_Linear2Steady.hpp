#ifndef CVFLOWRATE_LINEAR2STEADY_HPP
#define CVFLOWRATE_LINEAR2STEADY_HPP
// ==================================================================
// CVFlowRate_Linear2Steady.hpp
//
// Linear function starting from zero flow rate to a prescribed flow
// rate at a prescribed time.
//
// Author: Ju Liu
// Date Created: Oct. 1 2017
// ==================================================================
#include "Sys_Tools.hpp"
#include "ICVFlowRate.hpp"
#include "ALocal_Inflow_NodalBC.hpp"

class CVFlowRate_Linear2Steady : public ICVFlowRate
{
  public:
    // From time 0 to in_time, the flow rate = time * flrate / in_thred_time
    // From in_time to infty, flow_rate = flrate
    CVFlowRate_Linear2Steady( const int &input_num_nbc,
        const double &in_thred_time, const double &flrate );

    virtual ~CVFlowRate_Linear2Steady();

    // nbc_id is unused. The same inflow profile is prescribed for all inlets.
    virtual double get_flow_rate( const int &nbc_id, const double &time ) const;

    virtual void print_info() const;

  private:
    const int num_nbc;

    const double thred_time, target_flow_rate;

    // ------------------------------------------------------------------------
    // Generate a filename for inlet face nbc_id as Inlet_xxx_flowrate.txt
    // ------------------------------------------------------------------------
    virtual std::string gen_flowfile_name(const int &nbc_id) const
    {
      std::ostringstream ss;
      ss << "Inlet_";
      if( nbc_id/10 == 0 ) ss << "00";
      else if( nbc_id/100 == 0 ) ss << "0";

      ss << nbc_id << "_flowrate.txt";

      return ss.str();
    }
};

#endif
