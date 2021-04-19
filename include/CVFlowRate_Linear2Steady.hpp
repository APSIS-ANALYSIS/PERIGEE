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

class CVFlowRate_Linear2Steady : public ICVFlowRate
{
  public:
    // From time 0 to in_time, the flow rate = time * flrate / in_time
    // From in_time to infty, flow_rate = flrate
    CVFlowRate_Linear2Steady(const double &in_time,
        const double &flrate,
        const bool &prestress_flag = false);

    virtual ~CVFlowRate_Linear2Steady();

    virtual double get_flow_rate(const double &time) const;

    virtual void print_info() const;

  private:
    double thred_time;
    const double target_flow_rate;
};

#endif
