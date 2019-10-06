#ifndef CVFLOWRATE_STEADY_HPP
#define CVFLOWRATE_STEADY_HPP
// ==================================================================
// CVFlowRate_Steady.hpp
//
// Steady flow for inflow condition.
//
// flrate : flow rate
//
// Author: Ju Liu
// Date Created: Aug. 6 2017
// ==================================================================
#include "Sys_Tools.hpp"
#include "ICVFlowRate.hpp"

class CVFlowRate_Steady : public ICVFlowRate
{
  public:
    CVFlowRate_Steady(const double &in_flrate);

    virtual ~CVFlowRate_Steady();

    virtual double get_flow_rate(const double &time) const
    {return flrate;}

    virtual void print_info() const;

  private:
    const double flrate;
};

#endif
