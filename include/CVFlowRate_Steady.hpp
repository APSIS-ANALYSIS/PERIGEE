#ifndef CVFLOWRATE_STEADY_HPP
#define CVFLOWRATE_STEADY_HPP
// ============================================================================
// CVFlowRate_Steady.hpp
//
// Steady flow for inflow condition
//
// Author: Ju Liu
// Date Created: May 1 2021
// ============================================================================
#include "Sys_Tools.hpp"
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "ICVFlowRate.hpp"

class CVFlowRate_Steady : public ICVFlowRate
{
  public:
    CVFlowRate_Steady( const char * const &filename );

    CVFlowRate_Steady( const double &in_flowrate );

    virtual ~CVFlowRate_Steady();

    virtual double get_flow_rate(const double &time) const;

    virtual void print_info() const;

  private:
    double flowrate;
};

#endif
