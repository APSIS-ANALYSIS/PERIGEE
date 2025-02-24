#ifndef FLOWRATE_STEADY_HPP
#define FLOWRATE_STEADY_HPP
// ============================================================================
// FlowRate_Steady.hpp
//
// Steady flow for inflow condition
//
// Author: Ju Liu
// Date Created: May 1 2021
// ============================================================================
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "IFlowRate.hpp"

class FlowRate_Steady : public IFlowRate
{
  public:
    FlowRate_Steady( const std::string &filename );

    virtual ~FlowRate_Steady() = default;

    virtual double get_flow_rate(const int &nbc_id, const double &time) const;

    virtual int get_num_nbc() const { return num_nbc; }

    virtual void print_info() const;

  private:
    int num_nbc;
    
    std::vector<double> flowrate;
};

#endif
