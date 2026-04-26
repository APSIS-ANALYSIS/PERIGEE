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
#include <vector>
#include "IFlowRate.hpp"

class FlowRate_Steady final : public IFlowRate
{
  public:
    FlowRate_Steady( const std::string &filename );

    ~FlowRate_Steady() override = default;

    double get_flow_rate(int nbc_id, double time) const override;

    int get_num_nbc() const override { return num_nbc; }

    void print_info() const override;

  private:
    int num_nbc;
    
    std::vector<double> flowrate;
};

#endif
