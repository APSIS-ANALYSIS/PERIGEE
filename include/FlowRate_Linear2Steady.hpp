#ifndef FLOWRATE_LINEAR2STEADY_HPP
#define FLOWRATE_LINEAR2STEADY_HPP
// ============================================================================
// FlowRate_Linear2Steady.hpp
//
// Linear function starting from zero flow rate to a prescribed flow rate at a
// prescribed time.
//
// Author: Ju Liu
// Date Created: Oct. 1 2017
// ============================================================================
#include <string>
#include <vector>
#include "IFlowRate.hpp"

class FlowRate_Linear2Steady final : public IFlowRate
{
  public:
    // This constructor will set a uniform thred_time for all inlets, and the
    // target flow rate are determined from the file and set to be the flow rate
    // at the initial time (i.e. time = 0.0) by setting target_flow_rate to be
    // the sum of coef_a.
    FlowRate_Linear2Steady( const std::string &filename );

    ~FlowRate_Linear2Steady() override = default;

    // From time 0 to in_thred_time, the flow rate equals
    // start_rate + (target_rate - start_rate) * time / in_thred_time
    // From in_thred_time to infty, flow_rate = target_rate
    double get_flow_rate( int nbc_id, double time ) const override;

    // Get the turbulance intensity
    double get_flow_TI_std_dev( int nbc_id ) const override { return TI_std_dev[nbc_id]; }

    // Get the start rate if the users want to now it
    double get_flow_start_rate( int nbc_id ) const override { return start_flow_rate[nbc_id]; }

    int get_num_nbc() const override { return num_nbc; }

    void print_info() const override;

  private:
    int num_nbc;

    std::vector<double> thred_time;
    
    std::vector<double> target_flow_rate;

    std::vector<double> start_flow_rate;

    std::vector<double> TI_std_dev;
};

#endif
