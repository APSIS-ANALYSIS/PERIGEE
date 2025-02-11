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
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "IFlowRate.hpp"

class FlowRate_Linear2Steady : public IFlowRate
{
  public:
    // This constructor will set a uniform thred_time for all inlets, and the
    // target flow rate are determined from the file and set to be the flow rate
    // at the initial time (i.e. time = 0.0) by setting target_flow_rate to be
    // the sum of coef_a.
    FlowRate_Linear2Steady( const std::string &filename );

    virtual ~FlowRate_Linear2Steady() = default;

    // From time 0 to in_thred_time, the flow rate equals
    // start_rate + (target_rate - start_rate) * time / in_thred_time
    // From in_thred_time to infty, flow_rate = target_rate
    virtual double get_flow_rate( const int &nbc_id, const double &time ) const;

    // Get the turbulance intensity
    virtual double get_flow_TI_std_dev( const int &nbc_id ) const { return TI_std_dev[nbc_id]; }

    // Get the start rate if the users want to now it
    virtual double get_flow_start_rate( const int &nbc_id ) const { return start_flow_rate[nbc_id]; }

    virtual int get_num_nbc() const { return num_nbc; }

    virtual void print_info() const;

  private:
    int num_nbc;

    std::vector<double> thred_time;
    
    std::vector<double> target_flow_rate;

    std::vector<double> start_flow_rate;

    std::vector<double> TI_std_dev;
};

#endif
