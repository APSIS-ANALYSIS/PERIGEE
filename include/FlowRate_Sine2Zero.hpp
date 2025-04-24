#ifndef FLOWRATE_SINE2ZERO_HPP
#define FLOWRATE_SINE2ZERO_HPP
// ============================================================================
// FlowRate_Sine2Zero.hpp
//
// It is designed to match FlowRate_Cosine2Steady, using a sine function
// to prescibe the starting dot flow rate.
//
// Date Created: Nov. 04 2024
// ============================================================================
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "IFlowRate.hpp"

class FlowRate_Sine2Zero : public IFlowRate
{
  public:
    // This constructor will set a uniform thred_time for all inlets, and the
    // target flow rate are determined from the file and set to be the flow rate
    // at the initial time (i.e. time = 0.0) by setting target_flow_rate to be
    // the sum of coef_a.
    // Same as FlowRate_Cosine2Steady.
    FlowRate_Sine2Zero( const std::string &filename );

    virtual ~FlowRate_Sine2Zero() = default;

    // From time 0 to in_thred_time,
    // flow_rate =  start_rate + 0.5 * (target_rate - start_rate) 
    //                               * (1 -  cos (PI * time / in_thred_time)),
    // dot_flow_rate = 0.5 * (target_rate - start_rate) 
    //                 * PI / in_thred_time * sin (PI * time / in_thred_time),
    // From in_thred_time to infty, flow_rate = target_rate, dot_flow_rate = 0.0.
    virtual double get_flow_rate( const int &nbc_id, const double &time ) const;

    // Get the turbulance intensity
    virtual double get_flow_TI_std_dev( const int &nbc_id ) const 
    { return TI_std_dev[nbc_id]; }

    // Get the start rate if the users want to now it
    virtual double get_flow_start_rate( const int &nbc_id ) const 
    { return start_flow_rate[nbc_id]; }

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
