#ifndef FLOWRATE_COSINE2STEADY_HPP
#define FLOWRATE_COSINE2STEADY_HPP
// ============================================================================
// FlowRate_Cosine2Steady.hpp
//
// It is similar to FlowRate_Linear2Steady, except using a cosine function
// to prescibe the starting flow rate
//
// Date Created: Nov. 17 2023
// ============================================================================
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "IFlowRate.hpp"

class FlowRate_Cosine2Steady : public IFlowRate
{
  public:
    // This constructor will set a uniform thred_time for all inlets, and the
    // target flow rate are determined from the file and set to be the flow rate
    // at the initial time (i.e. time = 0.0) by setting target_flow_rate to be
    // the sum of coef_a.
    // Same as FlowRate_Linear2Steady.
    FlowRate_Cosine2Steady( const double &in_thred_time, const std::string &filename );

    virtual ~FlowRate_Cosine2Steady() = default;

    // From time 0 to in_time, the flow rate equals
    // 0.5 * flrate * (1 -  cos (PI * time / in_thred_time))
    // From in_time to infty, flow_rate = flrate
    virtual double get_flow_rate( const int &nbc_id, const double &time ) const;

    virtual double get_flow_TI_std_dev( const int &nbc_id ) const { return TI_std_dev[nbc_id]; }

    virtual double get_start_flow_rate( const int &nbc_id ) const { return start_flow_rate[nbc_id]; }

    virtual int get_num_nbc() const { return num_nbc; }

    virtual void print_info() const;

  private:
    int num_nbc;
    
    const double thred_time;
    
    std::vector<double> target_flow_rate;

    std::vector<double> start_flow_rate;

    std::vector<double> TI_std_dev;
};

#endif
