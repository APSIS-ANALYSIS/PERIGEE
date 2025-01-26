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
    // From time 0 to in_time, the flow rate = time * flrate / in_thred_time
    // From in_time to infty, flow_rate = flrate
    // This constructor will set a uniform thred_time and flrate for all inlets.
    FlowRate_Linear2Steady( const int &input_num_nbc,
        const double &in_thred_time, const double &flrate );

    // This constructor will set a uniform thred_time for all inlets, and the
    // target flow rate are determined from the file and set to be the flow rate
    // at the initial time (i.e. time = 0.0) by setting target_flow_rate to be
    // the sum of coef_a.
    FlowRate_Linear2Steady( const double &in_thred_time, 
        const std::string &filename );

    virtual ~FlowRate_Linear2Steady() = default;

    virtual double get_flow_rate( const int &nbc_id, const double &time ) const;

    virtual int get_num_nbc() const { return num_nbc; }

    virtual void print_info() const;

  private:
    const double thred_time;
    
    int num_nbc;

    std::vector<double> target_flow_rate;
};

#endif
