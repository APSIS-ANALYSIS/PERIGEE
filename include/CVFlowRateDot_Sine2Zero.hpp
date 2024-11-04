#ifndef CVFLOWRATEDOT_SINE2ZERO_HPP
#define CVFLOWRATEDOT_SINE2ZERO_HPP
// ============================================================================
// CVFlowRateDot_Sine2Zero.hpp
//
// It is designed to match CVFlowRate_Cosine2Steady, using a sine function
// to prescibe the starting dot flow rate
//
// Date Created: Nov. 04 2024
// ============================================================================
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "ICVFlowRate.hpp"

class CVFlowRateDot_Sine2Zero : public ICVFlowRate
{
  public:
    // This constructor will set a uniform thred_time and dot flrate for all inlets.
    // Same as CVFlowRate_Linear2Steady.
    CVFlowRateDot_Sine2Zero( const int &input_num_nbc,
        const double &in_thred_time, const double &flrate,
        const double &in_TI_std_dev );

    // This constructor will set a uniform thred_time for all inlets, and the
    // target flow rate are determined from the file and set to be the flow rate
    // at the initial time (i.e. time = 0.0) by setting target_flow_rate to be
    // the sum of coef_a.
    // Same as CVFlowRate_Linear2Steady.
    CVFlowRateDot_Sine2Zero( const double &in_thred_time, 
        const double &in_TI_std_dev, const std::string &filename );

    virtual ~CVFlowRateDot_Sine2Zero() = default;

    // From time 0 to in_time, the dot flow rate = 0.5 * flrate * PI / in_thred_time * sin (PI * time / in_thred_time))
    // From in_time to infty, dot_flow_rate = 0.0
    virtual double get_flow_rate( const int &nbc_id, const double &time ) const;

    virtual double get_flow_TI_std_dev( const int &nbc_id ) const { return TI_std_dev; }

    virtual int get_num_nbc() const { return num_nbc; }

    virtual void print_info() const;

  private:
    int num_nbc;
    
    const double thred_time, TI_std_dev;
    
    std::vector<double> target_flow_rate;
};

#endif
