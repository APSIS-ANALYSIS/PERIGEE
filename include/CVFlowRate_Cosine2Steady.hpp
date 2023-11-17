#ifndef CVFLOWRATE_COSINE2STEADY_HPP
#define CVFLOWRATE_COSINE2STEADY_HPP
// ============================================================================
// CVFlowRate_Cosine2Steady.hpp
//
// It is similar to CVFlowRate_Linear2Steady, except using a cosine function
// to prescibe the starting flow rate
//
// Date Created: Nov. 17 2023
// ============================================================================
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "ICVFlowRate.hpp"

class CVFlowRate_Cosine2Steady : public ICVFlowRate
{
  public:
    // This constructor will set a uniform thred_time and flrate for all inlets.
    // Same as CVFlowRate_Linear2Steady.
    CVFlowRate_Cosine2Steady( const int &input_num_nbc,
        const double &in_thred_time, const double &flrate );

    // This constructor will set a uniform thred_time for all inlets, and the
    // target flow rate are determined from the file and set to be the flow rate
    // at the initial time (i.e. time = 0.0) by setting target_flow_rate to be
    // the sum of coef_a.
    // Same as CVFlowRate_Linear2Steady.
    CVFlowRate_Cosine2Steady( const double &in_thred_time, 
        const std::string &filename );

    virtual ~CVFlowRate_Cosine2Steady() = default;

    // From time 0 to in_time, the flow rate = 0.5 * flrate * (1 -  cos (PI * time / in_thred_time))
    // From in_time to infty, flow_rate = flrate
    virtual double get_flow_rate( const int &nbc_id, const double &time ) const;

    virtual int get_num_nbc() const { return num_nbc; }

    virtual void print_info() const;

  private:
    const double thred_time;
    
    int num_nbc;

    std::vector<double> target_flow_rate;
};

#endif