#ifndef FLOWRATE_UNSTEADY_HPP
#define FLOWRATE_UNSTEADY_HPP
// ==================================================================
// FlowRate_Unsteady.hpp
//
// Unsteady flow for inflow condition.
// 
// The period value should be compatible with the original data file.
//
// I use time / period to determine the past full cycles. Then
// time - num_of_past_period * period gives the time in the local
// cycle, which can be used to evaluate the fourier wave.
//
// Output flow rate is defined as follows
//   a_0 + sum_i a_i cos(i w t) + b_i sin(i w t) 
//       for i = 1 : num_of_mode
//           t = t - [t/period] x period
//   
// Author: Ju Liu
// Date Created: Sept. 23 2017
// ==================================================================
#include <vector>
#include "IFlowRate.hpp"

class FlowRate_Unsteady final : public IFlowRate
{
  public:
    FlowRate_Unsteady( const std::string &filename );

    ~FlowRate_Unsteady() override = default;

    double get_flow_rate( int nbc_id, double time ) const override;
    double get_dot_flow_rate( int nbc_id, double time ) const override;

    int get_num_nbc() const override { return num_nbc; }

    void print_info() const override;

  private:
    int num_nbc;

    std::vector< std::vector<double> > coef_a, coef_b;

    std::vector<int> num_of_mode;

    std::vector<double> w, period;
};

#endif
