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
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "IFlowRate.hpp"

class FlowRate_Unsteady : public IFlowRate
{
  public:
    FlowRate_Unsteady( const std::string &filename );

    virtual ~FlowRate_Unsteady() = default;

    virtual double get_flow_rate( const int &nbc_id, const double &time ) const;

    virtual int get_num_nbc() const { return num_nbc; }

    virtual void print_info() const;

  private:
    int num_nbc;

    std::vector< std::vector<double> > coef_a, coef_b;

    std::vector<int> num_of_mode;

    std::vector<double> w, period;
};

#endif
