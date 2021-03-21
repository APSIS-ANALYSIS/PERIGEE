#ifndef CVFLOWRATE_UNSTEADY_HPP
#define CVFLOWRATE_UNSTEADY_HPP
// ==================================================================
// CVFlowRate_Unsteady.hpp
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
#include "Sys_Tools.hpp"
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "ICVFlowRate.hpp"

class CVFlowRate_Unsteady : public ICVFlowRate
{
  public:
    CVFlowRate_Unsteady( const char * const &filename );

    virtual ~CVFlowRate_Unsteady();

    virtual double get_flow_rate(const double &time) const;

    virtual void print_info() const;

  private:
    std::vector<double> coef_a, coef_b;

    int num_of_mode;

    double w, period;
};

#endif
