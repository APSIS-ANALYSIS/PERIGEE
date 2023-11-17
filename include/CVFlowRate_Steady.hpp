#ifndef CVFLOWRATE_STEADY_HPP
#define CVFLOWRATE_STEADY_HPP
// ============================================================================
// CVFlowRate_Steady.hpp
//
// Steady flow for inflow condition
//
// Author: Ju Liu
// Date Created: May 1 2021
// ============================================================================
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "ICVFlowRate.hpp"

class CVFlowRate_Steady : public ICVFlowRate
{
  public:
    CVFlowRate_Steady( const std::string &filename );

    CVFlowRate_Steady( const int &input_num_nbc, const double &in_flowrate );

    virtual ~CVFlowRate_Steady() = default;

    virtual double get_flow_rate(const int &nbc_id, const double &time) const;

    virtual int get_num_nbc() const { return num_nbc; }

    virtual void print_info() const;

  private:
    int num_nbc;
    
    std::vector<double> flowrate;
};

#endif
