#ifndef ICVFLOWRATE_HPP
#define ICVFLOWRATE_HPP
// ==================================================================
// ICVFlowRate.hpp
//
// Interface file for the cardiovascular inflow flow rate function.
//
// This function shall have two instantiations: steady case and 
// unsteady case.
//
// Author: Ju Liu
// Date created: Aug. 6 2017
// ==================================================================

class ICVFlowRate
{
  public:
    ICVFlowRate(){};

    virtual ~ICVFlowRate(){};

    virtual double get_flow_rate( const double &time ) const = 0;

    virtual void print_info() const = 0;
};

#endif
