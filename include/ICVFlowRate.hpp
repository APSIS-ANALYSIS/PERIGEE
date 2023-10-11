#ifndef ICVFLOWRATE_HPP
#define ICVFLOWRATE_HPP
// ============================================================================
// ICVFlowRate.hpp
//
// Interface file for the cardiovascular inflow flow rate function.
//
// This function shall have two instantiations: steady case and unsteady case.
//
// Author: Ju Liu
// Date created: Aug. 6 2017
// ============================================================================
#include <sstream>
#include <string>

class ICVFlowRate
{
  public:
    ICVFlowRate() = default;

    virtual ~ICVFlowRate() = default;

    virtual double get_flow_rate( const int &nbc_id, const double &time ) const = 0;

    virtual int get_num_nbc() const = 0;

    virtual void print_info() const = 0;

  protected:
    // ------------------------------------------------------------------------
    // Generate a filename for inlet face nbc_id as Inlet_xxx_flowrate.txt
    // ------------------------------------------------------------------------
    virtual std::string gen_flowfile_name(const int &nbc_id) const
    {
      std::ostringstream ss;
      ss << "Inlet_";
      if( nbc_id/10 == 0 ) ss << "00";
      else if( nbc_id/100 == 0 ) ss << "0";

      ss << nbc_id << "_precalculated_flowrate.txt";

      return ss.str();
    }
};

#endif
