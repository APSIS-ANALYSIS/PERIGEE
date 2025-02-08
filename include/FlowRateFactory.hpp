#ifndef FLOWRATEFACTORY_HPP
#define FLOWRATEFACTORY_HPP

#include "FlowRate_Steady.hpp"
#include "FlowRate_Unsteady.hpp"
#include "FlowRate_Linear2Steady.hpp"
#include "FlowRate_Cosine2Steady.hpp"
#include "FlowRate_Tools.hpp"

class FlowRateFactory
{
  public:
    static std::unique_ptr<IFlowRate> createFlowRate( const std::string &inflow_filename,
      const double &inflow_thd_time )
    {
      // Retrieve the file type
      const int file_type = FLOW_T::get_flowrate_file_type(inflow_filename);

      switch (file_type)
      {
        case 1:
          return SYS_T::make_unique<FlowRate_Steady>(inflow_filename);
        case 2:
          return SYS_T::make_unique<FlowRate_Unsteady>(inflow_filename);
        case 3:
          return SYS_T::make_unique<FlowRate_Linear2Steady>(inflow_thd_time, inflow_filename);
        case 4:
          return SYS_T::make_unique<FlowRate_Cosine2Steady>(inflow_thd_time, inflow_filename);

        default:
          SYS_T::print_fatal("Error: Inflow input file %s format cannot be recognized.\n", inflow_filename.c_str());
          return nullptr;
      } 
    }
};

#endif