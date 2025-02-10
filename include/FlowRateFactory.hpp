#ifndef FLOWRATEFACTORY_HPP
#define FLOWRATEFACTORY_HPP

#include "FlowRate_Steady.hpp"
#include "FlowRate_Unsteady.hpp"
#include "FlowRate_Linear2Steady.hpp"
#include "FlowRate_Cosine2Steady.hpp"

class FlowRateFactory
{
  public:
    static std::unique_ptr<IFlowRate> createFlowRate( const std::string &inflow_filename )
    {
      // Retrieve the file type
      const int file_type = get_flowrate_file_type(inflow_filename);

      switch (file_type)
      {
        case 1:
          return SYS_T::make_unique<FlowRate_Steady>(inflow_filename);
        case 2:
          return SYS_T::make_unique<FlowRate_Unsteady>(inflow_filename);
        case 3:
          return SYS_T::make_unique<FlowRate_Linear2Steady>(inflow_filename);
        case 4:
          return SYS_T::make_unique<FlowRate_Cosine2Steady>(inflow_filename);

        default:
          SYS_T::print_fatal("Error: Inflow input file %s format cannot be recognized.\n", inflow_filename.c_str());
          return nullptr;
      } 
    }

  private:
    static int get_flowrate_file_type( const std::string &inflow_filename )
    {
      // open the file
      SYS_T::file_check( inflow_filename );

      std::ifstream reader;
      reader.open( inflow_filename.c_str(), std::ifstream::in );

      // Read the fist non-comment line of the file
      std::istringstream sstrm;
      std::string sline, bc_type;

      while( std::getline(reader, sline) )
      {
        if( sline[0] != '#' && !sline.empty() )
        {
          sstrm.str(sline);
          sstrm >> bc_type;
          sstrm.clear();
          break;
        }
      }

      reader.close();

      // Now compare and return the corresponding type value
      if( bc_type.compare("Steady") == 0
          || bc_type.compare("steady") == 0
          || bc_type.compare("STEADY") == 0 )
      {
        return 1;
      }
      else if( bc_type.compare("Unsteady") ==0
          || bc_type.compare("unsteady") == 0
          || bc_type.compare("UNSTEADY") == 0 )
      {
        return 2;
      }
      else if( bc_type.compare("Linear") ==0
          || bc_type.compare("linear") == 0
          || bc_type.compare("LINEAR") == 0 )
      {
        return 3;
      }
      else if( bc_type.compare("Cosine") ==0
          || bc_type.compare("cosine") == 0
          || bc_type.compare("COSINE") == 0 )
      {
        return 4;
      }       
      else
      {
        return 0;
      }
    }
};

#endif
