#ifndef GENBCFACTORY_HPP
#define GENBCFACTORY_HPP
// ============================================================================
// GenBCFactory.hpp
//
// This is the factory that constructs a GenBC object based on the file type and
// time integration parameters.
//
// Author: Ju Liu, liujuy@gmail.com
// Date: Feb. 4 2025
// ============================================================================
#include "GenBC_Resistance.hpp"
#include "GenBC_RCR.hpp"
#include "GenBC_Inductance.hpp"
#include "GenBC_Coronary.hpp"
#include "GenBC_Pressure.hpp"

class GenBCFactory
{
  public:
    static std::unique_ptr<IGenBC> createGenBC(const std::string &lpn_file,
        const double &initial_time = 0.0, const double &initial_step = 0.1,
        const int &initial_index = 0, const int &num_steps_0d = 1000)
    {
      // Retrieve the file type
      const int file_type = GENBC_T::get_genbc_file_type(lpn_file);

      switch (file_type)
      {
        case 1:
          return SYS_T::make_unique<GenBC_Resistance>(lpn_file);
        case 2:
          return SYS_T::make_unique<GenBC_RCR>(lpn_file, num_steps_0d, initial_step);
        case 3:
          return SYS_T::make_unique<GenBC_Inductance>(lpn_file);
        case 4:
          return SYS_T::make_unique<GenBC_Coronary>(lpn_file, num_steps_0d, initial_step, initial_index);
        case 5:
          return SYS_T::make_unique<GenBC_Pressure>(lpn_file, initial_time);
        default:
          SYS_T::print_fatal("Error: GenBC input file %s format cannot be recognized.\n", lpn_file.c_str());
          return nullptr;
      }
    }
};

#endif
