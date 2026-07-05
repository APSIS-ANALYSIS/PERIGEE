#ifndef HIMODELFACTORY_HPP
#define HIMODELFACTORY_HPP
// ============================================================================
// HIModelFactory.hpp
// ============================================================================
#include <memory>
#include "Sys_Tools.hpp"
#include "HIModel_PowerLaw.hpp"

class HIModelFactory
{
  public:
    static std::unique_ptr<IHIModel> createModel(const std::string &model_name)
    {
      if(model_name == "none" || model_name == "NONE" || model_name == "None")
        return nullptr;

      auto model = HIModel_PowerLaw::createPreset(model_name);
      if(model != nullptr) return model;

      SYS_T::print_fatal("Error: unknown HI model name.\n");
      return nullptr;
    }
};

#endif
