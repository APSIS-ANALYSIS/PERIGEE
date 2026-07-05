#ifndef IHIMODEL_HPP
#define IHIMODEL_HPP
// ============================================================================
// IHIModel.hpp
//
// Interface for hemolysis-index source models used in scalar transport
// equations.
//
// Date: Jul. 5 2026
// ============================================================================
#include <algorithm>
#include <cmath>
#include <string>
#include "Sys_Tools.hpp"

class IHIModel
{
  public:
    IHIModel() = default;

    virtual ~IHIModel() = default;

    virtual void print_info() const = 0;

    virtual std::string get_model_name() const
    {
      SYS_T::commPrint("Warning: IHIModel::get_model_name() is not implemented.\n");
      return "unknown";
    }

    virtual double get_C() const = 0;

    virtual double get_alpha() const = 0;

    virtual double get_beta() const = 0;

    virtual double get_HI(const double &tau_s, const double &expo_time) const
    {
      const double tau_eff = std::max(0.0, tau_s);
      const double time_eff = std::max(0.0, expo_time);
      return get_C() * std::pow(time_eff, get_alpha()) * std::pow(tau_eff, get_beta());
    }

    virtual double get_source(const double &tau_s) const
    {
      const double tau_eff = std::max(0.0, tau_s);
      return std::pow(get_C(), 1.0 / get_beta()) *
             std::pow(tau_eff, get_alpha() / get_beta());
    }
};

#endif
