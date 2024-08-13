#ifndef MATERIALMODEL_NEOHOOKEAN_HPP
#define MATERIALMODEL_NEOHOOKEAN_HPP

#include "IMaterialModel_new.hpp"

class MaterialModel_NeoHookean : public IMaterialModel_new
{
  public:
    MaterialModel_NeoHookean( std::unique_ptr<IMaterialModel_vol> vmodel) 
      : IMaterialModel_new(std::move(vmodel)) {}

  private:
    double elastic_E, elastic_nu, elastic_mu;

};

#endif
