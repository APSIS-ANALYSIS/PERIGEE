#ifndef MATERIALMODELDATA_HPP
#define MATERIALMODELDATA_HPP
// --------------------------------------------------------------------------
// Centralized material parameters for examples/solids.
//
// Design goal:
//   Keep one single source of truth for material-model construction so that
//   drivers, local assembly routines, and post-processing use exactly the
//   same constitutive setup.
// --------------------------------------------------------------------------
#include "MaterialModel_Mixed_Elasticity.hpp"
#include "MaterialModel_ich_NeoHookean.hpp"
#include "MaterialModel_vol_Incompressible.hpp"

namespace MaterialModelData
{
  constexpr double density = 1.0e3;
  constexpr double shear_modulus = 6.666666666e4;

  inline std::unique_ptr<IMaterialModel_vol> create_vol_model()
  {
    return SYS_T::make_unique<MaterialModel_vol_Incompressible>( density );
  }

  inline std::unique_ptr<IMaterialModel_ich> create_isochoric_model()
  {
    return SYS_T::make_unique<MaterialModel_ich_NeoHookean>( shear_modulus );
  }

  inline std::unique_ptr<MaterialModel_Mixed_Elasticity> create_mixed_model()
  {
    return SYS_T::make_unique<MaterialModel_Mixed_Elasticity>(
        create_vol_model(), create_isochoric_model() );
  }
}

#endif
