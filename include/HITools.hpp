#ifndef HITOOLS_HPP
#define HITOOLS_HPP
// ============================================================================
// HITools.hpp
//
// Utilities for converting between the linearized transport variable D and
// the hemolysis index HI = D^beta.
//
// Date: Jul. 5 2026
// ============================================================================
#include <algorithm>
#include <cmath>
#include <vector>
#include "PDNSolution.hpp"

namespace HI_T
{
  inline double D_to_HI(const double &D, const double &beta)
  {
    return std::pow(std::max(0.0, D), beta);
  }

  inline double HI_to_D(const double &HI, const double &beta)
  {
    return std::pow(std::max(0.0, HI), 1.0 / beta);
  }

  inline void transform_scalar_solution(
      PDNSolution * const sol, const double &beta, const bool is_D_to_HI )
  {
    SYS_T::print_fatal_if(sol == nullptr, "Error: HI_T::transform_scalar_solution received nullptr.\n");
    SYS_T::print_fatal_if(sol->get_dof_num() != 1, "Error: HI_T::transform_scalar_solution expects a scalar solution.\n");

    const int nlocal = sol->get_nlocal();
    const std::vector<double> local_array = sol->GetLocalArray();

    std::vector<PetscInt> index(nlocal, 0);
    std::vector<PetscScalar> value(nlocal, 0.0);

    for(int ii=0; ii<nlocal; ++ii)
    {
      index[ii] = ii;
      value[ii] = is_D_to_HI ? D_to_HI(local_array[ii], beta) : HI_to_D(local_array[ii], beta);
    }

    VecSetValuesLocal(sol->solution, nlocal, index.data(), value.data(), INSERT_VALUES);
    sol->Assembly_GhostUpdate();
  }

  inline void convert_D_to_HI(PDNSolution * const sol, const double &beta)
  {
    transform_scalar_solution(sol, beta, true);
  }

  inline void convert_HI_to_D(PDNSolution * const sol, const double &beta)
  {
    transform_scalar_solution(sol, beta, false);
  }
}

#endif
