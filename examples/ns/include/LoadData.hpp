#ifndef NS_LOADDATA_HPP
#define NS_LOADDATA_HPP

#include "Vector_3.hpp"
#include "ALocal_InflowBC.hpp"
#include "IFlowRate.hpp"
#include "PDNSolution.hpp"
#include "Math_Tools.hpp"

namespace LoadData
{
  // --------------------------------------------------------------------------
  // body_force
  //   Body force density per unit mass, i.e. acceleration-like term b(x,t).
  //   It is multiplied by rho0 during assembly where needed.
  // --------------------------------------------------------------------------
  inline Vector_3 body_force( const Vector_3 &pt, const double &tt )
  {
    return Vector_3(0.0, 0.0, 0.0);
  }

  // --------------------------------------------------------------------------
  // ebc_traction
  //   Surface traction h(x,t,n) applied on an EBC face identified by ebc_id.
  //   The return value is the traction vector in physical coordinates.
  // --------------------------------------------------------------------------
  inline Vector_3 ebc_traction( const int &ebc_id, const Vector_3 &pt,
      const double &tt, const Vector_3 &n_out )
  {
    switch(ebc_id)
    {
      case 0:
      {
        const double p0 = 0.0;
        return Vector_3( p0*n_out.x(), p0*n_out.y(), p0*n_out.z() );
      }
      default:
        return Vector_3(0.0, 0.0, 0.0);
    }
  }

  // --------------------------------------------------------------------------
  // rescale_inflow_value
  //   Rescale the baseline inflow velocity profile using time-dependent
  //   flow-rate factors and optional turbulence-intensity perturbations.
  // --------------------------------------------------------------------------
  inline void rescale_inflow_value( const double &stime,
      const ALocal_InflowBC * const &infbc,
      const IFlowRate * const &flrate,
      const PDNSolution * const &sol_base,
      PDNSolution * const &sol )
  {
    const int num_nbc = infbc -> get_num_nbc();

    for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
    {
      const int numnode = infbc -> get_Num_LD( nbc_id );

      const double factor  = flrate -> get_flow_rate( nbc_id, stime );
      const double std_dev = flrate -> get_flow_TI_std_dev( nbc_id );

      for(int ii=0; ii<numnode; ++ii)
      {
        const int node_index = infbc -> get_LDN( nbc_id, ii );

        const int base_idx[3] = { node_index*4+1, node_index*4+2, node_index*4+3 };

        double base_vals[3];

        VecGetValues(sol_base->solution, 3, base_idx, base_vals);

        const double perturb_x = MATH_T::gen_double_rand_normal(0, std_dev);
        const double perturb_y = MATH_T::gen_double_rand_normal(0, std_dev);
        const double perturb_z = MATH_T::gen_double_rand_normal(0, std_dev);

        const double vals[3] = { base_vals[0] * factor * (1.0 + perturb_x),
          base_vals[1] * factor * (1.0 + perturb_y),
          base_vals[2] * factor * (1.0 + perturb_z) };

        VecSetValues(sol->solution, 3, base_idx, vals, INSERT_VALUES);
      }
    }

    sol->Assembly_GhostUpdate();
  }

}

#endif
