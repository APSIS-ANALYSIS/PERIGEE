#ifndef SOLIDS_LOADDATA_HPP
#define SOLIDS_LOADDATA_HPP
// ============================================================================
// LoadData.hpp
//
// Problem-dependent loading data for solid dynamics.
// ============================================================================
#include "Vector_3.hpp"

namespace LoadData
{
  // --------------------------------------------------------------------------
  // body_force
  //   Body force density per unit mass, i.e. acceleration-like term b(x,t).
  //   It is multiplied by rho0 during assembly where needed.
  // --------------------------------------------------------------------------
  inline Vector_3 body_force( const Vector_3 &pt, const double &tt )
  {
    (void)pt;
    (void)tt;
    return Vector_3(0.0, 0.0, 0.0);
  }

  // --------------------------------------------------------------------------
  // surface_traction
  //   Prescribed traction vector t(x,t) on natural boundary faces.
  //   ebc_id is the local natural-boundary tag from ALocal_EBC.
  // --------------------------------------------------------------------------
  inline Vector_3 surface_traction( const int &ebc_id, const Vector_3 &pt,
      const double &tt, const Vector_3 &n_out )
  {
    (void)pt;
    (void)tt;

    switch(ebc_id)
    {
      case 0:
      {
        const double p0 = 0.0;
        return Vector_3( p0 * n_out.x(), p0 * n_out.y(), p0 * n_out.z() );
      }
      default:
        return Vector_3(0.0, 0.0, 0.0);
    }
  }

  // --------------------------------------------------------------------------
  // disp_loading
  //   Prescribed displacement-driven loading for Dirichlet nodes.
  //   disp, velo, and acce are returned for x-, y-, and z-directions together.
  //   Their first and second time derivatives are used to set dot_disp and
  //   dot_velo consistently.
  // --------------------------------------------------------------------------
  inline void disp_loading( const double &tt,
      Vector_3 &disp, Vector_3 &velo, Vector_3 &acce )
  {
    disp = Vector_3(0.0, 0.0, tt);
    velo = Vector_3(0.0, 0.0, 1.0);
    acce = Vector_3(0.0, 0.0, 0.0);
  }
}

#endif
