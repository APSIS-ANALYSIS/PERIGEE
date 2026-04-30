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
  //   field = 1, 2, 3 denotes x-, y-, z-direction displacement, respectively.
  //   disp is the imposed displacement value; velo and acce are its first and
  //   second time derivatives used to set dot_disp and dot_velo consistently.
  // --------------------------------------------------------------------------
  inline void disp_loading( const int &field, const double &tt,
      double &disp, double &velo, double &acce )
  {
    switch(field)
    {
      // x direction
      case 1:
        disp = 0.0;
        velo = 0.0;
        acce = 0.0;
        break;

      // y direction
      case 2:
        disp = 0.0;
        velo = 0.0;
        acce = 0.0;
        break;

      // z direction
      case 3:
        disp = tt;
        velo = 1.0;
        acce = 0.0;
        break;

      default:
        disp = 0.0;
        velo = 0.0;
        acce = 0.0;
        break;
    }
  }
}

#endif
