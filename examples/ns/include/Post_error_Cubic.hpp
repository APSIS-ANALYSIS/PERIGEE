#ifndef POST_ERROR_CUBIC_HPP
#define POST_ERROR_CUBIC_HPP

#include "FEAElement.hpp"
#include "IQuadPts.hpp"
#include "Math_Tools.hpp"

namespace POST_ERROR_C
{
  inline double exact_sol_uz ( const double &x, const double &y, const double &z, const double &v_ave, const double &R );

  inline double exact_sol_duz_dx ( const double &x, const double &y, const double &z, const double &v_ave, const double &R );

  inline double exact_sol_duz_dy ( const double &x, const double &y, const double &z, const double &v_ave, const double &R );

  inline double exact_sol_duz_dz ( const double &x, const double &y, const double &z, const double &v_ave, const double &R )
  {return 0.0;}

  inline double exact_sol_dduz_dxdx ( const double &x, const double &y, const double &z, const double &v_ave, const double &R );

  inline double exact_sol_dduz_dydy ( const double &x, const double &y, const double &z, const double &v_ave, const double &R );

  inline double exact_sol_dduz_dxdy ( const double &x, const double &y, const double &z, const double &v_ave, const double &R );

  inline Vector_3 exact_sol_u ( const double &x, const double &y, const double &z, const double &v_ave, const double &R )
  {
    return Vector_3( 0.0, 0.0, exact_sol_uz(x,y,z,v_ave,R) );
  }

  inline Vector_3 exact_sol_u_dx ( const double &x, const double &y, const double &z, const double &v_ave, const double &R )
  {
    return Vector_3( 0.0, 0.0, exact_sol_duz_dx(x,y,z,v_ave,R) );
  }

  inline Vector_3 exact_sol_u_dy ( const double &x, const double &y, const double &z, const double &v_ave, const double &R )
  {
    return Vector_3( 0.0, 0.0, exact_sol_duz_dy(x,y,z,v_ave,R) );
  }

  inline Vector_3 exact_sol_u_dz( const double &x, const double &y, const double &z, const double &v_ave, const double &R )
  {
    return Vector_3( 0.0, 0.0, exact_sol_duz_dz(x,y,z,v_ave,R) );
  }

  double get_manu_sol_u_error(
      const double * const &solu_ux,
      const double * const &solu_uy,
      const double * const &solu_uz,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &V_average,
      const double &Radius );

  double get_manu_sol_u_errorH1(
      const double * const &solu_ux,
      const double * const &solu_uy,
      const double * const &solu_uz,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &V_average,
      const double &Radius );

  double get_exact_sol_u_normH2(
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &V_average,
      const double &Radius );
}

#endif