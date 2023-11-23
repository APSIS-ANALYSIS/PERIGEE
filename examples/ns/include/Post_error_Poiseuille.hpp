#ifndef POST_ERROR_POISEULLE_HPP
#define POST_ERROR_POISEULLE_HPP

#include "FEAElement.hpp"
#include "IQuadPts.hpp"
#include "Math_Tools.hpp"

namespace POST_ERROR_P
{
  // coef = P_diff / (4 * fl_mu * L)

  inline double exact_sol_ux ( const double &x, const double &y, const double &z, const double &coef, const double &R )
  {return 0.0;}

  inline double exact_sol_dux_dx ( const double &x, const double &y, const double &z, const double &coef )
  {return 0.0;}

  inline double exact_sol_dux_dy ( const double &x, const double &y, const double &z, const double &coef )
  {return 0.0;}

  inline double exact_sol_dux_dz ( const double &x, const double &y, const double &z, const double &coef )
  {return 0.0;}

  inline double exact_sol_uy ( const double &x, const double &y, const double &z, const double &coef, const double &R )
  {return 0.0;}

  inline double exact_sol_duy_dx ( const double &x, const double &y, const double &z, const double &coef )
  {return 0.0;}

  inline double exact_sol_duy_dy ( const double &x, const double &y, const double &z, const double &coef )
  {return 0.0;}

  inline double exact_sol_duy_dz ( const double &x, const double &y, const double &z, const double &coef )
  {return 0.0;}

  inline double exact_sol_uz ( const double &x, const double &y, const double &z, const double &coef, const double &R );

  inline double exact_sol_duz_dx ( const double &x, const double &y, const double &z, const double &coef );

  inline double exact_sol_duz_dy ( const double &x, const double &y, const double &z, const double &coef );

  inline double exact_sol_duz_dz ( const double &x, const double &y, const double &z, const double &coef )
  {return 0.0;}

  inline Vector_3 exact_sol_u ( const double &x, const double &y, const double &z, const double &coef, const double &R )
  {
    return Vector_3( exact_sol_ux(x,y,z,coef,R), exact_sol_uy(x,y,z,coef,R), exact_sol_uz(x,y,z,coef,R) );
  }

  inline Vector_3 exact_sol_u_dx ( const double &x, const double &y, const double &z, const double &coef )
  {
    return Vector_3( exact_sol_dux_dx(x,y,z,coef), exact_sol_duy_dx(x,y,z,coef), exact_sol_duz_dx(x,y,z,coef) );
  }

  inline Vector_3 exact_sol_u_dy ( const double &x, const double &y, const double &z, const double &coef )
  {
    return Vector_3( exact_sol_dux_dy(x,y,z,coef), exact_sol_duy_dy(x,y,z,coef), exact_sol_duz_dy(x,y,z,coef) );
  }

  inline Vector_3 exact_sol_u_dz( const double &x, const double &y, const double &z, const double &coef )
  {
    return Vector_3( exact_sol_dux_dz(x,y,z,coef), exact_sol_duy_dz(x,y,z,coef), exact_sol_duz_dz(x,y,z,coef) );
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
      const double &fl_mu,
      const double &Radius,
      const double &Length,
      const double &P_diff );

  double get_manu_sol_u_errorH1(
      const double * const &solu_ux,
      const double * const &solu_uy,
      const double * const &solu_uz,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &fl_mu,
      const double &Radius,
      const double &Length,
      const double &P_diff );

  double get_exact_sol_u_normH2(
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &fl_mu,
      const double &Radius,
      const double &Length,
      const double &P_diff );
}

#endif