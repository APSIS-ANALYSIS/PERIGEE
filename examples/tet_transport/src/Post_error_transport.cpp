#include "Post_error_transport.hpp"

double POST_ERROR_T::exact_sol( const double &x, const double &y, const double &z, const double &time )
{
  const double pi = MATH_T::PI;
  return time * time * sin(pi*x) * sin(pi*y) * sin(pi*z);
}


double POST_ERROR_T::exact_sol_x( const double &x, const double &y, const double &z, const double &time )
{
  return 0.0;
}


double POST_ERROR_T::exact_sol_y( const double &x, const double &y, const double &z, const double &time )
{
  return 0.0;
}


double POST_ERROR_T::exact_sol_z( const double &x, const double &y, const double &z, const double &time )
{
  return 0.0;
}


double POST_ERROR_T::get_manu_sol_error(
    const double &time,
    const double * const &solu,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const IQuadPts * const &quad )
{
  double error = 0.0;

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const std::vector<double> R = element -> get_R( qua );

    double sol = 0.0, coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    for(int ii=0; ii<element-> get_nLocBas(); ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];
      sol    += solu[ii]       * R[ii];
    }

    const double exact = exact_sol(coor_x, coor_y, coor_z, time );
    
    error += (sol - exact) * (sol - exact) * gwts;
  }

  return error;
}

// EOF
