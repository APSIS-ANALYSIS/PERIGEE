#include "Post_error_transport.hpp"

double POST_ERROR_T::exact_sol( const double &x, const double &y, const double &z, const double &time )
{
  //const double pi = MATH_T::PI;
  return time * time * time * time * x*(x-1)*y*(y-1)*z*(z-1);
}

double POST_ERROR_T::exact_sol_dx( const double &x, const double &y, const double &z, const double &time )
{
  //const double pi = MATH_T::PI;
  return time * time * time * time * (2.0*x-1)*y*(y-1)*z*(z-1);
}

double POST_ERROR_T::exact_sol_dy( const double &x, const double &y, const double &z, const double &time )
{
  //const double pi = MATH_T::PI;
  return time * time * time * time * x*(x-1)*(2.0*y-1)*z*(z-1);
}

double POST_ERROR_T::exact_sol_dz( const double &x, const double &y, const double &z, const double &time )
{
  //const double pi = MATH_T::PI;
  return time * time * time * time * x*(x-1)*y*(y-1)*(2.0*z-1);
}

double POST_ERROR_T::get_manu_sol_errorL2(
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

double POST_ERROR_T::get_manu_sol_errorH1(
    const double &time,
    const double * const &solu,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const IQuadPts * const &quad )
{
  double errorH1 = 0.0;

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const std::vector<double> R = element -> get_R( qua );

    const std::vector<double> R_dx = element -> get_dR_dx( qua );
    
    const std::vector<double> R_dy = element -> get_dR_dy( qua );
    
    const std::vector<double> R_dz = element -> get_dR_dz( qua );
    
    double sol = 0.0, sol_dx = 0.0, sol_dy = 0.0, sol_dz = 0.0;
    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    for(int ii=0; ii<element-> get_nLocBas(); ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];
      sol    += solu[ii]       * R[ii];
      sol_dx += solu[ii]    * R_dx[ii];
      sol_dy += solu[ii]    * R_dy[ii];
      sol_dz += solu[ii]    * R_dz[ii];
    }

    const double exact = exact_sol(coor_x, coor_y, coor_z, time );
    const double exact_dx = exact_sol_dx(coor_x, coor_y, coor_z, time );
    const double exact_dy = exact_sol_dy(coor_x, coor_y, coor_z, time );
    const double exact_dz = exact_sol_dz(coor_x, coor_y, coor_z, time );

    errorH1 += (sol_dx - exact_dx) * (sol_dx - exact_dx) * gwts 
	     + (sol_dy - exact_dy) * (sol_dy - exact_dy) * gwts 
	     + (sol_dz - exact_dz) * (sol_dz - exact_dz) * gwts 
	     + (sol - exact) * (sol - exact) * gwts;
  }

  return errorH1;
}

// EOF