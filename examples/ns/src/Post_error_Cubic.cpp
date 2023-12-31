#include "Post_error_Cubic.hpp"

double POST_ERROR_C::exact_sol_uz( const double &x, const double &y, const double &z, const double &v_ave, const double &R )
{
  const double r = std::sqrt(x * x + y * y);
  return (5.0 * v_ave / 3.0) * (1 - (x*x + y*y) * r / (R*R*R));
}

double POST_ERROR_C::exact_sol_duz_dx( const double &x, const double &y, const double &z, const double &v_ave, const double &R )
{
  const double r = std::sqrt(x * x + y * y);
  return -5.0 * v_ave * r * x / (R*R*R);
}

double POST_ERROR_C::exact_sol_duz_dy( const double &x, const double &y, const double &z, const double &v_ave, const double &R )
{
  const double r = std::sqrt(x * x + y * y);
  return -5.0 * v_ave * r * y / (R*R*R);
}

double POST_ERROR_C::exact_sol_dduz_dxdx( const double &x, const double &y, const double &z, const double &v_ave, const double &R )
{
  const double r = std::sqrt(x * x + y * y);
  return -5.0 * v_ave * (r + x*x / r) / (R*R*R);
}

double POST_ERROR_C::exact_sol_dduz_dydy( const double &x, const double &y, const double &z, const double &v_ave, const double &R )
{
  const double r = std::sqrt(x * x + y * y);
  return -5.0 * v_ave * (r + y*y / r) / (R*R*R);
}

double POST_ERROR_C::exact_sol_dduz_dxdy( const double &x, const double &y, const double &z, const double &v_ave, const double &R )
{
  const double r = std::sqrt(x * x + y * y);
  return -5.0 * v_ave * x * y / (R*R*R * r);
}

double POST_ERROR_C::get_manu_sol_u_error(
      const double * const &solu_ux,
      const double * const &solu_uy,
      const double * const &solu_uz,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &V_average,
      const double &Radius )
{
  double error = 0.0;

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const std::vector<double> R = element -> get_R( qua );

    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    Vector_3 sol (0.0, 0.0, 0.0);

    for(int ii=0; ii<element-> get_nLocBas(); ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];

      sol.x()+= solu_ux[ii]    * R[ii];
      sol.y()+= solu_uy[ii]    * R[ii];
      sol.z()+= solu_uz[ii]    * R[ii];
    }

    const Vector_3 exact = exact_sol_u( coor_x, coor_y, coor_z, V_average, Radius );
    
    error += (sol - exact).dot_product(sol - exact) * gwts;
  }

  return error;
}

double POST_ERROR_C::get_manu_sol_u_errorH1(
      const double * const &solu_ux,
      const double * const &solu_uy,
      const double * const &solu_uz,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &V_average,
      const double &Radius )
{
  double errorH1 = 0.0;

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const std::vector<double> R = element -> get_R( qua );

    const std::vector<double> R_dx = element -> get_dR_dx( qua );
    
    const std::vector<double> R_dy = element -> get_dR_dy( qua );
    
    const std::vector<double> R_dz = element -> get_dR_dz( qua );
    
    Vector_3 sol (0.0, 0.0, 0.0);
    Vector_3 sol_dx (0.0, 0.0, 0.0);
    Vector_3 sol_dy (0.0, 0.0, 0.0);
    Vector_3 sol_dz (0.0, 0.0, 0.0);

    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    for(int ii=0; ii<element-> get_nLocBas(); ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];

      sol.x()+= solu_ux[ii]    * R[ii];
      sol.y()+= solu_uy[ii]    * R[ii];
      sol.z()+= solu_uz[ii]    * R[ii];

      sol_dx.x() += solu_ux[ii] * R_dx[ii];
      sol_dx.y() += solu_uy[ii] * R_dx[ii];
      sol_dx.z() += solu_uz[ii] * R_dx[ii];

      sol_dy.x() += solu_ux[ii] * R_dy[ii];
      sol_dy.y() += solu_uy[ii] * R_dy[ii];
      sol_dy.z() += solu_uz[ii] * R_dy[ii];

      sol_dz.x() += solu_ux[ii] * R_dz[ii];
      sol_dz.y() += solu_uy[ii] * R_dz[ii];
      sol_dz.z() += solu_uz[ii] * R_dz[ii];
    }

    const Vector_3 exact = exact_sol_u( coor_x, coor_y, coor_z, V_average, Radius );
    const Vector_3 exact_dx = exact_sol_u_dx( coor_x, coor_y, coor_z, V_average, Radius );
    const Vector_3 exact_dy = exact_sol_u_dy( coor_x, coor_y, coor_z, V_average, Radius );
    const Vector_3 exact_dz = exact_sol_u_dz( coor_x, coor_y, coor_z, V_average, Radius );

    errorH1 += (sol_dx - exact_dx).dot_product(sol_dx - exact_dx) * gwts 
	     + (sol_dy - exact_dy).dot_product(sol_dy - exact_dy) * gwts 
	     + (sol_dz - exact_dz).dot_product(sol_dz - exact_dz) * gwts 
	     + (sol - exact).dot_product(sol - exact) * gwts;
  }

  return errorH1;
}

double POST_ERROR_C::get_exact_sol_u_normH2(
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const IQuadPts * const &quad,
    const double &V_average,
    const double &Radius )
{
  double normH2 = 0.0;

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const std::vector<double> R = element -> get_R( qua );
    
    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    for(int ii=0; ii<element-> get_nLocBas(); ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];
    }

    const Vector_3 exact = exact_sol_u( coor_x, coor_y, coor_z, V_average, Radius );
    const Vector_3 exact_dx = exact_sol_u_dx( coor_x, coor_y, coor_z, V_average, Radius );
    const Vector_3 exact_dy = exact_sol_u_dy( coor_x, coor_y, coor_z, V_average, Radius );
    const Vector_3 exact_dz = exact_sol_u_dz( coor_x, coor_y, coor_z, V_average, Radius );

    const double uz_dxdx = exact_sol_dduz_dxdx( coor_x, coor_y, coor_z, V_average, Radius );
    const double uz_dydy = exact_sol_dduz_dydy( coor_x, coor_y, coor_z, V_average, Radius );
    const double uz_dxdy = exact_sol_dduz_dxdy( coor_x, coor_y, coor_z, V_average, Radius );

    normH2 += exact_dx.dot_product(exact_dx) * gwts 
	     + exact_dy.dot_product(exact_dy) * gwts 
	     + exact_dz.dot_product(exact_dz) * gwts 
	     + exact.dot_product(exact) * gwts
       + (uz_dxdx * uz_dxdx + uz_dydy * uz_dydy + 2 * uz_dxdy * uz_dxdy) * gwts; // non-trivial second order derivative
  }

  return normH2;
}