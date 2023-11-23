#include "Post_error_Poiseuille.hpp"

double POST_ERROR_P::exact_sol_uz( const double &x, const double &y, const double &z, const double &coef, const double &R )
{
    return coef * R * R - coef * x * x - coef * y * y;
}

double POST_ERROR_P::exact_sol_duz_dx( const double &x, const double &y, const double &z, const double &coef )
{
    return -2 * coef * x;
}

double POST_ERROR_P::exact_sol_duz_dy( const double &x, const double &y, const double &z, const double &coef)
{
    return -2 * coef * y;
}

double POST_ERROR_P::get_manu_sol_u_error(
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
      const double &P_diff )
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

    const double coefficent = 0.25 * P_diff / (fl_mu * Length);

    const Vector_3 exact = exact_sol_u( coor_x, coor_y, coor_z, coefficent, Radius );
    
    error += (sol - exact).dot_product(sol - exact) * gwts;
  }

  return error;
}

double POST_ERROR_P::get_manu_sol_u_errorH1(
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
      const double &P_diff )
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

      sol_dy.x() += solu_ux[ii] * R_dy[ii];
      sol_dy.y() += solu_uy[ii] * R_dy[ii];
      sol_dy.z() += solu_uz[ii] * R_dy[ii];
    }

    const double coefficent = 0.25 * P_diff / (fl_mu * Length);

    const Vector_3 exact = exact_sol_u( coor_x, coor_y, coor_z, coefficent, Radius );
    const Vector_3 exact_dx = exact_sol_u_dx( coor_x, coor_y, coor_z, coefficent );
    const Vector_3 exact_dy = exact_sol_u_dy( coor_x, coor_y, coor_z, coefficent );
    const Vector_3 exact_dz = exact_sol_u_dz( coor_x, coor_y, coor_z, coefficent );

    errorH1 += (sol_dx - exact_dx).dot_product(sol_dx - exact_dx) * gwts 
	     + (sol_dy - exact_dy).dot_product(sol_dy - exact_dy) * gwts 
	     + (sol_dz - exact_dz).dot_product(sol_dz - exact_dz) * gwts 
	     + (sol - exact).dot_product(sol - exact) * gwts;
  }

  return errorH1;
}

double POST_ERROR_P::get_exact_sol_u_normH2(
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const IQuadPts * const &quad,
    const double &fl_mu,
    const double &Radius,
    const double &Length,
    const double &P_diff )
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

    const double coefficent = 0.25 * P_diff / (fl_mu * Length);

    const Vector_3 exact = exact_sol_u( coor_x, coor_y, coor_z, coefficent, Radius );
    const Vector_3 exact_dx = exact_sol_u_dx( coor_x, coor_y, coor_z, coefficent );
    const Vector_3 exact_dy = exact_sol_u_dy( coor_x, coor_y, coor_z, coefficent );
    const Vector_3 exact_dz = exact_sol_u_dz( coor_x, coor_y, coor_z, coefficent );

    normH2 += exact_dx.dot_product(exact_dx) * gwts 
	     + exact_dy.dot_product(exact_dy) * gwts 
	     + exact_dz.dot_product(exact_dz) * gwts 
	     + exact.dot_product(exact) * gwts
       + 8 * coefficent * coefficent * gwts; // non-trivial second order derivative
  }

  return normH2;
}