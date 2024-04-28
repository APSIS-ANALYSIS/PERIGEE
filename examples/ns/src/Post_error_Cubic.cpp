#include "Post_error_Cubic.hpp"

Vector_3 POST_ERROR_C::exact_velocity ( const Vector_3 &pt, const double &Qt, const double &RR )
{
    const double Area = MATH_T::PI * RR * RR; 
    const double RR_3 = RR * RR * RR;

    const double yy = pt.y(), zz = pt.z();
    const double rr = std::sqrt(yy * yy + zz * zz);
    const double rr_3 = rr * rr * rr;

    const double fr  = 5 / (3 * Area) * (1 - rr_3 / RR_3 );
    
    Vector_3 velocity (Qt * fr, 0.0, 0.0);
    return velocity;
}

double  POST_ERROR_C::exact_pressure ( const Vector_3 &pt, const double &dp_dx, const double &Length )
{
    double pressure = 0.0 + dp_dx * (pt.x() - Length);
    return pressure;
}

Vector_3 POST_ERROR_C::exact_grad_u ( const Vector_3 &pt, const double &Qt, const double &RR )
{
    const double Area = MATH_T::PI * RR * RR; 
    const double RR_3 = RR * RR * RR;
    
    const double yy = pt.y(), zz = pt.z();
    const double rr = std::sqrt(yy * yy + zz * zz);

    const double df_dx = 0.0;
    const double df_dy = 5 / (3 * Area) * (-3.0 * rr * yy / RR_3 );
    const double df_dz = 5 / (3 * Area) * (-3.0 * rr * zz / RR_3 );

    Vector_3 grad_u (Qt * df_dx, Qt * df_dy, Qt * df_dz);
    return grad_u;
}

Vector_3 POST_ERROR_C::exact_grad_v ( const Vector_3 &pt, const double &Qt, const double &RR )
{
    Vector_3 grad_v (0.0, 0.0, 0.0);
    return grad_v;
}

Vector_3 POST_ERROR_C::exact_grad_w ( const Vector_3 &pt, const double &Qt, const double &RR )
{
    Vector_3 grad_w (0.0, 0.0, 0.0);
    return grad_w;
}

double POST_ERROR_C::get_manu_sol_u_errorL2(
      const double * const &sol_u,
      const double * const &sol_v,
      const double * const &sol_w,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &Qt,
      const double &RR )
{
  double error = 0.0;

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const std::vector<double> NR = element -> get_R( qua );

    Vector_3 pt (0.0, 0.0, 0.0);
    Vector_3 sol_velo (0.0, 0.0, 0.0);
    
    for(int ii=0; ii<element-> get_nLocBas(); ++ii)
    {
      pt.x() += ectrlPts_x[ii] * NR[ii];
      pt.y() += ectrlPts_y[ii] * NR[ii];
      pt.z() += ectrlPts_z[ii] * NR[ii];

      sol_velo.x() += sol_u[ii] * NR[ii];
      sol_velo.y() += sol_v[ii] * NR[ii];
      sol_velo.z() += sol_w[ii] * NR[ii];
    }

    const Vector_3 exact_velo = exact_velocity(pt, Qt, RR);
    
    // L2-norm of error
    error += (sol_velo - exact_velo).dot_product(sol_velo - exact_velo) * gwts;
  }

  return error;
}

double POST_ERROR_C::get_manu_sol_u_L2(
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &Qt,
      const double &RR )
{
  double norm = 0.0;

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const std::vector<double> NR = element -> get_R( qua );

    Vector_3 pt (0.0, 0.0, 0.0);
    Vector_3 sol_velo (0.0, 0.0, 0.0);
    
    for(int ii=0; ii<element-> get_nLocBas(); ++ii)
    {
      pt.x() += ectrlPts_x[ii] * NR[ii];
      pt.y() += ectrlPts_y[ii] * NR[ii];
      pt.z() += ectrlPts_z[ii] * NR[ii];
    }

    const Vector_3 exact_velo = exact_velocity(pt, Qt, RR);
    
    // L2-norm of manu_velocity
    norm += exact_velo.dot_product(exact_velo) * gwts;
  }

  return norm;
}

double POST_ERROR_C::get_manu_sol_u_errorH1(
      const double * const &sol_u,
      const double * const &sol_v,
      const double * const &sol_w,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &Qt,
      const double &RR )
{
  double errorH1 = 0.0;

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const std::vector<double> NR = element -> get_R( qua );
    const std::vector<double> NR_x = element -> get_dR_dx( qua );
    const std::vector<double> NR_y = element -> get_dR_dy( qua );
    const std::vector<double> NR_z = element -> get_dR_dz( qua );
    
    Vector_3 sol_velo (0.0, 0.0, 0.0);
    Vector_3 sol_grad_u (0.0, 0.0, 0.0);
    Vector_3 sol_grad_v (0.0, 0.0, 0.0);
    Vector_3 sol_grad_w (0.0, 0.0, 0.0);

    Vector_3 pt (0.0, 0.0, 0.0);

    for(int ii=0; ii < element->get_nLocBas(); ++ii)
    {
      pt.x() += ectrlPts_x[ii] * NR[ii];
      pt.y() += ectrlPts_y[ii] * NR[ii];
      pt.z() += ectrlPts_z[ii] * NR[ii];

      sol_velo.x() += sol_u[ii] * NR[ii];
      sol_velo.y() += sol_v[ii] * NR[ii];
      sol_velo.z() += sol_w[ii] * NR[ii];

      sol_grad_u.x() += sol_u[ii] * NR_x[ii];
      sol_grad_u.y() += sol_u[ii] * NR_y[ii];
      sol_grad_u.z() += sol_u[ii] * NR_z[ii];

      sol_grad_v.x() += sol_v[ii] * NR_x[ii];
      sol_grad_v.y() += sol_v[ii] * NR_y[ii];
      sol_grad_v.z() += sol_v[ii] * NR_z[ii];

      sol_grad_w.x() += sol_w[ii] * NR_x[ii];
      sol_grad_w.y() += sol_w[ii] * NR_y[ii];
      sol_grad_w.z() += sol_w[ii] * NR_z[ii];
    
    }

    const Vector_3 exa_velo = exact_velocity( pt, Qt, RR );
    const Vector_3 exa_grad_u = exact_grad_u( pt, Qt, RR );
    const Vector_3 exa_grad_v = exact_grad_v( pt, Qt, RR );
    const Vector_3 exa_grad_w = exact_grad_w( pt, Qt, RR );

    // H1 norm of error

    errorH1 += (sol_grad_u - exa_grad_u).dot_product( sol_grad_u - exa_grad_u) * gwts
	         + (sol_grad_v - exa_grad_v).dot_product( sol_grad_v - exa_grad_v) * gwts
	         + (sol_grad_w - exa_grad_w).dot_product( sol_grad_w - exa_grad_w) * gwts 
	         + (sol_velo - exa_velo).dot_product(sol_velo - exa_velo) * gwts;
  }

  return errorH1;
}

double POST_ERROR_C::get_manu_sol_u_H1(
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      const double &Qt,
      const double &RR )
{
  double norm = 0.0;

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const std::vector<double> NR = element -> get_R( qua );
    const std::vector<double> NR_x = element -> get_dR_dx( qua );
    const std::vector<double> NR_y = element -> get_dR_dy( qua );
    const std::vector<double> NR_z = element -> get_dR_dz( qua );
    
    Vector_3 sol_velo (0.0, 0.0, 0.0);
    Vector_3 sol_grad_u (0.0, 0.0, 0.0);
    Vector_3 sol_grad_v (0.0, 0.0, 0.0);
    Vector_3 sol_grad_w (0.0, 0.0, 0.0);

    Vector_3 pt (0.0, 0.0, 0.0);

    for(int ii=0; ii < element->get_nLocBas(); ++ii)
    {
      pt.x() += ectrlPts_x[ii] * NR[ii];
      pt.y() += ectrlPts_y[ii] * NR[ii];
      pt.z() += ectrlPts_z[ii] * NR[ii];
    }

    const Vector_3 exa_velo = exact_velocity( pt, Qt, RR );
    const Vector_3 exa_grad_u = exact_grad_u( pt, Qt, RR );
    const Vector_3 exa_grad_v = exact_grad_v( pt, Qt, RR );
    const Vector_3 exa_grad_w = exact_grad_w( pt, Qt, RR );

    // H1 norm of error

    norm += exa_grad_u.dot_product(exa_grad_u) * gwts
             + exa_grad_v.dot_product(exa_grad_v) * gwts
             + exa_grad_w.dot_product(exa_grad_w) * gwts 
             + exa_velo.dot_product(exa_velo) * gwts;
  }

  return norm;
}