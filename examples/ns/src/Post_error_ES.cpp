#include "Post_error_ES.hpp"

Vector_3 POST_ERROR_C::exact_velocity ( const double &x, const double &y, const double &z, const double &t )
{   
    const double pi = MATH_T::PI; 
    
    const double u = 0.5 * (1 - std::cos(2.0 * pi * t)) * 2.0 / pi * (1 - (y * y + z * z));
    const double v = 0.0;
    const double w = 0.0;

    const Vector_3 velocity (u, v, w);    

    return velocity;
}

double  POST_ERROR_C::exact_pressure ( const Vector_3 &pt, const double &t )
{
    const double pi = MATH_T::PI; 

    const double x = pt.x();
    const double y = pt.y();
    const double z = pt.z();

    const double pressure = -100.0 * x * 0.5 * (1 - std::cos(2.0 * pi * t));

    return pressure;
}

Vector_3 POST_ERROR_C::exact_grad_u ( const double &x, const double &y, const double &z, const double &t )
{
    const double pi = MATH_T::PI; 
    const double a = pi/4.0;
    const double d = pi/2.0;
    const double nu = 0.1/1.0;
    
    const double ux = a*std::exp(-d*d*nu*t)*(a*std::sin(a*x + d*y)*std::exp(a*z) - a*std::sin(a*y + d*z)*std::exp(a*x));
    const double uy = -a*std::exp(-d*d*nu*t)*(a*std::cos(a*y + d*z)*std::exp(a*x) - d*std::sin(a*x + d*y)*std::exp(a*z));   
    const double uz = -a*std::exp(-d*d*nu*t)*(a*std::cos(a*x + d*y)*std::exp(a*z) + d*std::cos(a*y + d*z)*std::exp(a*x));      

    const Vector_3 grad_u (ux, uy, uz); 

    return grad_u;
}

Vector_3 POST_ERROR_C::exact_grad_v ( const double &x, const double &y, const double &z, const double &t )
{
    const double pi = MATH_T::PI; 
    const double a = pi/4.0;
    const double d = pi/2.0;
    const double nu = 0.1/1.0;
    
    const double vx = -a*std::exp(-d*d*nu*t)*(a*std::cos(a*y + d*z)*std::exp(a*x) + d*std::cos(a*z + d*x)*std::exp(a*y));
    const double vy = -a*std::exp(-d*d*nu*t)*(a*std::sin(a*z + d*x)*std::exp(a*y) - a*std::sin(a*y + d*z)*std::exp(a*x));   
    const double vz = -a*std::exp(-d*d*nu*t)*(a*std::cos(a*z + d*x)*std::exp(a*y) - d*std::sin(a*y + d*z)*std::exp(a*x));  

    const Vector_3 grad_v (vx, vy, vz);
    return grad_v;
}

Vector_3 POST_ERROR_C::exact_grad_w ( const double &x, const double &y, const double &z, const double &t )
{
  const double pi = MATH_T::PI; 
  const double a = pi/4.0;
  const double d = pi/2.0;
  const double nu = 0.1/1.0;

  const double wx = -a*std::exp(-d*d*nu*t)*(a*std::cos(a*x + d*y)*std::exp(a*z) - d*std::sin(a*z + d*x)*std::exp(a*y));
  const double wy = -a*std::exp(-d*d*nu*t)*(a*std::cos(a*z + d*x)*std::exp(a*y) + d*std::cos(a*x + d*y)*std::exp(a*z));
  const double wz = -a*std::exp(-d*d*nu*t)*(a*std::sin(a*x + d*y)*std::exp(a*z) - a*std::sin(a*z + d*x)*std::exp(a*y));
   
    Vector_3 grad_w (wx, wy, wz);
    return grad_w;
}

double POST_ERROR_C::get_manu_sol_p_L2(
      const double &time,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad )
      // const double &Qt,
      // const double &RR )
{
  double norm = 0.0;

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const std::vector<double> NR = element -> get_R( qua );

    Vector_3 pt (0.0, 0.0, 0.0);
    double sol_pres = 0.0;
    
    for(int ii=0; ii<element-> get_nLocBas(); ++ii)
    {
      pt.x() += ectrlPts_x[ii] * NR[ii];
      pt.y() += ectrlPts_y[ii] * NR[ii];
      pt.z() += ectrlPts_z[ii] * NR[ii];
    }

    const double exact_pres = exact_pressure( pt, time );
    
    // L2-norm of manu_velocity
    norm += exact_pres * exact_pres * gwts;
  }

  return norm;
}

double POST_ERROR_C::get_manu_sol_p_errorL2(
      const double &time,
      const double * const &sol_p,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad )
      // const double &Qt,
      // const double &RR )
{
  double error = 0.0;

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const std::vector<double> NR = element -> get_R( qua );

    Vector_3 pt (0.0, 0.0, 0.0);
    double sol_pres = 0.0;
    
    for(int ii=0; ii<element-> get_nLocBas(); ++ii)
    {
      pt.x() += ectrlPts_x[ii] * NR[ii];
      pt.y() += ectrlPts_y[ii] * NR[ii];
      pt.z() += ectrlPts_z[ii] * NR[ii];

      sol_pres += sol_p[ii] * NR[ii];
    }

    const double exact_pres = exact_pressure( pt, time );
    
    // L2-norm of error
    // (sol_velo - exact_velo).print();

    error += (sol_pres - exact_pres) * (sol_pres - exact_pres) * gwts;
  }

  return error;
}

double POST_ERROR_C::get_manu_sol_u_errorL2(
      const double &time,
      const double * const &sol_u,
      const double * const &sol_v,
      const double * const &sol_w,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad )
      // const double &Qt,
      // const double &RR )
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

    const Vector_3 exact_velo = exact_velocity( pt.x(), pt.y(), pt.z(), time );
    
    // L2-norm of error
    // (sol_velo - exact_velo).print();

    error += (sol_velo - exact_velo).dot_product(sol_velo - exact_velo) * gwts;
  }

  return error;
}

double POST_ERROR_C::get_manu_sol_u_L2(
      const double &time,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad )
      // const double &Qt,
      // const double &RR )
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

    const Vector_3 exact_velo = exact_velocity( pt.x(), pt.y(), pt.z(), time );
    
    // L2-norm of manu_velocity
    norm += exact_velo.dot_product(exact_velo) * gwts;
  }

  return norm;
}

double POST_ERROR_C::get_manu_sol_u_errorH1(
      const double &time,
      const double * const &sol_u,
      const double * const &sol_v,
      const double * const &sol_w,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad )
      // const double &Qt,
      // const double &RR )
{
  double errorH1 = 0.0;

  const int nLocBas = element->get_nLocBas();

  std::vector<double> NR(nLocBas, 0.0), NR_x(nLocBas, 0.0), NR_y(nLocBas, 0.0), NR_z(nLocBas, 0.0);

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    element->get_R_gradR( qua, &NR[0], &NR_x[0], &NR_y[0], &NR_z[0] );
    
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

    const Vector_3 exa_velo = exact_velocity( pt.x(), pt.y(), pt.z(), time );
    const Vector_3 exa_grad_u = exact_grad_u( pt.x(), pt.y(), pt.z(), time );
    const Vector_3 exa_grad_v = exact_grad_v( pt.x(), pt.y(), pt.z(), time );
    const Vector_3 exa_grad_w = exact_grad_w( pt.x(), pt.y(), pt.z(), time );

    // H1 norm of error

    errorH1 += (sol_grad_u - exa_grad_u).dot_product( sol_grad_u - exa_grad_u) * gwts
	         + (sol_grad_v - exa_grad_v).dot_product( sol_grad_v - exa_grad_v) * gwts
	         + (sol_grad_w - exa_grad_w).dot_product( sol_grad_w - exa_grad_w) * gwts 
	         + (sol_velo - exa_velo).dot_product(sol_velo - exa_velo) * gwts;
  }

  return errorH1;
}

double POST_ERROR_C::get_manu_sol_u_H1(
      const double &time,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad )
      // const double &Qt,
      // const double &RR )
{
  double norm = 0.0;

  const int nLocBas = element->get_nLocBas();

  std::vector<double> NR(nLocBas, 0.0), NR_x(nLocBas, 0.0), NR_y(nLocBas, 0.0), NR_z(nLocBas, 0.0);

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    element->get_R_gradR( qua, &NR[0], &NR_x[0], &NR_y[0], &NR_z[0] );
    
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

    const Vector_3 exa_velo = exact_velocity( pt.x(), pt.y(), pt.z(), time );
    const Vector_3 exa_grad_u = exact_grad_u( pt.x(), pt.y(), pt.z(), time );
    const Vector_3 exa_grad_v = exact_grad_v( pt.x(), pt.y(), pt.z(), time );
    const Vector_3 exa_grad_w = exact_grad_w( pt.x(), pt.y(), pt.z(), time );


    // H1 norm of error

    norm += exa_grad_u.dot_product(exa_grad_u) * gwts
             + exa_grad_v.dot_product(exa_grad_v) * gwts
             + exa_grad_w.dot_product(exa_grad_w) * gwts 
             + exa_velo.dot_product(exa_velo) * gwts;
  }

  return norm;
}

//EOF
