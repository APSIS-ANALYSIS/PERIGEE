#include "Post_error.hpp"

double POST_T::exact_scalar(const double &x, const double &y, const double &z,
		const double &t)
{
	const double pi = MATH_T::PI;
  
  //const double beta = 0.001 * pi;

  //const double u = t*t*( x * cos(beta*z) - y * sin(beta*z) - x );
  //const double u = t*t*( x * sin(beta*z) + y * cos(beta*z) - y );
  //const double u = 0.0;

  //const double J = 2*t*t*cos(beta*z) - 2*t*t*t*t*cos(beta*z) - 2*t*t + 2*t*t*t*t + 1.0;

  //const double p = 0.5 * 4.009433333333775e+06 * (J - 1.0 / J);

  return t*t*sin(pi*x)*sin(pi*y)*sin(pi*z);
}


double POST_T::exact_scalar(const double &x, const double &y, const double &t)
{
  const double pi = MATH_T::PI;
  
  const double theta = sin(pi*t)*sin(pi*x) * sin(pi*y);

  return theta;
}


void POST_T::exact_3dvector(const double &x, const double &y, const double &z,
		const double &t, double &val_x, double &val_y, double &val_z )
{
  const double pi = MATH_T::PI;

  val_x = pi*t*t*cos(pi*x)*sin(pi*y)*sin(pi*z);
  val_y = pi*t*t*cos(pi*y)*sin(pi*x)*sin(pi*z);
  val_z = pi*t*t*cos(pi*z)*sin(pi*x)*sin(pi*y);
}


void POST_T::exact_2dvector(const double &x, const double &y,
    const double &t, double &val_x, double &val_y )
{
  const double pi = MATH_T::PI;

  const double u_x = sin(pi*t)*cos(4*pi*x)*sin(4*pi*y);
  const double u_y = sin(pi*t)*sin(4*pi*x)*cos(4*pi*y);

  val_x = u_x;
  val_y = u_y;
}


double POST_T::get_manu_scalar_l2_error(
    const double * const &sol,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const AInt_Weight * const &weight,
    double * const &R,
    const double &curr )
{
  double error = 0.0;

  double exa_qua, sol_qua;
  double coor_x, coor_y, coor_z;

  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  for(int qua=0; qua<nqp; ++qua)
  {
    double gwts = element->get_detJac(qua) * weight->get_weight(qua);

    sol_qua = 0.0;
    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    element->get_R(qua, R);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];
      sol_qua += sol[ii] * R[ii];
    } 

    exa_qua = exact_scalar(coor_x, coor_y, coor_z, curr);

    error += (sol_qua - exa_qua)*(sol_qua - exa_qua) * gwts;
  }

  return error;
}


double POST_T::get_manu_scalar_l2_error(
    const double * const &sol,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const AInt_Weight * const &weight,
    double * const &R,
    const double &curr )
{
  double error = 0.0;

  double exa_qua, sol_qua;
  double coor_x, coor_y;

  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  for(int qua=0; qua<nqp; ++qua)
  {
    double gwts = element->get_detJac(qua) * weight->get_weight(qua);

    sol_qua = 0.0;
    coor_x = 0.0; coor_y = 0.0;

    element->get_R(qua, R);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      sol_qua += sol[ii] * R[ii];
    } 

    exa_qua = exact_scalar(coor_x, coor_y, curr);

    error += (sol_qua - exa_qua)*(sol_qua - exa_qua) * gwts;
  }

  return error;
}


double POST_T:: get_manu_scalar_h1_error(
    const double * const &solu,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const AInt_Weight * const &weight, 
    double * const &R,
    double * const &dR_dx,
    double * const &dR_dy,
    double * const &dR_dz,
    const double &curr )
{
  double error = 0.0;
  double exa_x, exa_y, exa_z, sol_x, sol_y, sol_z, exa, sol;
  double coor_x, coor_y, coor_z;
  double gwts;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  for(int qua=0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    sol = 0.0; sol_x = 0.0; sol_y = 0.0; sol_z = 0.0;
    exa = 0.0; exa_x = 0.0; exa_y = 0.0; exa_z = 0.0;
    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    element->get_R(qua, R);
    element->get_gradR(qua, dR_dx, dR_dy, dR_dz);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];
      sol    += solu[ii] * R[ii];
      sol_x  += solu[ii] * dR_dx[ii];
      sol_y  += solu[ii] * dR_dy[ii];
      sol_z  += solu[ii] * dR_dz[ii];
    }

    exa = exact_scalar(coor_x, coor_y, coor_z, curr);
    exact_3dvector(coor_x, coor_y, coor_z, curr, exa_x, exa_y, exa_z);

    error += (sol - exa) * (sol - exa) * gwts;
    error += (sol_x - exa_x) * (sol_x - exa_x) * gwts;
    error += (sol_y - exa_y) * (sol_y - exa_y) * gwts;
    error += (sol_z - exa_z) * (sol_z - exa_z) * gwts;
  }

  return error;
}


double POST_T:: get_manu_scalar_h1_error(
    const double * const &solu,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const AInt_Weight * const &weight, 
    double * const &R,
    double * const &dR_dx,
    double * const &dR_dy,
    const double &curr )
{
  double error = 0.0;
  double exa_x, exa_y, sol_x, sol_y, exa, sol;
  double coor_x, coor_y;
  double gwts;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  for(int qua=0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    sol = 0.0; sol_x = 0.0; sol_y = 0.0;
    exa = 0.0; exa_x = 0.0; exa_y = 0.0;
    coor_x = 0.0; coor_y = 0.0;

    element->get_R_gradR(qua, R, dR_dx, dR_dy);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      sol    += solu[ii] * R[ii];
      sol_x  += solu[ii] * dR_dx[ii];
      sol_y  += solu[ii] * dR_dy[ii];
    }

    exa = exact_scalar(coor_x, coor_y, curr);
    exact_2dvector(coor_x, coor_y, curr, exa_x, exa_y);

    error += (sol - exa) * (sol - exa) * gwts;
    error += (sol_x - exa_x) * (sol_x - exa_x) * gwts;
    error += (sol_y - exa_y) * (sol_y - exa_y) * gwts;
  }

  return error;
}

// EOF
