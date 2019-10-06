#include "Post_error_tet4_ns.hpp"

double POST_T_TET4_NS::exact_pres(const double &x, const double &y, 
    const double &z, const double &t)
{
  const double pi = MATH_T::PI;
  return t*t*sin(pi*x)*sin(pi*y)*sin(pi*z);
} 


void POST_T_TET4_NS::exact_grad_pres(const double &x, const double &y, 
    const double &z, const double &t,
    double &val_x, double &val_y, double &val_z)
{
  const double pi = MATH_T::PI;
  val_x = t*t*pi*cos(pi*x)*sin(pi*y)*sin(pi*z);
  val_y = t*t*pi*sin(pi*x)*cos(pi*y)*sin(pi*z);
  val_z = t*t*pi*sin(pi*x)*sin(pi*y)*cos(pi*z);
}



void POST_T_TET4_NS::exact_velo(const double &x, const double &y, 
    const double &z, const double &t, 
    double &val_x, double &val_y, double &val_z)
{
  val_x = 2*t*t*x*x*y*z*(x-1)*(x-1)*(z-1)*(2*y*y-3*y + 1); 
  val_y = -4*t*t*x*y*y*z*(x-z)*(x-1)*(y-1)*(y-1)*(z-1);
  val_z = -2*t*t*x*y*z*z*(x-1)*(z-1)*(z-1)*(2*y*y-3*y+1);
}


double POST_T_TET4_NS::get_pres_l2_error(
    const double * const &sol,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const AInt_Weight * const &weight,
    double * const &R,
    const double &curr )
{
  double exa_qua, sol_qua;
  double coor_x, coor_y, coor_z;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();
  double error = 0.0;
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
    exa_qua = exact_pres(coor_x, coor_y, coor_z, curr);
    error += (sol_qua - exa_qua)*(sol_qua - exa_qua) * gwts;
  }

  return error;
}


double POST_T_TET4_NS::get_pres_h1_error(
    const double * const &sol,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const AInt_Weight * const &weight,
    double * const &R,
    double * const &Rx,
    double * const &Ry,
    double * const &Rz,
    const double &curr )
{
  double exa_x, exa_y, exa_z, sol_x, sol_y, sol_z;
  double coor_x, coor_y, coor_z;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();
  double error = 0.0;
  for(int qua=0; qua<nqp; ++qua)
  {
    double gwts = element->get_detJac(qua) * weight->get_weight(qua);
    sol_x = 0.0; sol_y = 0.0; sol_z = 0.0;
    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    element->get_R_gradR(qua, R, Rx, Ry, Rz);
    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];
      sol_x  += sol[ii] * Rx[ii];
      sol_y  += sol[ii] * Ry[ii];
      sol_z  += sol[ii] * Rz[ii];
    }
    exact_grad_pres(coor_x, coor_y, coor_z, curr, exa_x,
        exa_y, exa_z);
    error += (sol_x - exa_x)*(sol_x - exa_x) * gwts;
    error += (sol_y - exa_y)*(sol_y - exa_y) * gwts;
    error += (sol_z - exa_z)*(sol_z - exa_z) * gwts;
  }

  return error;
}


double POST_T_TET4_NS::get_velo_l2_error(
    const double * const &solu,
    const double * const &solv,
    const double * const &solw,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const AInt_Weight * const &weight,
    double * const &R,
    const double &curr )
{
  double error = 0.0;
  double exa_u, exa_v, exa_w, sol_u, sol_v, sol_w;
  double coor_x, coor_y, coor_z, gwts;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();
  for(int qua=0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    sol_u = 0.0; sol_v = 0.0; sol_w = 0.0;
    exa_u = 0.0; exa_v = 0.0; exa_w = 0.0;
    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    element->get_R(qua, R);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];
      sol_u  += solu[ii] * R[ii];
      sol_v  += solv[ii] * R[ii];
      sol_w  += solw[ii] * R[ii];
    }
    exact_velo(coor_x, coor_y, coor_z, curr, exa_u, exa_v, exa_w);

    error += (sol_u - exa_u) * (sol_u - exa_u) * gwts;
    error += (sol_v - exa_v) * (sol_v - exa_v) * gwts;
    error += (sol_w - exa_w) * (sol_w - exa_w) * gwts;
  }

  return error;
}


// EOF
