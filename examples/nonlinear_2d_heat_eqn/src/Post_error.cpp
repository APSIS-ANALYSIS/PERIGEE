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


double POST_T::get_INSK_mass( const double * const &sol,
    const FEAElement * const &element,
    const AInt_Weight * const &weight,
    double * const &R )
{
  double mass = 0.0; // mass in this element
  double gwts;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();	

  for(int qua=0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    element->get_R(qua,R);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      mass += sol[ii] * R[ii] * gwts;
    }
  }
  return mass;
}


double POST_T::get_SINSK_Lyapunov(
    const double &theta,
    const double &We,
    const double * const &sol_rho,
    const double * const &sol_u,
    const double * const &sol_v,
    const double * const &sol_w,
    const FEAElement * const &element,
    const AInt_Weight * const &weight,
    double * const &R,
    double * const &R_x,
    double * const &R_y,
    double * const &R_z )
{
  double lya = 0.0;
  double rho, u, v, w, rho_x, rho_y, rho_z;
  double gwts, lya_1, lya_2, lya_3;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  double theta_8_27 = theta * 8.0 / 27.0;
  double invWe = 1.0 / We;

  for(int qua=0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    element->get_R_gradR(qua, R, R_x, R_y, R_z);

    rho = 0.0; rho_x = 0.0; rho_y = 0.0; rho_z = 0.0;
    u = 0.0; v = 0.0; w = 0.0;

    for(int ii=0; ii<nLocBas; ++ii)
    {
      rho += sol_rho[ii] * R[ii];
      u   += sol_u[ii] * R[ii];
      v   += sol_v[ii] * R[ii];
      w   += sol_w[ii] * R[ii];
      rho_x += sol_rho[ii] * R_x[ii];
      rho_y += sol_rho[ii] * R_y[ii];
      rho_z += sol_rho[ii] * R_z[ii];
    }

    lya_1 = theta_8_27 * rho * log(rho/(1.0-rho)) - rho * rho;
    lya_2 = 0.5 * invWe * (rho_x * rho_x + rho_y * rho_y + rho_z * rho_z);
    lya_3 = 0.5 * rho * (u*u+v*v+w*w);

    lya += gwts * (lya_1 + lya_2 + lya_3);
  }

  return lya;
}


double POST_T::get_SINSK_free_energy(
    const double &theta,
    const double * const &sol_rho,
    const FEAElement * const &element,
    const AInt_Weight * const &weight,
    double * const &R )
{
  double fe = 0.0; // free energy scalar in this element
  double rho, gwts, lya_1;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  double theta_8_27 = theta * 8.0 / 27.0;

  for(int qua=0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    element->get_R(qua, R);
    rho = 0.0;

    for(int ii=0; ii<nLocBas; ++ii)
      rho += sol_rho[ii] * R[ii];

    lya_1 = theta_8_27 * rho * log(rho/(1.0-rho)) - rho * rho;

    fe += gwts * ( lya_1 );
  }

  return fe;
}



double POST_T::get_TNSK_uz_theta(
    const double * const &uz_theta,
    const double * const &inv_theta,
    const FEAElement * const &element,
    const AInt_Weight * const &weight,
    double * const &R )
{
  double ut = 0.0, uz, the, gwts;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  for(int qua=0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    element->get_R(qua, R);

    uz = 0.0; the = 0.0;
    for(int ii=0; ii<nLocBas; ++ii)
    {
      uz  += uz_theta[ii] * R[ii];
      the += inv_theta[ii] * R[ii];
    }

    ut   += gwts * uz / (the * the);
  }

  return ut;
}


double POST_T::get_volume( const FEAElement * const &element, 
    const AInt_Weight * const &weight )
{
  double area = 0.0, gwts;
  int nqp = weight->get_num();
  for(int qua =0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);
    area += gwts;
  }
  return area;
}



double POST_T::get_TNSK_rhou(
    const double * const &rho_vec,
    const double * const &ux_theta,
    const double * const &uy_theta,
    const double * const &uz_theta,
    const double * const &inv_theta,
    const FEAElement * const &element,
    const AInt_Weight * const &weight,
    double * const &R )
{
  double rhou2 = 0.0, rho, ux, uy, uz, theta, gwts;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  for(int qua = 0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);
    element->get_R(qua, R);

    rho = 0.0; ux = 0.0; uy = 0.0; uz = 0.0; theta = 0.0;

    for(int ii=0; ii<nLocBas; ++ii)
    {
      rho += rho_vec[ii] * R[ii];
      ux  += ux_theta[ii] * R[ii];
      uy  += uy_theta[ii] * R[ii];
      uz  += uz_theta[ii] * R[ii];
      theta += inv_theta[ii] * R[ii];
    }

    rhou2 += gwts * rho * rho * (ux*ux + uy*uy + uz*uz) / (theta * theta);
  }

  return sqrt(rhou2);
}

double POST_T::get_TNSK_rhou(
    const double * const &rho_vec,
    const double * const &ux_theta,
    const double * const &uy_theta,
    const double * const &inv_theta,
    const FEAElement * const &element,
    const AInt_Weight * const &weight,
    double * const &R )
{
  double rhou2 = 0.0, rho, ux, uy, theta, gwts;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  for(int qua = 0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);
    element->get_R(qua, R);

    rho = 0.0; ux = 0.0; uy = 0.0; theta = 0.0;

    // Interpolation at the quadrature points
    for(int ii=0; ii<nLocBas; ++ii)
    {
      rho += rho_vec[ii] * R[ii];
      ux  += ux_theta[ii] * R[ii];
      uy  += uy_theta[ii] * R[ii];
      theta += inv_theta[ii] * R[ii];
    }

    rhou2 += gwts * rho * rho * (ux*ux + uy*uy) / (theta * theta);
  }

  return sqrt(rhou2);
}

// EOF
