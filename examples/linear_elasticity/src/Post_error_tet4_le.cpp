#include "Post_error_tet4_le.hpp"

void POST_T_TET4_LE::exact_disp( const double &x, const double &y, 
    const double &z, double &val_x, double &val_y, double &val_z )
{
  const double pi = MATH_T::PI;
  
  val_x = sin(pi*x) * sin(2*pi*y) * sin(pi*z);
  val_y = sin(2*pi*x) * sin(pi*y) * sin(pi*z);
  val_z = sin(pi*x) * sin(pi*y) * sin(2*pi*z);
}


void POST_T_TET4_LE::exact_3dmatrix( const double &x, const double &y,
          const double &z, Matrix_3x3 &val )
{
  const double pi = MATH_T::PI;
  
  // Copy these values from the analysis driver code
  const double lambda = 0.5769231;
  const double mu = 0.3846154;
  
  val(0) = lambda*(pi*cos(pi*x)*sin(2*pi*y)*sin(pi*z) + pi*cos(pi*y)*sin(2*pi*x)*sin(pi*z) + 2*pi*cos(2*pi*z)*sin(pi*x)*sin(pi*y)) + 2*mu*pi*cos(pi*x)*sin(2*pi*y)*sin(pi*z);

  val(1) = mu*(2*pi*cos(2*pi*x)*sin(pi*y)*sin(pi*z) + 2*pi*cos(2*pi*y)*sin(pi*x)*sin(pi*z));

  val(2) = mu*(pi*cos(pi*x)*sin(pi*y)*sin(2*pi*z) + pi*cos(pi*z)*sin(pi*x)*sin(2*pi*y));

  val(3) = val(1);

  val(4) = lambda*(pi*cos(pi*x)*sin(2*pi*y)*sin(pi*z) + pi*cos(pi*y)*sin(2*pi*x)*sin(pi*z) + 2*pi*cos(2*pi*z)*sin(pi*x)*sin(pi*y)) + 2*mu*pi*cos(pi*y)*sin(2*pi*x)*sin(pi*z);

  val(5) = mu*(pi*cos(pi*y)*sin(pi*x)*sin(2*pi*z) + pi*cos(pi*z)*sin(2*pi*x)*sin(pi*y));

  val(6) = val(2);

  val(7) = val(5);

  val(8) = lambda*(pi*cos(pi*x)*sin(2*pi*y)*sin(pi*z) + pi*cos(pi*y)*sin(2*pi*x)*sin(pi*z) + 2*pi*cos(2*pi*z)*sin(pi*x)*sin(pi*y)) + 4*mu*pi*cos(2*pi*z)*sin(pi*x)*sin(pi*y);

}


double POST_T_TET4_LE:: get_manu_disp_error(
    const double * const &solu,
    const double * const &solv,
    const double * const &solw,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const AInt_Weight * const &weight,
    double * const &R )
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

    exact_disp(coor_x, coor_y, coor_z, exa_u, exa_v, exa_w);

    error += (sol_u - exa_u) * (sol_u - exa_u) * gwts;
    error += (sol_v - exa_v) * (sol_v - exa_v) * gwts;
    error += (sol_w - exa_w) * (sol_w - exa_w) * gwts;
  }

  return error;
}


double POST_T_TET4_LE::get_manu_stress_error(
    const double * const &sol_x,
    const double * const &sol_y,
    const double * const &sol_z,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const AInt_Weight * const &weight,
    double * const &R,
    double * const &dR_dx,
    double * const &dR_dy,
    double * const &dR_dz )
{
  Matrix_3x3 sigma_exa, sigma_sol;
  double u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z;
  double coor_x, coor_y, coor_z, gwts;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  // Copy these values from the analysis driver code
  const double lambda = 0.5769231;
  const double mu = 0.3846154;

  double error = 0.0;
  
  for(int qua=0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    element->get_R(qua, R);
    element->get_gradR(qua, dR_dx, dR_dy, dR_dz);

    u_x = 0.0; u_y = 0.0; u_z = 0.0;
    v_x = 0.0; v_y = 0.0; v_z = 0.0;
    w_x = 0.0; w_y = 0.0; w_z = 0.0;

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];

      u_x += sol_x[ii] * dR_dx[ii];
      u_y += sol_x[ii] * dR_dy[ii];
      u_z += sol_x[ii] * dR_dz[ii];

      v_x += sol_y[ii] * dR_dx[ii];
      v_y += sol_y[ii] * dR_dy[ii];
      v_z += sol_y[ii] * dR_dz[ii];

      w_x += sol_z[ii] * dR_dx[ii];
      w_y += sol_z[ii] * dR_dy[ii];
      w_z += sol_z[ii] * dR_dz[ii];
    }
    
    const double divu = u_x + v_y + w_z;

    sigma_sol(0,0) = 2.0 * mu * u_x + lambda * divu;
    sigma_sol(0,1) = mu * ( u_y + v_x );
    sigma_sol(0,2) = mu * ( u_z + w_x );

    sigma_sol(1,0) = sigma_sol(0,1);
    sigma_sol(1,1) = 2.0 * mu * v_y + lambda * divu;
    sigma_sol(1,2) = mu * ( v_z + w_y );

    sigma_sol(2,0) = sigma_sol(0,2);
    sigma_sol(2,1) = sigma_sol(1,2);
    sigma_sol(2,2) = 2.0 * mu * w_z + lambda * divu;

    exact_3dmatrix(coor_x, coor_y, coor_z, sigma_exa);

    for(int ii=0; ii<9; ++ii)
      error += gwts * (sigma_sol(ii) - sigma_exa(ii)) * (sigma_sol(ii) - sigma_exa(ii));
  }

  return error;
}

// EOF
