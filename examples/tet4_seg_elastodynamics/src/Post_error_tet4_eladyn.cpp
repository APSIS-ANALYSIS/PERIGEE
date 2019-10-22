#include "Post_error_tet4_eladyn.hpp"

double POST_T_TET4_ELADYN::exact_scalar(const double &x, const double &y, 
    const double &z, const double &t)
{
	const double pi = MATH_T::PI;
  const double beta = 0.5 * pi;

  //const double J = 2*t*t*cos(beta*z) - 2*t*t*t*t*cos(beta*z) - 2*t*t + 2*t*t*t*t + 1.0;

  //const double kappa = 1.111111e2; 
  //const double p = (-0.5) * kappa * (J - 1.0 / J);
  //const double p = alpha * t * t * ( sin(beta*x) * sin(beta*y) * sin(beta*z) ); 
  
  //const double beta = 0.2 * pi;
  const double p = sin(beta*x)*sin(beta*y)*sin(beta*z)*t*t;

  return p;
}


void POST_T_TET4_ELADYN::exact_disp(const double &x, const double &y, 
    const double &z, const double &t, 
    double &val_x, double &val_y, double &val_z )
{
  const double pi = MATH_T::PI;
  //const double beta = 0.001* pi;
  //val_x = t*t*( x * cos(beta*z) - y * sin(beta*z) - x );
  //val_y = t*t*( x * sin(beta*z) + y * cos(beta*z) - y );
  //val_z = 0.0; 
  
  const double gamma = 0.1 * pi;
  val_x = t*t*sin(gamma*z)*sin(gamma*y);
  val_y = 0.0; 
  val_z = 0.0;
}


void POST_T_TET4_ELADYN::exact_velo(const double &x, const double &y, 
    const double &z, const double &t, 
    double &val_x, double &val_y, double &val_z )
{
  const double pi = MATH_T::PI;
  //const double beta = 0.001* pi;
  //val_x = 2.0*t*( x * cos(beta*z) - y * sin(beta*z) - x );
  //val_y = 2.0*t*( x * sin(beta*z) + y * cos(beta*z) - y );
  //val_z = 0.0;
  
  const double gamma = 0.1 * pi;
  val_x = 2.0*t*sin(gamma*z)*sin(gamma*y);
  val_y = 0.0; 
  val_z = 0.0;
}


void POST_T_TET4_ELADYN::exact_3dmatrix( const double &x, const double &y,
    const double &z, const double &t, Matrix_3x3 &val )
{
  const double pi = MATH_T::PI;
  const double beta = 0.001 * pi;
  Matrix_3x3 F;

  F(0,0) = (cos(beta*z) - 1)*t*t + 1.0;
  F(0,1) = -t*t*sin(beta*z);
  F(0,2) = -beta*t*t*(y*cos(beta*z) + x*sin(beta*z));

  F(1,0) = t*t*sin(beta*z);
  F(1,1) = (cos(beta*z) - 1)*t*t + 1;
  F(1,2) = beta*t*t*(x*cos(beta*z) - y*sin(beta*z));

  F(2,0) = 0.0;
  F(2,1) = 0.0;
  F(2,2) = 1.0;

  for(int ii=0; ii<9; ++ii) val(ii) = F(ii);
}


double POST_T_TET4_ELADYN::get_manu_scalar_l2_error(
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

    exa_qua = exact_scalar(coor_x, coor_y, coor_z, curr);

    error += (sol_qua - exa_qua)*(sol_qua - exa_qua) * gwts;
  }

  return error;
}


double POST_T_TET4_ELADYN:: get_manu_disp_error(
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

    exact_disp(coor_x, coor_y, coor_z, curr, exa_u, exa_v, exa_w);

    error += (sol_u - exa_u) * (sol_u - exa_u) * gwts;
    error += (sol_v - exa_v) * (sol_v - exa_v) * gwts;
    error += (sol_w - exa_w) * (sol_w - exa_w) * gwts;
  }

  return error;
}


double POST_T_TET4_ELADYN:: get_manu_velo_error(
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


double POST_T_TET4_ELADYN::get_manu_matrix_l2_error(
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
    double * const &dR_dz,
    const double &curr )
{
  double error = 0.0;
  Matrix_3x3 Fexa, Fsol;
  double coor_x, coor_y, coor_z, gwts;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  for(int qua=0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    element->get_R(qua, R);
    element->get_gradR(qua, dR_dx, dR_dy, dR_dz);

    Fsol.gen_id();

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];

      Fsol(0,0) += sol_x[ii] * dR_dx[ii];
      Fsol(0,1) += sol_x[ii] * dR_dy[ii];
      Fsol(0,2) += sol_x[ii] * dR_dz[ii];

      Fsol(1,0) += sol_y[ii] * dR_dx[ii];
      Fsol(1,1) += sol_y[ii] * dR_dy[ii];
      Fsol(1,2) += sol_y[ii] * dR_dz[ii];

      Fsol(2,0) += sol_z[ii] * dR_dx[ii];
      Fsol(2,1) += sol_z[ii] * dR_dy[ii];
      Fsol(2,2) += sol_z[ii] * dR_dz[ii]; 
    }

    exact_3dmatrix(coor_x, coor_y, coor_z, curr, Fexa);

    for(int ii=0; ii<9; ++ii)
      error += gwts * (Fsol(ii) - Fexa(ii)) * (Fsol(ii) - Fexa(ii));
  }

  return error;
}


double POST_T_TET4_ELADYN::get_manu_sigma_l2_error(
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
    double * const &dR_dz,
    const double &curr )
{
  double error = 0.0;
  Matrix_3x3 Fexa, Fsol;
  double coor_x, coor_y, coor_z, gwts;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  for(int qua=0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    element->get_R(qua, R);
    element->get_gradR(qua, dR_dx, dR_dy, dR_dz);

    Fsol.gen_id();

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];

      Fsol(0,0) += sol_x[ii] * dR_dx[ii];
      Fsol(0,1) += sol_x[ii] * dR_dy[ii];
      Fsol(0,2) += sol_x[ii] * dR_dz[ii];

      Fsol(1,0) += sol_y[ii] * dR_dx[ii];
      Fsol(1,1) += sol_y[ii] * dR_dy[ii];
      Fsol(1,2) += sol_y[ii] * dR_dz[ii];

      Fsol(2,0) += sol_z[ii] * dR_dx[ii];
      Fsol(2,1) += sol_z[ii] * dR_dy[ii];
      Fsol(2,2) += sol_z[ii] * dR_dz[ii]; 
    }

    exact_3dmatrix(coor_x, coor_y, coor_z, curr, Fexa);

    Matrix_3x3 bsol; bsol.MatMultTransposeRight( Fsol );
    Matrix_3x3 bexa; bexa.MatMultTransposeRight( Fexa );

    double je = Fexa.det();
    double js = Fsol.det();

    double trbs = bsol.tr();
    double trbe = bexa.tr();

    const double mu = 3.703704e+01;

    Matrix_3x3 sigma_e; 
    sigma_e.scale(-1.0 * mu * std::pow(je, -2.0/3.0) * trbe / (3.0 * je));

    Matrix_3x3 sigma_s; 
    sigma_s.scale(-1.0 * mu * std::pow(js, -2.0/3.0) * trbs / (3.0 * js));

    sigma_e.AXPY( std::pow(je, -2.0/3.0) * mu / je, bexa );
    sigma_s.AXPY( std::pow(js, -2.0/3.0) * mu / js, bsol );

    for(int ii=0; ii<9; ++ii)
      error += gwts * (sigma_e(ii) - sigma_s(ii)) * (sigma_e(ii) - sigma_s(ii));
  }

  return error;
}


// EOF
