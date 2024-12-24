#include "Post_error_elastodynamics.hpp"

Vector_3 POST_ERROR_E::exact_disp( const double &x, const double &y, const double &z, const double &time )
{
  const double ux = aa*x*sin(w*x)*sin(w*y)*sin(w*z)*time;
  const double uy = aa*y*sin(w*x)*sin(w*y)*sin(w*z)*time;
  const double uz = aa*z*sin(w*x)*sin(w*y)*sin(w*z)*time;
  return Vector_3(ux, uy, uz);
}

Tensor2_3D POST_ERROR_E::exact_grad_disp( const double &x, const double &y, const double &z, const double &time )
{
  const double ux_x = aa*time*sin(w*x)*sin(w*y)*sin(w*z) + aa*time*w*x*cos(w*x)*sin(w*y)*sin(w*z);
  const double ux_y = aa*time*w*x*cos(w*y)*sin(w*x)*sin(w*z);
  const double ux_z = aa*time*w*x*cos(w*z)*sin(w*x)*sin(w*y);
  const double uy_x = aa*time*w*y*cos(w*x)*sin(w*y)*sin(w*z);
  const double uy_y = aa*time*sin(w*x)*sin(w*y)*sin(w*z) + aa*time*w*y*cos(w*y)*sin(w*x)*sin(w*z);
  const double uy_z = aa*time*w*y*cos(w*z)*sin(w*x)*sin(w*y);
  const double uz_x = aa*time*w*z*cos(w*x)*sin(w*y)*sin(w*z);
  const double uz_y = aa*time*w*z*cos(w*y)*sin(w*x)*sin(w*z);
  const double uz_z = aa*time*sin(w*x)*sin(w*y)*sin(w*z) + aa*time*w*z*cos(w*z)*sin(w*x)*sin(w*y);
  return Tensor2_3D(ux_x, ux_y, ux_z, uy_x, uy_y, uy_z, uz_x, uz_y, uz_z);
}

double POST_ERROR_E::get_manu_sol_errorL2(
    const double &time,
    const double * const &solux,
    const double * const &soluy,
    const double * const &soluz,
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

    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;
    Vector_3 sol(0.0, 0.0, 0.0);

    for(int ii=0; ii<element->get_nLocBas(); ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];

      sol.x() += solux[ii] * R[ii];
      sol.y() += soluy[ii] * R[ii];
      sol.z() += soluz[ii] * R[ii];
    }

    const Vector_3 exact = exact_disp(coor_x, coor_y, coor_z, time );
    
    error += (sol.x() - exact.x()) * (sol.x() - exact.x()) * gwts;
    error += (sol.y() - exact.y()) * (sol.y() - exact.y()) * gwts;
    error += (sol.z() - exact.z()) * (sol.z() - exact.z()) * gwts;
  }

  return error;
}

double POST_ERROR_E::get_manu_sol_errorH1(
    const double &time,
    const double * const &solux,
    const double * const &soluy,
    const double * const &soluz,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const IQuadPts * const &quad )
{
  double errorH1 = 0.0;

  double * R = new double [element->get_nLocBas()];
  double * Rx = new double [element->get_nLocBas()];
  double * Ry = new double [element->get_nLocBas()];
  double * Rz = new double [element->get_nLocBas()];

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    element->get_R_gradR(qua, R, Rx, Ry, Rz);
    
    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;
    Vector_3 sol(0.0, 0.0, 0.0);
    Tensor2_3D grad_sol = Ten2::gen_zero();

    for(int ii=0; ii<element-> get_nLocBas(); ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];

      sol.x() += solux[ii] * R[ii];
      sol.y() += soluy[ii] * R[ii];
      sol.z() += soluz[ii] * R[ii];

      grad_sol(0,0) += solux[ii] * Rx[ii];
      grad_sol(0,1) += solux[ii] * Ry[ii];
      grad_sol(0,2) += solux[ii] * Rz[ii];

      grad_sol(1,0) += soluy[ii] * Rx[ii];
      grad_sol(1,1) += soluy[ii] * Ry[ii];
      grad_sol(1,2) += soluy[ii] * Rz[ii];

      grad_sol(2,0) += soluz[ii] * Rx[ii];
      grad_sol(2,1) += soluz[ii] * Ry[ii];
      grad_sol(2,2) += soluz[ii] * Rz[ii];
    }

    const Vector_3 exact = exact_disp(coor_x, coor_y, coor_z, time );
    const Tensor2_3D grad_exact = exact_grad_disp(coor_x, coor_y, coor_z, time);

    for(int ii=0; ii<9; ++ii)
    {
      errorH1 += (grad_sol(ii) - grad_exact(ii)) * (grad_sol(ii) - grad_exact(ii)) * gwts;
    }
    errorH1 += (sol.x() - exact.x()) * (sol.x() - exact.x()) * gwts;
    errorH1 += (sol.y() - exact.y()) * (sol.y() - exact.y()) * gwts;
    errorH1 += (sol.z() - exact.z()) * (sol.z() - exact.z()) * gwts;
  }

  delete [] R; R = nullptr;
  delete [] Rx; Rx = nullptr;
  delete [] Ry; Ry = nullptr;
  delete [] Rz; Rz = nullptr;

  return errorH1;
}

// EOF
