#include "Post_error_elastodynamics.hpp"

Vector_3 POST_ERROR_E::exact_disp( const double &x, const double &y, const double &z, const double &time )
{
  return Vector_3(0.0, 0.0, 0.0);
}

Tensor2_3D POST_ERROR_E::exact_grad_disp( const double &x, const double &y, const double &z, const double &time )
{
  return Tensor2_3D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}

double POST_ERROR_E::get_manu_sol_error(
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

    for(int ii=0; ii<element-> get_nLocBas(); ++ii)
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

  for(int qua=0; qua<quad->get_num_quadPts(); ++qua)
  {
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const std::vector<double> R = element -> get_R( qua );
    const std::vector<double> Rx = element -> get_dR_dx( qua );
    const std::vector<double> Ry = element -> get_dR_dy( qua );
    const std::vector<double> Rz = element -> get_dR_dz( qua );
    
    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;
    Vector_3 sol(0.0, 0.0, 0.0);
    Tensor2_3D grad_sol; grad_sol.gen_zero();

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

  return errorH1;
}

// EOF
