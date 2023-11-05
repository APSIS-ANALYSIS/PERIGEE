#include "PLocAssem_FSI_Mesh_Elastostatic.hpp"

PLocAssem_FSI_Mesh_Elastostatic::PLocAssem_FSI_Mesh_Elastostatic(
    const double &in_mat_E, const double &in_mat_nu )
: E(in_mat_E), nu(in_mat_nu), lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  nLocBas(4), vec_size(12)
{
  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  Zero_Tangent_Residual();

  print_info();
}

PLocAssem_FSI_Mesh_Elastostatic::~PLocAssem_FSI_Mesh_Elastostatic()
{
  delete [] Tangent; delete [] Residual; Tangent = nullptr; Residual = nullptr;
}

void PLocAssem_FSI_Mesh_Elastostatic::print_info() const
{
  SYS_T::commPrint("  Three-dimensional Elastostatic equation: \n");
  SYS_T::commPrint("  Spatial: Galerkin Finite element \n");
  SYS_T::commPrint("  This solver is for mesh motion in the fluid sub-domain for FSI problems.\n");
  SYS_T::commPrint("  Young's Modulus E  = %e \n", E);
  SYS_T::commPrint("  Possion's ratio nu = %e \n", nu);
  SYS_T::commPrint("  Lame coeff lambda  = %e \n", lambda);
  SYS_T::commPrint("  Shear modulus mu   = %e \n", mu);
  SYS_T::commPrint("  Bulk modulus kappa = %e \n", kappa);
  SYS_T::commPrint("  Note: Element stiffening is applied. \n");
  SYS_T::print_sep_line();
}


void PLocAssem_FSI_Mesh_Elastostatic::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 0.0;
}

void PLocAssem_FSI_Mesh_Elastostatic::Zero_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
}

void PLocAssem_FSI_Mesh_Elastostatic::Assem_Estimate()
{
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 1.0;
}

void PLocAssem_FSI_Mesh_Elastostatic::get_currPts(
    const double * const &ept_x, const double * const &ept_y,
    const double * const &ept_z, const double * const &sol,
    double * const &cur_x, double * const &cur_y,
    double * const &cur_z )
{
  for(int ii=0; ii<nLocBas; ++ii)
  {
    cur_x[ii] = ept_x[ii] + sol[3*ii];
    cur_y[ii] = ept_y[ii] + sol[3*ii+1];
    cur_z[ii] = ept_z[ii] + sol[3*ii+2];
  }
}

void PLocAssem_FSI_Mesh_Elastostatic::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  double curPt_x[4], curPt_y[4], curPt_z[4];

  // vec_a passes the previous time step displacement.
  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, vec_a, curPt_x, curPt_y, curPt_z );

  const int nqp = quad -> get_num_quadPts();

  // Because the call of get_currPts, this basis is at tilde(x)
  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  const double l2mu = lambda + 2.0 * mu;

  Zero_Tangent_Residual();

  double dR_dx[4], dR_dy[4], dR_dz[4];

  for( int qua=0; qua<nqp; ++qua )
  {
    element->get_gradR(qua, dR_dx, dR_dy, dR_dz);
    const double detJac = element->get_detJac(qua);

    double ux = 0.0, uy = 0.0, uz = 0.0;
    double vx = 0.0, vy = 0.0, vz = 0.0;
    double wx = 0.0, wy = 0.0, wz = 0.0;

    for(int ii=0; ii<nLocBas; ++ii)
    {
      // Let vec_a be the previous disp, vec_b be the current disp
      ux += (vec_b[ii*3+0] - vec_a[ii*3+0]) * dR_dx[ii];
      uy += (vec_b[ii*3+0] - vec_a[ii*3+0]) * dR_dy[ii];
      uz += (vec_b[ii*3+0] - vec_a[ii*3+0]) * dR_dz[ii];

      vx += (vec_b[ii*3+1] - vec_a[ii*3+1]) * dR_dx[ii];
      vy += (vec_b[ii*3+1] - vec_a[ii*3+1]) * dR_dy[ii];
      vz += (vec_b[ii*3+1] - vec_a[ii*3+1]) * dR_dz[ii];

      wx += (vec_b[ii*3+2] - vec_a[ii*3+2]) * dR_dx[ii];
      wy += (vec_b[ii*3+2] - vec_a[ii*3+2]) * dR_dy[ii];
      wz += (vec_b[ii*3+2] - vec_a[ii*3+2]) * dR_dz[ii];
    }
    
    // element stiffening: we multiply the inverse of the element volume given
    // by detJac to enhance the mesh moving algorithm's robustness, as the small
    // element will be stiffer and thus is resistant to mesh distortion.
    // Reference: T.E. Tezduyar and A.A. Johnson (1994) CMAME 119 (1994) 73–94
    const double gwts = detJac * quad->get_qw(qua) / detJac;

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      Residual[3*A] += gwts * (
          NA_x * (l2mu * ux + lambda * vy + lambda * wz)
          + NA_y * mu * (uy + vx)
          + NA_z * mu * (uz + wx) );

      Residual[3*A+1] += gwts * (
          NA_x * mu * (uy + vx)
          + NA_y * (lambda * ux + l2mu * vy + lambda * wz)
          + NA_z * mu * (vz + wy) );

      Residual[3*A+2] += gwts * (
          NA_x * mu * (uz + wx)
          + NA_y * mu * (vz + wy)
          + NA_z * (lambda * ux + lambda * vy + l2mu * wz ) );

      for(int B=0; B<nLocBas; ++B)
      {
        const double NB_x = dR_dx[B], NB_y = dR_dy[B], NB_z = dR_dz[B];

        Tangent[ 3*nLocBas*(3*A+0) + 3*B + 0 ] += gwts * (l2mu * NA_x * NB_x + mu * NA_y * NB_y + mu * NA_z * NB_z );

        Tangent[ 3*nLocBas*(3*A+0) + 3*B + 1 ] += gwts * (lambda * NA_x * NB_y + mu * NA_y * NB_x);

        Tangent[ 3*nLocBas*(3*A+0) + 3*B + 2 ] += gwts * (lambda * NA_x * NB_z + mu * NA_z * NB_x);

        Tangent[ 3*nLocBas*(3*A+1) + 3*B + 0 ] += gwts * (mu * NA_x * NB_y + lambda * NA_y * NB_x);

        Tangent[ 3*nLocBas*(3*A+1) + 3*B + 1 ] += gwts * (mu * NA_x * NB_x + l2mu * NA_y * NB_y + mu * NA_z * NB_z);

        Tangent[ 3*nLocBas*(3*A+1) + 3*B + 2 ] += gwts * (lambda * NA_y * NB_z + mu * NA_z * NB_y);

        Tangent[ 3*nLocBas*(3*A+2) + 3*B + 0 ] += gwts * (lambda * NA_z * NB_x + mu * NA_x * NB_z);

        Tangent[ 3*nLocBas*(3*A+2) + 3*B + 1 ] += gwts * (lambda * NA_z * NB_y + mu * NA_y * NB_z);

        Tangent[ 3*nLocBas*(3*A+2) + 3*B + 2 ] += gwts * ( mu * NA_x * NB_x + mu * NA_y * NB_y + l2mu * NA_z * NB_z );
      }
    }
  }
}

// EOF
