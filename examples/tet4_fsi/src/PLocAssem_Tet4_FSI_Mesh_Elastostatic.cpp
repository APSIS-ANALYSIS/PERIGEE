#include "PLocAssem_Tet4_FSI_Mesh_Elastostatic.hpp"

PLocAssem_Tet4_FSI_Mesh_Elastostatic::PLocAssem_Tet4_FSI_Mesh_Elastostatic(
    const double &in_mat_E, const double &in_mat_nu )
: E(in_mat_E), nu(in_mat_nu), lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  nLocBas(4), dof_per_node(7), vec_size(12)
{
  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  Zero_Tangent_Residual();

  print_info();
}


PLocAssem_Tet4_FSI_Mesh_Elastostatic::~PLocAssem_Tet4_FSI_Mesh_Elastostatic()
{
  delete [] Tangent; delete [] Residual; Tangent = NULL; Residual = NULL;
}


void PLocAssem_Tet4_FSI_Mesh_Elastostatic::print_info() const
{
  SYS_T::print_sep_line();
  PetscPrintf(PETSC_COMM_WORLD, "  Three-dimensional Elastostatic equation: \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Spatial: Galerkin Finite element \n");
  PetscPrintf(PETSC_COMM_WORLD, "  This solver is for mesh motion in the fluid sub-domain for FSI problems.\n");
  PetscPrintf(PETSC_COMM_WORLD, "  Young's Modulus E  = %e \n", E);
  PetscPrintf(PETSC_COMM_WORLD, "  Possion's ratio nu = %e \n", nu);
  PetscPrintf(PETSC_COMM_WORLD, "  Lame coeff lambda  = %e \n", lambda);
  PetscPrintf(PETSC_COMM_WORLD, "  Shear modulus mu   = %e \n", mu);
  PetscPrintf(PETSC_COMM_WORLD, "  Bulk modulus kappa = %e \n", kappa);
  PetscPrintf(PETSC_COMM_WORLD, "  Note: Element stiffening is applied. \n");
  SYS_T::print_sep_line();
}


void PLocAssem_Tet4_FSI_Mesh_Elastostatic::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 0.0;
}


void PLocAssem_Tet4_FSI_Mesh_Elastostatic::Zero_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
}


void PLocAssem_Tet4_FSI_Mesh_Elastostatic::Assem_Estimate()
{
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 1.0;
}


void PLocAssem_Tet4_FSI_Mesh_Elastostatic::get_currPts(
    const double * const &ept_x, const double * const &ept_y,
    const double * const &ept_z, const double * const &sol )
{
  for(int ii=0; ii<nLocBas; ++ii)
  {
    curPt_x[ii] = ept_x[ii] + sol[7*ii];
    curPt_y[ii] = ept_y[ii] + sol[7*ii+1];
    curPt_z[ii] = ept_z[ii] + sol[7*ii+2];
  }
}


void PLocAssem_Tet4_FSI_Mesh_Elastostatic::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  // vec_a passes the previous time step displacement.
  get_currPts(eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z, vec_a );

  const int nqp = quad -> get_num_quadPts();

  // Because the call of get_currPts, this basis is at tilde(x)
  element->buildBasis( quad, curPt_x, curPt_y, curPt_z );

  int ii, jj, qua, A, B, index, ii7;
  double ux, uy, uz, vx, vy, vz, wx, wy, wz;
  double gwts, NA_x, NA_y, NA_z, NB_x, NB_y, NB_z, detJac;
  
  const double l2mu = lambda + 2.0 * mu;

  Zero_Tangent_Residual();

  Zero_Sub_Tan();

  for( qua=0; qua<nqp; ++qua )
  {
    element->get_gradR(qua, dR_dx, dR_dy, dR_dz);
    detJac = element->get_detJac(qua);

    ux = 0.0; uy = 0.0; uz = 0.0;
    vx = 0.0; vy = 0.0; vz = 0.0;
    wx = 0.0; wy = 0.0; wz = 0.0;

    for(ii=0; ii<nLocBas; ++ii)
    {
      ii7 = ii * 7;

      // Let vec_a be the previous disp, vec_b be the current disp
      ux += (vec_b[ii7+0] - vec_a[ii7+0]) * dR_dx[ii];
      uy += (vec_b[ii7+0] - vec_a[ii7+0]) * dR_dy[ii];
      uz += (vec_b[ii7+0] - vec_a[ii7+0]) * dR_dz[ii];

      vx += (vec_b[ii7+1] - vec_a[ii7+1]) * dR_dx[ii];
      vy += (vec_b[ii7+1] - vec_a[ii7+1]) * dR_dy[ii];
      vz += (vec_b[ii7+1] - vec_a[ii7+1]) * dR_dz[ii];

      wx += (vec_b[ii7+2] - vec_a[ii7+2]) * dR_dx[ii];
      wy += (vec_b[ii7+2] - vec_a[ii7+2]) * dR_dy[ii];
      wz += (vec_b[ii7+2] - vec_a[ii7+2]) * dR_dz[ii];
    }

    gwts = detJac * quad->get_qw(qua) / detJac; // element stiffening

    for(A=0; A<nLocBas; ++A)
    {
      NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

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

      for(B=0; B<nLocBas; ++B)
      {
        index = A * nLocBas + B;
        NB_x = dR_dx[B]; NB_y = dR_dy[B]; NB_z = dR_dz[B];

        Sub_Tan[0][index] += gwts * (l2mu * NA_x * NB_x
            + mu * NA_y * NB_y + mu * NA_z * NB_z );

        Sub_Tan[1][index] += gwts * (lambda * NA_x * NB_y + mu * NA_y * NB_x);

        Sub_Tan[2][index] += gwts * (lambda * NA_x * NB_z + mu * NA_z * NB_x);

        Sub_Tan[3][index] += gwts * (mu * NA_x * NB_y + lambda * NA_y * NB_x);

        Sub_Tan[4][index] += gwts * (mu * NA_x * NB_x + l2mu * NA_y * NB_y
            + mu * NA_z * NB_z);

        Sub_Tan[5][index] += gwts * (lambda * NA_y * NB_z + mu * NA_z * NB_y);

        Sub_Tan[6][index] += gwts * (lambda * NA_z * NB_x + mu * NA_x * NB_z);

        Sub_Tan[7][index] += gwts * (lambda * NA_z * NB_y + mu * NA_y * NB_z);

        Sub_Tan[8][index] += gwts * ( mu * NA_x * NB_x + mu * NA_y * NB_y
            + l2mu * NA_z * NB_z );
      }
    }
  }

  for(ii=0; ii<3; ++ii)
  {
    for(jj=0; jj<3; ++jj)
    {
      for(A=0; A<nLocBas; ++A)
      {
        for(B=0; B<nLocBas; ++B)
        {
          Tangent[ 3*nLocBas*(3*A+ii) + 3*B + jj ]
            = Sub_Tan[ii*3+jj][A*nLocBas + B];
        }
      }
    }
  }
}

// EOF
