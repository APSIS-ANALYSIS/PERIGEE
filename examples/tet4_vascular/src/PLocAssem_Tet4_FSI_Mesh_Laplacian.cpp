#include "PLocAssem_Tet4_FSI_Mesh_Laplacian.hpp"

PLocAssem_Tet4_FSI_Mesh_Laplacian::PLocAssem_Tet4_FSI_Mesh_Laplacian(
    const int &in_nlocbas, const int &in_snlocbas )
: num_ebc_fun(0), nLocBas(4), dof_per_node(7), vec_size(12),
  snLocBas(in_snlocbas)
{
  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  Zero_Tangent_Residual();

  if( num_ebc_fun == 0 ) flist = NULL;
  else flist = new locassem_fsi_mesh_lap_funs [num_ebc_fun];

  print_info();
}


PLocAssem_Tet4_FSI_Mesh_Laplacian::~PLocAssem_Tet4_FSI_Mesh_Laplacian()
{
  delete [] Tangent; delete [] Residual; Tangent = NULL; Residual = NULL;
  if(num_ebc_fun > 0) delete [] flist;
}


void PLocAssem_Tet4_FSI_Mesh_Laplacian::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, " Three-dimensional Laplacian equation: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t Spatial: Galerkin Finite element \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t This solver is for mesh motion in the fluid sub-domain for FSI problems.\n");
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PLocAssem_Tet4_FSI_Mesh_Laplacian::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 0.0;
}



void PLocAssem_Tet4_FSI_Mesh_Laplacian::Zero_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
}



void PLocAssem_Tet4_FSI_Mesh_Laplacian::Assem_Estimate()
{
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 1.0;
}


void PLocAssem_Tet4_FSI_Mesh_Laplacian::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  const int nqp = quad -> get_num_quadPts();

  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y,
      eleCtrlPts_z );

  int ii, jj, qua, A, B, index, ii7;
  double ux, uy, uz, vx, vy, vz, wx, wy, wz;
  double gwts, NA_x, NA_y, NA_z, NB_x, NB_y, NB_z;

  Zero_Tangent_Residual();

  Zero_Sub_Tan();

  for(qua=0; qua<nqp; ++qua)
  {
    element->get_gradR(qua, dR_dx, dR_dy, dR_dz);
    
    ux = 0.0; uy = 0.0; uz = 0.0;
    vx = 0.0; vy = 0.0; vz = 0.0;
    wx = 0.0; wy = 0.0; wz = 0.0;

    for(ii=0; ii<nLocBas; ++ii)
    {
      ii7 = ii * 7;

      ux += vec_b[ii7+0] * dR_dx[ii];
      uy += vec_b[ii7+0] * dR_dy[ii];
      uz += vec_b[ii7+0] * dR_dz[ii];

      vx += vec_b[ii7+1] * dR_dx[ii];
      vy += vec_b[ii7+1] * dR_dy[ii];
      vz += vec_b[ii7+1] * dR_dz[ii];

      wx += vec_b[ii7+2] * dR_dx[ii];
      wy += vec_b[ii7+2] * dR_dy[ii];
      wz += vec_b[ii7+2] * dR_dz[ii];
    }
    
    gwts = element->get_detJac(qua) * quad->get_qw(qua);

    for(A=0; A<nLocBas; ++A)
    {
      NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      Residual[3*A] += gwts * ( NA_x * ux + NA_y * uy + NA_z * uz );
      
      Residual[3*A+1] += gwts * ( NA_x * vx + NA_y * vy + NA_z * vz );

      Residual[3*A+2] += gwts * ( NA_x * wx + NA_y * wy + NA_z * wz );

      for(B=0; B<nLocBas; ++B)
      {
        index = A * nLocBas + B;
        NB_x = dR_dx[B]; NB_y = dR_dy[B]; NB_z = dR_dz[B];

        Sub_Tan[0][index] += gwts * (NA_x * NB_x + NA_y * NB_y + NA_z * NB_z);

        Sub_Tan[4][index] += gwts * (NA_x * NB_x + NA_y * NB_y + NA_z * NB_z);

        Sub_Tan[8][index] += gwts * (NA_x * NB_x + NA_y * NB_y + NA_z * NB_z);
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


void PLocAssem_Tet4_FSI_Mesh_Laplacian::Assem_Residual(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  const int nqp = quad -> get_num_quadPts();

  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  int ii, qua, A, ii7;
  double ux, uy, uz, vx, vy, vz, wx, wy, wz;
  double gwts, NA_x, NA_y, NA_z;

  double coor_x, coor_y, coor_z, fx, fy, fz;

  Zero_Residual();

  const double curr = time;

  for(qua=0; qua<nqp; ++qua)
  {
    element->get_R_gradR(qua, R, dR_dx, dR_dy, dR_dz);
    
    ux = 0.0; uy = 0.0; uz = 0.0;
    vx = 0.0; vy = 0.0; vz = 0.0;
    wx = 0.0; wy = 0.0; wz = 0.0;

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    for(ii=0; ii<nLocBas; ++ii)
    {
      ii7 = ii * 7;

      ux += vec_b[ii7+0] * dR_dx[ii];
      uy += vec_b[ii7+0] * dR_dy[ii];
      uz += vec_b[ii7+0] * dR_dz[ii];

      vx += vec_b[ii7+1] * dR_dx[ii];
      vy += vec_b[ii7+1] * dR_dy[ii];
      vz += vec_b[ii7+1] * dR_dz[ii];

      wx += vec_b[ii7+2] * dR_dx[ii];
      wy += vec_b[ii7+2] * dR_dy[ii];
      wz += vec_b[ii7+2] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }
    
    gwts = element->get_detJac(qua) * quad->get_qw(qua);

    get_f(coor_x, coor_y, coor_z, curr, fx, fy, fz);

    for(A=0; A<nLocBas; ++A)
    {
      NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      Residual[3*A] += gwts * ( NA_x * ux + NA_y * uy + NA_z * uz - R[A] * fx );
      
      Residual[3*A+1] += gwts * ( NA_x * vx + NA_y * vy + NA_z * vz - R[A] * fy );

      Residual[3*A+2] += gwts * ( NA_x * wx + NA_y * wy + NA_z * wz - R[A] * fz );
    }
  }
}

// EOF
