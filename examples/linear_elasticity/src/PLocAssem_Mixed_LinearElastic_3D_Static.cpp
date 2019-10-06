#include "PLocAssem_Mixed_LinearElastic_3D_Static.hpp"

PLocAssem_Mixed_LinearElastic_3D_Static::PLocAssem_Mixed_LinearElastic_3D_Static( 
    const double &in_mat_E, const double &in_mat_nu, const int &in_nlocbas,
    const int &in_nqp, const int &in_snlocbas )
: num_ebc_fun( 4 ), E( in_mat_E ), nu( in_mat_nu ), 
  lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  nLocBas( in_nlocbas ), dof_per_node(4),
  vec_size(nLocBas * dof_per_node), nqp(in_nqp), snLocBas(in_snlocbas)
{
  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  Zero_Tangent_Residual();

  R = new double [nLocBas];
  dR_dx = new double [nLocBas];
  dR_dy = new double [nLocBas];
  dR_dz = new double [nLocBas];

  Sub_Tan = new double * [16];

  for(int ii=0; ii<16; ++ii) Sub_Tan[ii] = new double [nLocBas * nLocBas];

  if( num_ebc_fun == 0 ) flist = NULL;
  else flist = new locassem_mle3d_funs [num_ebc_fun];

  flist[0] = &PLocAssem_Mixed_LinearElastic_3D_Static::get_lef_H;
  flist[1] = &PLocAssem_Mixed_LinearElastic_3D_Static::get_rig_H;
  flist[2] = &PLocAssem_Mixed_LinearElastic_3D_Static::get_fro_H;
  flist[3] = &PLocAssem_Mixed_LinearElastic_3D_Static::get_bac_H;

  print_info();
}


PLocAssem_Mixed_LinearElastic_3D_Static::~PLocAssem_Mixed_LinearElastic_3D_Static()
{
  delete [] Tangent; Tangent = NULL;
  delete [] Residual; Residual = NULL;
  delete [] R;
  delete [] dR_dx; delete [] dR_dy; delete [] dR_dz;
  for(int ii=0; ii<16; ++ii) delete [] Sub_Tan[ii];

  delete [] Sub_Tan;

  if(num_ebc_fun > 0) delete [] flist;
}


void PLocAssem_Mixed_LinearElastic_3D_Static::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Three-dimensional Linear Elastostatic model Mixed FEM formulation: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Spatial: Finite element \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Young's Modulus E  = %e \n", E);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Possion's ratio nu = %e \n", nu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Lame coeff lambda  = %e \n", lambda);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Shear modulus mu   = %e \n", mu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Bulk modulus kappa = %e \n", kappa);
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PLocAssem_Mixed_LinearElastic_3D_Static::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double inv_lambda = 1.0 / lambda;

  int ii, jj, A, B, qua, index, ii4;
  double ux, uy, uz, vx, vy, vz, wx, wy, wz, p;
  double fx, fy, fz, gwts, coor_x, coor_y, coor_z;
  double NA, NA_x, NA_y, NA_z, NB, NB_x, NB_y, NB_z;

  Zero_Tangent_Residual();

  Zero_Sub_Tan();

  for(qua=0; qua<nqp; ++qua)
  {
    element -> get_R_gradR(qua, R, dR_dx, dR_dy, dR_dz);

    ux = 0.0; uy = 0.0; uz = 0.0;
    vx = 0.0; vy = 0.0; vz = 0.0;
    wx = 0.0; wy = 0.0; wz = 0.0;
    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0; p = 0.0;
    for(ii=0; ii<nLocBas; ++ii)
    {
      ii4 = 4 * ii;
      p  += vec_a[ii4] * R[ii];

      ux += vec_a[ii4+1] * dR_dx[ii];
      uy += vec_a[ii4+1] * dR_dy[ii];
      uz += vec_a[ii4+1] * dR_dz[ii];

      vx += vec_a[ii4+2] * dR_dx[ii];
      vy += vec_a[ii4+2] * dR_dy[ii];
      vz += vec_a[ii4+2] * dR_dz[ii];

      wx += vec_a[ii4+3] * dR_dx[ii];
      wy += vec_a[ii4+3] * dR_dy[ii];
      wz += vec_a[ii4+3] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    gwts = element -> get_detJac(qua) * quad -> get_qw(qua);

    get_f(coor_x, coor_y, coor_z, fx, fy, fz);

    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A]; NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      Residual[4*A] -= gwts * ( NA * (ux + vy + wz) + NA * p * inv_lambda );

      Residual[4*A+1] += gwts * ( NA_x * 2.0 * mu * ux + NA_y * mu * (uy + vx) 
          + NA_z * mu * (uz + wx) - NA_x * p - NA * fx );

      Residual[4*A+2] += gwts * ( NA_x * mu * (uy+vx) + NA_y * 2.0 * mu * vy
          + NA_z * mu * (vz + wy) - NA_y * p - NA * fy );

      Residual[4*A+3] += gwts * ( NA_x * mu * (uz+wx) + NA_y * mu * (vz+wy)
          + NA_z * 2.0 * mu * wz - NA_z * p - NA * fz );

      for(B=0; B<nLocBas; ++B)
      {
        index = A * nLocBas + B;
        NB = dR_dx[B]; NB_x = dR_dx[B]; NB_y = dR_dy[B]; NB_z = dR_dz[B];

        Sub_Tan[0][index] -= gwts * NA * NB * inv_lambda;

        Sub_Tan[1][index] -= gwts * NA * NB_x;

        Sub_Tan[2][index] -= gwts * NA * NB_y;

        Sub_Tan[3][index] -= gwts * NA * NB_z;

        Sub_Tan[4][index] -= gwts * NA_x * NB;

        Sub_Tan[5][index] += gwts * ( NA_x * 2.0 * mu * NB_x
            + NA_y * mu * NB_y + NA_z * mu * NB_z );

        Sub_Tan[6][index] += gwts * NA_y * mu * NB_x;

        Sub_Tan[7][index] += gwts * NA_z * mu * NB_x;

        Sub_Tan[8][index] -= gwts * NA_y * NB;

        Sub_Tan[9][index] += gwts * NA_x * mu * NB_y;

        Sub_Tan[10][index] += gwts * ( NA_x * mu * NB_x + NA_y * 2.0 * mu * NB_y
            + NA_z * mu * NB_z );

        Sub_Tan[11][index] += gwts * NA_z * mu * NB_y;

        Sub_Tan[12][index] -= gwts * NA_z * NB;

        Sub_Tan[13][index] += gwts * NA_x * mu * NB_z;

        Sub_Tan[14][index] += gwts * NA_y * mu * NB_z;

        Sub_Tan[15][index] += gwts * ( NA_x * mu * NB_x + NA_y * mu * NB_y
            + NA_z * 2.0 * mu * NB_z );
      }
    }
  }

  for(ii=0; ii<4; ++ii)
  {
    for(jj=0; jj<4; ++jj)
    {
      for(A=0; A<nLocBas; ++A)
      {
        for(B=0; B<nLocBas; ++B)
        {
          Tangent[ 4*nLocBas*(4*A+ii) + 4*B + jj ]
            = Sub_Tan[ii*4+jj][A*nLocBas + B];
        }
      }
    }
  }
}


void PLocAssem_Mixed_LinearElastic_3D_Static::Assem_Residual_EBC(
    const int &ebc_id,
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  int ii, qua, A;
  double gwts, coor_x, coor_y, coor_z, gx, gy, gz, nx, ny, nz, surface_area;

  Zero_Residual();

  for(qua=0; qua<face_nqp; ++qua)
  {
    element -> get_R(qua, R);
    element -> get_2d_normal_out(qua, nx, ny, nz, surface_area);
    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    for(ii=0; ii<snLocBas; ++ii)
    {
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    get_ebc_fun( ebc_id, coor_x, coor_y, coor_z, nx, ny, nz,
        gx, gy, gz );

    gwts = surface_area * quad -> get_qw(qua);

    for(A=0; A<snLocBas; ++A)
    {
      Residual[4*A+1] -= gwts * R[A] * gx;
      Residual[4*A+2] -= gwts * R[A] * gy;
      Residual[4*A+3] -= gwts * R[A] * gz;
    }
  }
}


// EOF
