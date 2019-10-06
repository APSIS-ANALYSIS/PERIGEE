#include "PLocAssem_LinearElastic_3D_Static.hpp"

PLocAssem_LinearElastic_3D_Static::PLocAssem_LinearElastic_3D_Static( 
    const double &in_mat_E, const double &in_mat_nu,
    const int &in_nlocbas, const int &in_nqp, const int &in_snlocbas )
: num_ebc_fun( 4 ),
  E( in_mat_E ), nu( in_mat_nu ), lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ), 
  nLocBas( in_nlocbas ), dof_per_node(3),
  vec_size(nLocBas * dof_per_node), nqp(in_nqp), snLocBas(in_snlocbas)
{
  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  Zero_Tangent_Residual();

  R = new double [nLocBas];
  dR_dx = new double [nLocBas];
  dR_dy = new double [nLocBas];
  dR_dz = new double [nLocBas];

  Sub_Tan = new double * [9];

  for(int ii=0; ii<9; ++ii) Sub_Tan[ii] = new double [nLocBas * nLocBas];

  if( num_ebc_fun == 0 ) flist = NULL;
  else flist = new locassem_le3d_fem_funs [num_ebc_fun];

  flist[0] = &PLocAssem_LinearElastic_3D_Static::get_lef_H;
  flist[1] = &PLocAssem_LinearElastic_3D_Static::get_rig_H;
  flist[2] = &PLocAssem_LinearElastic_3D_Static::get_fro_H;
  flist[3] = &PLocAssem_LinearElastic_3D_Static::get_bac_H;

  print_info();
}


PLocAssem_LinearElastic_3D_Static::~PLocAssem_LinearElastic_3D_Static()
{
  delete [] Tangent; Tangent = NULL;
  delete [] Residual; Residual = NULL;
  delete [] R;
  delete [] dR_dx; delete [] dR_dy; delete [] dR_dz;
  for(int ii=0; ii<9; ++ii) delete [] Sub_Tan[ii];

  delete [] Sub_Tan;

  if(num_ebc_fun > 0) delete [] flist;
}


void PLocAssem_LinearElastic_3D_Static::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Three-dimensional Linear Elastostatic model FEM formulation: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Spatial: Finite element \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Young's Modulus E  = %e \n", E);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Possion's ratio nu = %e \n", nu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Lame coeff lambda  = %e \n", lambda);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Shear modulus mu   = %e \n", mu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Bulk modulus kappa = %e \n", kappa);
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PLocAssem_LinearElastic_3D_Static::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 0.0;
}


void PLocAssem_LinearElastic_3D_Static::Zero_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
}


void PLocAssem_LinearElastic_3D_Static::Assem_Estimate()
{
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 1.0;
}


void PLocAssem_LinearElastic_3D_Static::Assem_Tangent_Residual(
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
  
  const double l2mu = lambda + 2.0 * mu;

  int ii, jj, A, B, qua, index, ii3;
  double ux, uy, uz, vx, vy, vz, wx, wy, wz;
  double fx, fy, fz, gwts, coor_x, coor_y, coor_z;
  double NA, NA_x, NA_y, NA_z, NB_x, NB_y, NB_z;

  Zero_Tangent_Residual();

  Zero_Sub_Tan();

  for(qua=0; qua<nqp; ++qua)
  {
    element->get_R_gradR(qua, R, dR_dx, dR_dy, dR_dz);
  
    ux = 0.0; uy = 0.0; uz = 0.0;
    vx = 0.0; vy = 0.0; vz = 0.0;
    wx = 0.0; wy = 0.0; wz = 0.0; 
    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0; 
    for(ii=0; ii<nLocBas; ++ii)
    {
      ii3 = 3 * ii;
      ux += vec_a[ii3] * dR_dx[ii];
      uy += vec_a[ii3] * dR_dy[ii];
      uz += vec_a[ii3] * dR_dz[ii];

      vx += vec_a[ii3+1] * dR_dx[ii];
      vy += vec_a[ii3+1] * dR_dy[ii];
      vz += vec_a[ii3+1] * dR_dz[ii];

      wx += vec_a[ii3+2] * dR_dx[ii];
      wy += vec_a[ii3+2] * dR_dy[ii];
      wz += vec_a[ii3+2] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    gwts = element->get_detJac(qua) * quad->get_qw(qua);

    get_f(coor_x, coor_y, coor_z, fx, fy, fz);
    
    // Residual here is the right-hand side forcing term
    for(A=0; A<nLocBas; ++A)
    {
      NA = R[A];
      NA_x = dR_dx[A]; NA_y = dR_dy[A]; NA_z = dR_dz[A];

      Residual[3*A] += gwts * ( 
          NA_x * (l2mu * ux + lambda * vy + lambda * wz)
        + NA_y * mu * (uy + vx)
        + NA_z * mu * (uz + wx) - NA * fx );
      
      Residual[3*A+1] += gwts * (
          NA_x * mu * (uy + vx)
          + NA_y * (lambda * ux + l2mu * vy + lambda * wz)
          + NA_z * mu * (vz + wy) - NA * fy );
     
      Residual[3*A+2] += gwts * (
          NA_x * mu * (uz + wx)
          + NA_y * mu * (vz + wy)
          + NA_z * (lambda * ux + lambda * vy + l2mu * wz )
          - NA * fz ); 
      
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


void PLocAssem_LinearElastic_3D_Static::Assem_Residual_EBC(
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

  for(qua = 0; qua < face_nqp; ++qua )
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
      Residual[3*A+0] -= gwts * R[A] * gx;
      Residual[3*A+1] -= gwts * R[A] * gy;
      Residual[3*A+2] -= gwts * R[A] * gz;
    }
  }
}


// EOF
