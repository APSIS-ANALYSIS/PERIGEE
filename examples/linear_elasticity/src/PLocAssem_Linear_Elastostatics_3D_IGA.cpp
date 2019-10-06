#include "PLocAssem_Linear_Elastostatics_3D_IGA.hpp"

PLocAssem_Linear_Elastostatics_3D_IGA::PLocAssem_Linear_Elastostatics_3D_IGA(
    const double &in_mat_E, const double &in_mat_nu,
    const int &in_nlocbas, const int &in_nqp,
    const int &in_sdeg, const int &in_tdeg, const int &in_udeg,
    const int &in_nqpx, const int &in_nqpy, const int &in_nqpz )
: E( in_mat_E ), nu( in_mat_nu ), lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  nLocBas(in_nlocbas), dof_per_node(3), vec_size(nLocBas * dof_per_node),
  sdeg(in_sdeg), tdeg(in_tdeg), udeg(in_udeg),
  nqp(in_nqp), nqpx(in_nqpx), nqpy(in_nqpy), nqpz(in_nqpz)
{
  if( nLocBas != (sdeg+1) * (tdeg+1) * (udeg+1) ) SYS_T::print_fatal("Error: PLocAssem_Hyperelastic_3D_GenAlpha : nLocBas != (sdeg+1)*(tdeg+1)*(udeg+1). \n");

  if( nqp != (nqpx * nqpy * nqpz) ) SYS_T::print_fatal("Error: PLocAssem_Hyperelastic_3D_GenAlpha : nqp != nqpx * nqpy * nqpz. \n");

  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  Zero_Tangent_Residual();

  R = new double [nLocBas];
  dR_dx = new double [nLocBas];
  dR_dy = new double [nLocBas];
  dR_dz = new double [nLocBas];

  Sub_Tan = new double * [9];

  for(int ii=0; ii<9; ++ii) Sub_Tan[ii] = new double [nLocBas * nLocBas];

  print_info();
}


PLocAssem_Linear_Elastostatics_3D_IGA::~PLocAssem_Linear_Elastostatics_3D_IGA()
{
  delete [] Tangent; delete [] Residual; Tangent = NULL; Residual = NULL;
  delete [] R; delete [] dR_dx; delete [] dR_dy; delete [] dR_dz;
  for(int ii=0; ii<9; ++ii) delete [] Sub_Tan[ii];

  delete [] Sub_Tan;
}


void PLocAssem_Linear_Elastostatics_3D_IGA::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "  Three-dimensional Linear Elastostatic model IGA formulation: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Spatial: IGA \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Young's Modulus E  = %e \n", E);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Possion's ratio nu = %e \n", nu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Lame coeff lambda  = %e \n", lambda);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Shear modulus mu   = %e \n", mu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Bulk modulus kappa = %e \n", kappa);
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PLocAssem_Linear_Elastostatics_3D_IGA::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 0.0;
}

void PLocAssem_Linear_Elastostatics_3D_IGA::Zero_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
}


void PLocAssem_Linear_Elastostatics_3D_IGA::Assem_Estimate()
{
  for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 1.0;
}


void PLocAssem_Linear_Elastostatics_3D_IGA::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    FEAElement * const &element,
    const double &hx, const double &hy, const double &hz,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const BernsteinBasis_Array * const &bu,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &eleCtrlPts_w,
    const double * const &ext_x,
    const double * const &ext_y,
    const double * const &ext_z,
    const class AInt_Weight * const &weight )
{
  element->buildBasis( hx, hy, hz, bs, bt, bu, eleCtrlPts_x, eleCtrlPts_y,
      eleCtrlPts_z, eleCtrlPts_w, ext_x, ext_y, ext_z );

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

    gwts = element->get_detJac(qua) * weight->get_weight(qua);

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


void PLocAssem_Linear_Elastostatics_3D_IGA::Assem_Residual_TopFace(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    FEAElement * const &element,
    const double &hx, const double &hy, const double &hz,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &eleCtrlPts_w,
    const double * const &ext_x,
    const double * const &ext_y,
    const double * const &ext_z )
{
  IQuadPts * face_quad_z = new QuadPts_bc1();
  IQuadPts * face_quad_y = new QuadPts_Gauss(nqpy);
  IQuadPts * face_quad_x = new QuadPts_Gauss(nqpx);

  AInt_Weight * face_Int_w = new AInt_Weight(face_quad_x, face_quad_y,
      face_quad_z);

  BernsteinBasis_Array * ba_s = new BernsteinBasis_Array(sdeg, face_quad_x);
  BernsteinBasis_Array * ba_t = new BernsteinBasis_Array(tdeg, face_quad_y);
  BernsteinBasis_Array * ba_u = new BernsteinBasis_Array(udeg, face_quad_z);

  element->buildBasis( hx, hy, hz, ba_s, ba_t, ba_u, eleCtrlPts_x,
      eleCtrlPts_y, eleCtrlPts_z, eleCtrlPts_w, ext_x, ext_y, ext_z );

  const int face_nqp = face_Int_w->get_num();

  int ii, qua, A;
  double gwts, coor_x, coor_y, coor_z, gx, gy, gz, nx, ny, nz, surface_area;

  Zero_Residual();

  for(qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, R);
    element->get_3d_normal_top(qua, nx, ny, nz, surface_area);

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    for(ii=0; ii<nLocBas; ++ii)
    {
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    get_top_H(coor_x, coor_y, coor_z, 0.0, nx, ny, nz, gx, gy, gz);

    gwts = surface_area * hx * hy * face_Int_w->get_weight(qua);

    for(A=0; A<nLocBas; ++A)
    {
      Residual[3*A]   -= gwts * R[A] * gx;
      Residual[3*A+1] -= gwts * R[A] * gy;
      Residual[3*A+2] -= gwts * R[A] * gz;
    }
  }

  delete ba_s; delete ba_t; delete ba_u; delete face_Int_w;
  delete face_quad_x; delete face_quad_y; delete face_quad_z;
}


void PLocAssem_Linear_Elastostatics_3D_IGA::Assem_Residual_BotFace(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    FEAElement * const &element,
    const double &hx, const double &hy, const double &hz,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &eleCtrlPts_w,
    const double * const &ext_x,
    const double * const &ext_y,
    const double * const &ext_z )
{
  IQuadPts * face_quad_z = new QuadPts_bc0();
  IQuadPts * face_quad_y = new QuadPts_Gauss(nqpy);
  IQuadPts * face_quad_x = new QuadPts_Gauss(nqpx);

  AInt_Weight * face_Int_w = new AInt_Weight(face_quad_x, face_quad_y, 
      face_quad_z);

  BernsteinBasis_Array * ba_s = new BernsteinBasis_Array(sdeg, face_quad_x);
  BernsteinBasis_Array * ba_t = new BernsteinBasis_Array(tdeg, face_quad_y);
  BernsteinBasis_Array * ba_u = new BernsteinBasis_Array(udeg, face_quad_z);

  element->buildBasis( hx, hy, hz, ba_s, ba_t, ba_u, eleCtrlPts_x, 
      eleCtrlPts_y, eleCtrlPts_z, eleCtrlPts_w, ext_x, ext_y, ext_z );

  const int face_nqp = face_Int_w->get_num();

  int ii, qua, A;
  double gwts, coor_x, coor_y, coor_z, gx, gy, gz, nx, ny, nz, surface_area;

  Zero_Residual();

  for(qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, R);
    element->get_3d_normal_bottom(qua, nx, ny, nz, surface_area);

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    for(ii=0; ii<nLocBas; ++ii)
    {
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }
    get_bot_H(coor_x, coor_y, coor_z, 0.0, nx, ny, nz, gx, gy, gz);

    gwts = surface_area * hx * hy * face_Int_w->get_weight(qua);

    for(A=0; A<nLocBas; ++A)
    {
      Residual[3*A]   -= gwts * R[A] * gx;
      Residual[3*A+1] -= gwts * R[A] * gy;
      Residual[3*A+2] -= gwts * R[A] * gz; 
    }
  }

  delete ba_s; delete ba_t; delete ba_u; delete face_Int_w; 
  delete face_quad_x; delete face_quad_y; delete face_quad_z;
}


void PLocAssem_Linear_Elastostatics_3D_IGA::Assem_Residual_LefFace(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    FEAElement * const &element,
    const double &hx, const double &hy, const double &hz,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &eleCtrlPts_w,
    const double * const &ext_x,
    const double * const &ext_y,
    const double * const &ext_z )
{
  IQuadPts * face_quad_z = new QuadPts_Gauss(nqpz);
  IQuadPts * face_quad_y = new QuadPts_bc0();
  IQuadPts * face_quad_x = new QuadPts_Gauss(nqpx);

  AInt_Weight * face_Int_w = new AInt_Weight(face_quad_x, face_quad_y, 
      face_quad_z);

  BernsteinBasis_Array * ba_s = new BernsteinBasis_Array(sdeg, face_quad_x);
  BernsteinBasis_Array * ba_t = new BernsteinBasis_Array(tdeg, face_quad_y);
  BernsteinBasis_Array * ba_u = new BernsteinBasis_Array(udeg, face_quad_z);

  element->buildBasis( hx, hy, hz, ba_s, ba_t, ba_u, eleCtrlPts_x, 
      eleCtrlPts_y, eleCtrlPts_z, eleCtrlPts_w, ext_x, ext_y, ext_z );

  const int face_nqp = face_Int_w->get_num();

  int ii, qua, A;
  double gwts, coor_x, coor_y, coor_z, gx, gy, gz, nx, ny, nz, surface_area;

  Zero_Residual();

  for(qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, R);
    element->get_3d_normal_left(qua, nx, ny, nz, surface_area);

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    for(ii=0; ii<nLocBas; ++ii)
    {
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    get_lef_H(coor_x, coor_y, coor_z, 0.0, nx, ny, nz, gx, gy, gz);

    gwts = surface_area * hx * hz * face_Int_w->get_weight(qua);

    for(A=0; A<nLocBas; ++A)
    {
      Residual[3*A]   -= gwts * R[A] * gx;
      Residual[3*A+1] -= gwts * R[A] * gy;
      Residual[3*A+2] -= gwts * R[A] * gz; 
    }
  }

  delete ba_s; delete ba_t; delete ba_u; delete face_Int_w; 
  delete face_quad_x; delete face_quad_y; delete face_quad_z;
}


void PLocAssem_Linear_Elastostatics_3D_IGA::Assem_Residual_RigFace(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    FEAElement * const &element,
    const double &hx, const double &hy, const double &hz,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &eleCtrlPts_w,
    const double * const &ext_x,
    const double * const &ext_y,
    const double * const &ext_z )
{
  IQuadPts * face_quad_z = new QuadPts_Gauss(nqpz);
  IQuadPts * face_quad_y = new QuadPts_bc1();
  IQuadPts * face_quad_x = new QuadPts_Gauss(nqpx);

  AInt_Weight * face_Int_w = new AInt_Weight(face_quad_x, face_quad_y, 
      face_quad_z);

  BernsteinBasis_Array * ba_s = new BernsteinBasis_Array(sdeg, face_quad_x);
  BernsteinBasis_Array * ba_t = new BernsteinBasis_Array(tdeg, face_quad_y);
  BernsteinBasis_Array * ba_u = new BernsteinBasis_Array(udeg, face_quad_z);

  element->buildBasis( hx, hy, hz, ba_s, ba_t, ba_u, eleCtrlPts_x, 
      eleCtrlPts_y, eleCtrlPts_z, eleCtrlPts_w, ext_x, ext_y, ext_z );

  const int face_nqp = face_Int_w->get_num();

  int ii, qua, A;
  double gwts, coor_x, coor_y, coor_z, gx, gy, gz, nx, ny, nz, surface_area;

  Zero_Residual();

  for(qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, R);
    element->get_3d_normal_right(qua, nx, ny, nz, surface_area);

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    for(ii=0; ii<nLocBas; ++ii)
    {
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }
    get_rig_H(coor_x, coor_y, coor_z, 0.0, nx, ny, nz, gx, gy, gz);

    gwts = surface_area * hx * hz * face_Int_w->get_weight(qua);

    for(A=0; A<nLocBas; ++A)
    {
      Residual[3*A]   -= gwts * R[A] * gx;
      Residual[3*A+1] -= gwts * R[A] * gy;
      Residual[3*A+2] -= gwts * R[A] * gz; 
    }
  }

  delete ba_s; delete ba_t; delete ba_u; delete face_Int_w; 
  delete face_quad_x; delete face_quad_y; delete face_quad_z;
}


void PLocAssem_Linear_Elastostatics_3D_IGA::Assem_Residual_FroFace(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    FEAElement * const &element,
    const double &hx, const double &hy, const double &hz,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &eleCtrlPts_w,
    const double * const &ext_x,
    const double * const &ext_y,
    const double * const &ext_z )
{
  IQuadPts * face_quad_z = new QuadPts_Gauss(nqpz);
  IQuadPts * face_quad_y = new QuadPts_Gauss(nqpy);
  IQuadPts * face_quad_x = new QuadPts_bc1();

  AInt_Weight * face_Int_w = new AInt_Weight(face_quad_x, face_quad_y, 
      face_quad_z);

  BernsteinBasis_Array * ba_s = new BernsteinBasis_Array(sdeg, face_quad_x);
  BernsteinBasis_Array * ba_t = new BernsteinBasis_Array(tdeg, face_quad_y);
  BernsteinBasis_Array * ba_u = new BernsteinBasis_Array(udeg, face_quad_z);

  element->buildBasis( hx, hy, hz, ba_s, ba_t, ba_u, eleCtrlPts_x, 
      eleCtrlPts_y, eleCtrlPts_z, eleCtrlPts_w, ext_x, ext_y, ext_z );

  const int face_nqp = face_Int_w->get_num();

  int ii, qua, A;
  double gwts, coor_x, coor_y, coor_z, gx, gy, gz, nx, ny, nz, surface_area;

  Zero_Residual();

  for(qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, R);
    element->get_3d_normal_front(qua, nx, ny, nz, surface_area);

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    for(ii=0; ii<nLocBas; ++ii)
    {
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }
    get_fro_H(coor_x, coor_y, coor_z, 0.0, nx, ny, nz, gx, gy, gz);

    gwts = surface_area * hz * hy * face_Int_w->get_weight(qua);

    for(A=0; A<nLocBas; ++A)
    {
      Residual[3*A]   -= gwts * R[A] * gx;
      Residual[3*A+1] -= gwts * R[A] * gy;
      Residual[3*A+2] -= gwts * R[A] * gz; 
    }
  }

  delete ba_s; delete ba_t; delete ba_u; delete face_Int_w; 
  delete face_quad_x; delete face_quad_y; delete face_quad_z;
}


void PLocAssem_Linear_Elastostatics_3D_IGA::Assem_Residual_BacFace(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    FEAElement * const &element,
    const double &hx, const double &hy, const double &hz,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &eleCtrlPts_w,
    const double * const &ext_x,
    const double * const &ext_y,
    const double * const &ext_z )
{
  IQuadPts * face_quad_z = new QuadPts_Gauss(nqpz);
  IQuadPts * face_quad_y = new QuadPts_Gauss(nqpy);
  IQuadPts * face_quad_x = new QuadPts_bc0();

  AInt_Weight * face_Int_w = new AInt_Weight(face_quad_x, face_quad_y, 
      face_quad_z);

  BernsteinBasis_Array * ba_s = new BernsteinBasis_Array(sdeg, face_quad_x);
  BernsteinBasis_Array * ba_t = new BernsteinBasis_Array(tdeg, face_quad_y);
  BernsteinBasis_Array * ba_u = new BernsteinBasis_Array(udeg, face_quad_z);

  element->buildBasis( hx, hy, hz, ba_s, ba_t, ba_u, eleCtrlPts_x, 
      eleCtrlPts_y, eleCtrlPts_z, eleCtrlPts_w, ext_x, ext_y, ext_z );

  const int face_nqp = face_Int_w->get_num();

  int ii, qua, A;
  double gwts, coor_x, coor_y, coor_z, gx, gy, gz, nx, ny, nz, surface_area;

  Zero_Residual();

  for(qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, R);

    element->get_3d_normal_back(qua, nx, ny, nz, surface_area);

    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;
    for(ii=0; ii<nLocBas; ++ii)
    {
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    get_bac_H(coor_x, coor_y, coor_z, 0.0, nx, ny, nz, gx, gy, gz);

    gwts = surface_area * hz * hy * face_Int_w->get_weight(qua);

    for(A=0; A<nLocBas; ++A)
    {
      Residual[3*A]   -= gwts * R[A] * gx;
      Residual[3*A+1] -= gwts * R[A] * gy;
      Residual[3*A+2] -= gwts * R[A] * gz; 
    }
  }

  delete ba_s; delete ba_t; delete ba_u; delete face_Int_w; 
  delete face_quad_x; delete face_quad_y; delete face_quad_z;
}

// EOF
