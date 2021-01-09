#include "FEAElement_Triangle3_membrane.hpp"

FEAElement_Triangle3_membrane::FEAElement_Triangle3_membrane( 
    const int &in_nqua )
: nLocBas( 3 ), numQuapts( in_nqua )
{
  ctrl_xl   = new double [ nLocBas ];
  ctrl_yl   = new double [ nLocBas ];
  ctrl_zl   = new double [ nLocBas ];

  R         = new double [ nLocBas * numQuapts ];
  dR_dx     = new double [ nLocBas ];
  dR_dy     = new double [ nLocBas ];
}


FEAElement_Triangle3_membrane::~FEAElement_Triangle3_membrane()
{
  clearBasisCache();
}


void FEAElement_Triangle3_membrane::clearBasisCache()
{
  delete [] ctrl_xl;   ctrl_xl   = nullptr;
  delete [] ctrl_yl;   ctrl_yl   = nullptr;
  delete [] ctrl_zl;   ctrl_zl   = nullptr;

  delete [] R;         R = nullptr;
  delete [] dR_dx; dR_dx = nullptr;
  delete [] dR_dy; dR_dy = nullptr;
}


void FEAElement_Triangle3_membrane::resize_container()
{
  clearBasisCache();
  R = new double [ nLocBas * numQuapts ];
}


void FEAElement_Triangle3_membrane::reset_numQuapts( 
    const int &new_num_qua )
{
  numQuapts = new_num_qua;
  resize_container();
}


void FEAElement_Triangle3_membrane::print() const
{
  SYS_T::commPrint("Triangle3_membrane: ");
  SYS_T::commPrint("3-node triangle element in the local lamina coordinate system. \n ");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d. \n", get_Type());
  SYS_T::commPrint("Note: This element is designed for the coupled momentum method. \n ");
}


double FEAElement_Triangle3_membrane::get_memory_usage() const
{
  double double_size = nLocBas * numQuapts + 10.0;
  double int_size = 2;
  return double_size * 8.0 + int_size * 4.0;
}


void FEAElement_Triangle3_membrane::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z )
{
  for( int qua = 0; qua < numQuapts; ++qua )
  {
    const double qua_r = quad -> get_qp(qua, 0);
    const double qua_s = quad -> get_qp(qua, 1);
    R[qua*3 + 0] = 1.0 - qua_r - qua_s;
    R[qua*3 + 1] = qua_r;
    R[qua*3 + 2] = qua_s;
  }

  //const double dN0_dr = -1.0; const double dN0_ds = -1.0;
  //const double dN1_dr = 1.0;  const double dN1_ds = 0.0;
  //const double dN2_dr = 0.0;  const double dN2_ds = 1.0;
  
  const double dx_dr = ctrl_x[0] * (-1.0) + ctrl_x[1];
  const double dy_dr = ctrl_y[0] * (-1.0) + ctrl_y[1];
  const double dz_dr = ctrl_z[0] * (-1.0) + ctrl_z[1];
  const double dx_ds = ctrl_x[0] * (-1.0) + ctrl_x[2];
  const double dy_ds = ctrl_y[0] * (-1.0) + ctrl_y[2];
  const double dz_ds = ctrl_z[0] * (-1.0) + ctrl_z[2];

  // vec(un) = vec(dx_dr) x vec(dx_ds)
  MATH_T::cross3d(dx_dr, dy_dr, dz_dr, dx_ds, dy_ds, dz_ds,
      unx, uny, unz);

  MATH_T::normalize3d( unx, uny, unz );

  // ======= Global-to-local rotation matrix =======
  const double inv_len_er = 1.0 / MATH_T::norm2( dx_dr, dy_dr, dz_dr );
  e_r[0] = dx_dr * inv_len_er;
  e_r[1] = dy_dr * inv_len_er; 
  e_r[2] = dz_dr * inv_len_er; 

  const double inv_len_es = 1.0 / MATH_T::norm2( dx_ds, dy_ds, dz_ds );
  e_s[0] = dx_ds * inv_len_es; 
  e_s[1] = dy_ds * inv_len_es; 
  e_s[2] = dz_ds * inv_len_es; 

  // e_a = 0.5*(e_r + e_s) / || 0.5*(e_r + e_s) ||
  double e_a[3];
  e_a[0] = 0.5 * ( e_r[0] + e_s[0] );
  e_a[1] = 0.5 * ( e_r[1] + e_s[1] );
  e_a[2] = 0.5 * ( e_r[2] + e_s[2] );

  MATH_T::normalize3d( e_a[0], e_a[1], e_a[2] );

  // e_b = vec(un) x e_a / || vec(un) x e_a ||
  double e_b[3];
  MATH_T::cross3d(unx, uny, unz, e_a[0], e_a[1], e_a[2],
      e_b[0], e_b[1], e_b[2]);
  MATH_T::normalize3d( e_b[0], e_b[1], e_b[2] );

  // e_l1 = sqrt(2)/2 * (e_a - e_b)
  // e_l2 = sqrt(2)/2 * (e_a + e_b)
  e_l1[0] = std::sqrt(2.0) * 0.5 * ( e_a[0] - e_b[0] );
  e_l2[0] = std::sqrt(2.0) * 0.5 * ( e_a[0] + e_b[0] );

  e_l1[1] = std::sqrt(2.0) * 0.5 * ( e_a[1] - e_b[1] );
  e_l2[1] = std::sqrt(2.0) * 0.5 * ( e_a[1] + e_b[1] );

  e_l1[2] = std::sqrt(2.0) * 0.5 * ( e_a[2] - e_b[2] );
  e_l2[2] = std::sqrt(2.0) * 0.5 * ( e_a[2] + e_b[2] );

  // Q = transpose([ e_l1, e_l2, un ])
  Q = Matrix_3x3(e_l1[0], e_l1[1], e_l1[2],
                 e_l2[0], e_l2[1], e_l2[2],
                     unx,     uny,     unz );

  // Rotated local coordinates
  for(int ii = 0; ii < nLocBas; ++ii)
  {
    double ctrl_xyzl[3] = {0.0, 0.0, 0.0};
    Q.VecMult( ctrl_x[ii], ctrl_y[ii], ctrl_z[ii], &ctrl_xyzl[0]);
    ctrl_xl[ii] = ctrl_xyzl[0];
    ctrl_yl[ii] = ctrl_xyzl[1];
    ctrl_zl[ii] = ctrl_xyzl[2];
  }

  // Rotated local 2D Jacobian components
  Jac[0] = ctrl_xl[0] * (-1.0) + ctrl_xl[1]; // dxl_dr 
  Jac[1] = ctrl_xl[0] * (-1.0) + ctrl_xl[2]; // dxl_ds

  Jac[2] = ctrl_yl[0] * (-1.0) + ctrl_yl[1]; // dyl_dr
  Jac[3] = ctrl_yl[0] * (-1.0) + ctrl_yl[2]; // dyl_ds

  detJac = Jac[0] * Jac[3] - Jac[1] * Jac[2];

  double inv_detJac = 1.0 / detJac;

  Jac[4] = Jac[3] * inv_detJac;        // dr_dx
  Jac[5] = -1.0 * Jac[1] * inv_detJac; // dr_dy
  Jac[6] = -1.0 * Jac[2] * inv_detJac; // ds_dx
  Jac[7] = Jac[0] * inv_detJac;        // ds_dy

  dR_dx[0] = (-1.0) * Jac[4] - Jac[6];
  dR_dx[1] = Jac[4];
  dR_dx[2] = Jac[6];

  dR_dy[0] = (-1.0) * Jac[5] - Jac[7];
  dR_dy[1] = Jac[5];
  dR_dy[2] = Jac[7];
}


void FEAElement_Triangle3_membrane::get_R( const int &quaindex, 
    double * const &basis ) const
{
  assert(quaindex>=0 && quaindex < numQuapts);
  const int offset = quaindex * nLocBas;
  basis[0] = R[offset];
  basis[1] = R[offset+1];
  basis[2] = R[offset+2];
}


void FEAElement_Triangle3_membrane::get_gradR( const int &quaindex,
    double * const &basis_x, double * const &basis_y, double * const &basis_z ) const
{
  // TODO
}


void FEAElement_Triangle3_membrane::get_R_gradR( const int &quaindex,
    double * const &basis, double * const &basis_x, double * const &basis_y,
    double * const &basis_z ) const
{
  // TODO
}


void FEAElement_Triangle3_membrane::get_2d_normal_out( const int &quaindex,
    double &nx, double &ny, double &nz, double &area ) const
{
  assert(quaindex>=0 && quaindex < numQuapts);
  nx = unx;
  ny = uny;
  nz = unz;
  area = detJac;
}


void FEAElement_Triangle3_membrane::get_rotationMatrix( const int &quaindex,
    Matrix_3x3 &rot_mat) const
{
  rot_mat = Q;
}


void FEAElement_Triangle3_membrane::get_normal_out( const int &quaindex,
    const double &sur_pt_x, const double &sur_pt_y, const double &sur_pt_z,
    const double &intpt_x, const double &intpt_y, const double &intpt_z,
    double &nx, double &ny, double &nz, double &area ) const
{
  // Construct a vector from the interior point to the triangle first node
  const double mx = sur_pt_x - intpt_x;
  const double my = sur_pt_y - intpt_y;
  const double mz = sur_pt_z - intpt_z;

  // Dot product of the defined vector with the calculated normal vector
  const double mdotn = mx * unx + my * uny + mz * unz;

  SYS_T::print_fatal_if( std::abs(mdotn) < 1.0e-10, "Warning: FEAElement_Triangle3_membrane::get_normal_out, the element might be ill-shaped.\n");

  // If dot product is negative, adjust the normal vector
  if(mdotn < 0)
  {
    nx = (-1.0) * unx;
    ny = (-1.0) * uny;
    nz = (-1.0) * unz;
  }
  else
  {
    nx = unx;
    ny = uny;
    nz = unz;
  }

  area = detJac;
}

// EOF
