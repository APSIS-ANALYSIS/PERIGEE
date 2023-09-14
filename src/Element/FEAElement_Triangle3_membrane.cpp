#include "FEAElement_Triangle3_membrane.hpp"

FEAElement_Triangle3_membrane::FEAElement_Triangle3_membrane( 
    const int &in_nqua ) : nLocBas( 3 ), numQuapts( in_nqua )
{
  R = new double [nLocBas * numQuapts];
}

FEAElement_Triangle3_membrane::~FEAElement_Triangle3_membrane()
{
  delete [] R; R = nullptr;
}

void FEAElement_Triangle3_membrane::print_info() const
{
  SYS_T::commPrint("Triangle3_membrane: ");
  SYS_T::commPrint("3-node triangle element in the local lamina coordinate system. \n ");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d. \n", get_Type());
  SYS_T::commPrint("Note: This element is designed for the coupled momentum method. \n ");
}

double FEAElement_Triangle3_membrane::get_memory_usage() const
{
  double double_size = nLocBas * numQuapts + 5*nLocBas + 33;
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

  const Vector_3 dx_dr(ctrl_x[0] * (-1.0) + ctrl_x[1],
    ctrl_y[0] * (-1.0) + ctrl_y[1],
    ctrl_z[0] * (-1.0) + ctrl_z[1]);

  const Vector_3 dx_ds(ctrl_x[0] * (-1.0) + ctrl_x[2],
    ctrl_y[0] * (-1.0) + ctrl_y[2],
    ctrl_z[0] * (-1.0) + ctrl_z[2]);

  // vec(un) = vec(dx_dr) x vec(dx_ds)
  un = cross_product( dx_dr, dx_ds );
  un.normalize();

  // ======= Global-to-local rotation matrix =======
  const double inv_len_er = 1.0 / dx_dr.norm2();
  const Vector_3 e_r = inv_len_er * dx_dr;

  const double inv_len_es = 1.0 / dx_ds.norm2();
  const Vector_3 e_s = inv_len_es * dx_ds;

  // e_a = 0.5*(e_r + e_s) / || 0.5*(e_r + e_s) ||
  Vector_3 e_a = 0.5 * (e_r + e_s);
  e_a.normalize();

  // e_b = vec(un) x e_a / || vec(un) x e_a ||
  Vector_3 e_b = cross_product( un, e_a );
  e_b.normalize();

  // e_l1 = sqrt(2)/2 * (e_a - e_b)
  // e_l2 = sqrt(2)/2 * (e_a + e_b)
  const Vector_3 e_l1 = std::sqrt(2.0) * 0.5 * ( e_a - e_b );
  
  const Vector_3 e_l2 = std::sqrt(2.0) * 0.5 * ( e_a + e_b );

  // Q = transpose([ e_l1, e_l2, un ])
  Q = Matrix_3x3( e_l1, e_l2, un );
  Q.transpose();
  
  // Rotated lamina coordinates
  double ctrl_xl [nLocBas], ctrl_yl [nLocBas];

  for(int ii = 0; ii < nLocBas; ++ii)
  {
    double temp_val;
    Q.VecMult( ctrl_x[ii], ctrl_y[ii], ctrl_z[ii], ctrl_xl[ii], ctrl_yl[ii], temp_val );
  }

  // Rotated lamina 2D Jacobian & inverse Jacobian components
  Jac[0] = ctrl_xl[0] * (-1.0) + ctrl_xl[1]; // dxl_dr 
  Jac[1] = ctrl_xl[0] * (-1.0) + ctrl_xl[2]; // dxl_ds

  Jac[2] = ctrl_yl[0] * (-1.0) + ctrl_yl[1]; // dyl_dr
  Jac[3] = ctrl_yl[0] * (-1.0) + ctrl_yl[2]; // dyl_ds

  detJac = Jac[0] * Jac[3] - Jac[1] * Jac[2];

  double inv_detJac = 1.0 / detJac;

  Jac[4] = Jac[3] * inv_detJac;              // dr_dxl
  Jac[5] = -1.0 * Jac[1] * inv_detJac;       // dr_dyl
  Jac[6] = -1.0 * Jac[2] * inv_detJac;       // ds_dxl
  Jac[7] = Jac[0] * inv_detJac;              // ds_dyl

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
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle3_membrane::get_R function error.\n" );
  const int offset = quaindex * nLocBas;
  basis[0] = R[offset];
  basis[1] = R[offset + 1];
  basis[2] = R[offset + 2];
}

std::vector<double> FEAElement_Triangle3_membrane::get_R( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle3_membrane::get_R function error.\n" );
  const int offset = quaindex * nLocBas;
  return { R[offset], R[offset + 1], R[offset + 2] };
}

void FEAElement_Triangle3_membrane::get_gradR( const int &quaindex,
    double * const &basis_x, double * const &basis_y ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle3_membrane::get_gradR function error.\n" );
  for( int ii=0; ii<nLocBas; ++ii )
  {
    basis_x[ii] = dR_dx[ii];
    basis_y[ii] = dR_dy[ii];
  }
}

void FEAElement_Triangle3_membrane::get_R_gradR( const int &quaindex,
    double * const &basis, double * const &basis_x, double * const &basis_y ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle3_membrane::get_R_gradR function error.\n" );
  const int offset = quaindex * nLocBas;
  for( int ii=0; ii<nLocBas; ++ii )
  {
    basis[ii]   = R[offset + ii];
    basis_x[ii] = dR_dx[ii];
    basis_y[ii] = dR_dy[ii];
  }
}

Vector_3 FEAElement_Triangle3_membrane::get_2d_normal_out( const int &quaindex, 
    double &area ) const
{
  ASSERT(quaindex>=0 && quaindex < numQuapts, "FEAElement_Triangle3_membrane::get_2d_normal_out function error.\n" );
  area = detJac;
  return un;
}

Matrix_3x3 FEAElement_Triangle3_membrane::get_rotationMatrix( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle3_membrane::get_rotationMatrix function error.\n" );
  return Q;
}

Vector_3 FEAElement_Triangle3_membrane::get_normal_out( const int &quaindex,
    const Vector_3 &sur_pt, const Vector_3 &int_pt, double &area ) const
{
  // Construct a vector from the interior point to the triangle first node
  const Vector_3 mm = sur_pt - int_pt;

  // Dot product of the defined vector with the calculated normal vector
  const double mdotn = dot_product( mm, un );

  SYS_T::print_fatal_if( std::abs(mdotn) < 1.0e-10, "Warning: FEAElement_Triangle3_membrane::get_normal_out, the element might be ill-shaped.\n");

  area = detJac;

  // If dot product is negative, adjust the normal vector
  if(mdotn < 0) return (-1.0) * un;
  else return un;
}

// EOF
