#include "FEAElement_Triangle6_membrane.hpp"

FEAElement_Triangle6_membrane::FEAElement_Triangle6_membrane( const int &in_nqua )
: numQuapts( in_nqua )
{
  R     = new double [6 * numQuapts];
  dR_dx = new double [6 * numQuapts];
  dR_dy = new double [6 * numQuapts];
  
  Jac    = new double [8 * numQuapts];
  detJac = new double [numQuapts];

  un.resize(numQuapts);
  Q.resize(numQuapts);
}

FEAElement_Triangle6_membrane::~FEAElement_Triangle6_membrane()
{
  delete []     R;     R = nullptr;
  delete [] dR_dx; dR_dx = nullptr;
  delete [] dR_dy; dR_dy = nullptr;

  delete []    Jac;    Jac = nullptr;
  delete [] detJac; detJac = nullptr;
}

void FEAElement_Triangle6_membrane::print_info() const
{
  SYS_T::commPrint("Triangle6_membrane: ");
  SYS_T::commPrint("6-node triangle element in the local lamina coordinate system. \n ");
  SYS_T::commPrint("elemType: %d. \n", get_Type());
  SYS_T::commPrint("Note: This element is designed for the coupled momentum method. \n ");
}

double FEAElement_Triangle6_membrane::get_memory_usage() const
{
  const double dsize = 16 * numQuapts;
  return dsize * 8.0 + 4.0;
}

void FEAElement_Triangle6_membrane::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z )
{
  ASSERT(quad->get_dim() == 3, "FEAElement_Triangle6_membrane::buildBasis function error.\n" );

  for( int qua = 0; qua < numQuapts; ++qua )
  {
    const double qua_r = quad -> get_qp( qua, 0 );
    const double qua_s = quad -> get_qp( qua, 1 );
    const double qua_t = quad -> get_qp( qua, 2 );

    const int offset = 6 * qua;

    R[offset + 0] = qua_t * (2.0 * qua_t - 1.0);
    R[offset + 1] = qua_r * (2.0 * qua_r - 1.0);
    R[offset + 2] = qua_s * (2.0 * qua_s - 1.0);
    R[offset + 3] = 4.0 * qua_r * qua_t;
    R[offset + 4] = 4.0 * qua_r * qua_s;
    R[offset + 5] = 4.0 * qua_s * qua_t;
 
    const double Rr[6] = { 4.0 * qua_r + 4.0 * qua_s - 3.0,
      4.0 * qua_r - 1.0,
      0.0,
      4.0 - 8.0 * qua_r - 4.0 * qua_s,
      4.0 * qua_s,
      -4.0 * qua_s };
    
    const double Rs[6] = { 4.0 * qua_r + 4.0 * qua_s - 3.0,
      0.0,
      4.0 * qua_s - 1.0,
      -4.0 * qua_r,
      4.0 * qua_r,
      4.0 - 4.0 * qua_r - 8.0 * qua_s };
    
    Vector_3 dx_dr( 0.0, 0.0, 0.0 );
    Vector_3 dx_ds( 0.0, 0.0, 0.0 );

    for( int ii=0; ii<6; ++ii )
    {
      dx_dr += Vector_3( ctrl_x[ii] * Rr[ii], ctrl_y[ii] * Rr[ii], ctrl_z[ii] * Rr[ii] );
      dx_ds += Vector_3( ctrl_x[ii] * Rs[ii], ctrl_y[ii] * Rs[ii], ctrl_z[ii] * Rs[ii] );
    }

    // vec(un) = vec(dx_dr) x vec(dx_ds)
    un[qua] = VEC3_T::cross_product( dx_dr, dx_ds );
    un[qua].normalize();

    // ======= Global-to-local rotation matrix =======
    const double inv_len_er = 1.0 / dx_dr.norm2();
    const Vector_3 e_r = inv_len_er * dx_dr;

    const double inv_len_es = 1.0 / dx_ds.norm2();
    const Vector_3 e_s = inv_len_es * dx_ds;

    // e_a = 0.5*(e_r + e_s) / || 0.5*(e_r + e_s) ||
    Vector_3 e_a = 0.5 * ( e_r + e_s );
    e_a.normalize();

    // e_b = vec(un) x e_a / || vec(un) x e_a ||
    Vector_3 e_b = VEC3_T::cross_product( un[qua], e_a );
    e_b.normalize();

    // e_l1 = sqrt(2)/2 * (e_a - e_b)
    // e_l2 = sqrt(2)/2 * (e_a + e_b)
    const Vector_3 e_l1 = std::sqrt(2.0) * 0.5 * ( e_a - e_b );
    const Vector_3 e_l2 = std::sqrt(2.0) * 0.5 * ( e_a + e_b );

    // Q = transpose([ e_l1, e_l2, un ])
    Q[qua] = Tensor2_3D( e_l1, e_l2, un[qua] );
    Q[qua].transpose();

    // Rotated lamina coordinates
    double ctrl_xl [6], ctrl_yl [6];

    for(int ii = 0; ii < 6; ++ii)
    {
      double temp_val;
      Q[qua].VecMult( ctrl_x[ii], ctrl_y[ii], ctrl_z[ii], ctrl_xl[ii], ctrl_yl[ii], temp_val );
    }

    // Rotated lamina 2D Jacobian & inverse Jacobian components
    double dxl_dr = 0.0, dxl_ds = 0.0, dyl_dr = 0.0, dyl_ds = 0.0;
    for( int ii = 0; ii < 6; ++ii )
    {
      dxl_dr += ctrl_xl[ii] * Rr[ii];
      dxl_ds += ctrl_xl[ii] * Rs[ii];

      dyl_dr += ctrl_yl[ii] * Rr[ii];
      dyl_ds += ctrl_yl[ii] * Rs[ii];
    }

    Jac[4*qua+0] = dxl_dr;
    Jac[4*qua+1] = dxl_ds;
    Jac[4*qua+2] = dyl_dr;
    Jac[4*qua+3] = dyl_ds;

    detJac[qua] = dxl_dr * dyl_ds - dxl_ds * dyl_dr;

    const double inv_detJac = 1.0 / detJac[qua];

    const double dr_dxl = Jac[4*qua+3] * inv_detJac;
    const double dr_dyl = (-1.0) * Jac[4*qua+1] * inv_detJac;
    const double ds_dxl = (-1.0) * Jac[4*qua+2] * inv_detJac;
    const double ds_dyl = Jac[4*qua+0] * inv_detJac;

    Jac[4*numQuapts + 4*qua + 0] = dr_dxl;
    Jac[4*numQuapts + 4*qua + 1] = dr_dyl;
    Jac[4*numQuapts + 4*qua + 2] = ds_dxl;
    Jac[4*numQuapts + 4*qua + 3] = ds_dyl;

    for(int ii=0; ii<6; ++ii)
    {
      dR_dx[offset+ii] = Rr[ii] * dr_dxl + Rs[ii] * ds_dxl;
      dR_dy[offset+ii] = Rr[ii] * dr_dyl + Rs[ii] * ds_dyl;
    }

  } // end qua loop
}

void FEAElement_Triangle6_membrane::get_R( 
    const int &quaindex, double * const &basis ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle6_membrane::get_R function error.\n" );
  const int offset = quaindex * 6;
  basis[0] = R[offset];
  basis[1] = R[offset + 1];
  basis[2] = R[offset + 2];
  basis[3] = R[offset + 3];
  basis[4] = R[offset + 4];
  basis[5] = R[offset + 5];
}

std::vector<double> FEAElement_Triangle6_membrane::get_R( 
    const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle6_membrane::get_R function error.\n" );
  const int offset = quaindex * 6;
  return { R[offset], R[offset + 1], R[offset + 2],
   R[offset + 3], R[offset + 4], R[offset + 5] };
}

void FEAElement_Triangle6_membrane::get_gradR( const int &quaindex,
    double * const &basis_x, double * const &basis_y ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle6_membrane::get_gradR function error.\n" );
  const int offset = quaindex * 6;
  for( int ii=0; ii<6; ++ii )
  {
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
  }
}

void FEAElement_Triangle6_membrane::get_R_gradR( const int &quaindex, 
    double * const &basis, double * const &basis_x, double * const &basis_y ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle6_membrane::get_R_gradR function error.\n" );
  const int offset = quaindex * 6;
  for( int ii=0; ii<6; ++ii )
  {
    basis[ii]   = R[offset + ii];
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
  }
}

std::vector<double> FEAElement_Triangle6_membrane::get_dR_dx( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle6_membrane::get_dR_dx function error.\n" );
  const int offset = quaindex * 6;
  return { dR_dx[offset], dR_dx[offset+1], dR_dx[offset+2],
    dR_dx[offset+3], dR_dx[offset+4], dR_dx[offset+5] }; 
}

std::vector<double> FEAElement_Triangle6_membrane::get_dR_dy( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle6_membrane::get_dR_dy function error.\n" );
  const int offset = quaindex * 6;
  return { dR_dy[offset], dR_dy[offset+1], dR_dy[offset+2],
    dR_dy[offset+3], dR_dy[offset+4], dR_dy[offset+5] }; 
}

Vector_3 FEAElement_Triangle6_membrane::get_2d_normal_out( const int &qua,
    double &area ) const
{
  area = detJac[qua];
  return un[qua];
}

Tensor2_3D FEAElement_Triangle6_membrane::get_rotationMatrix( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle6_membrane::get_rotationMatrix function error.\n" );
  return Q[quaindex];
}

Vector_3 FEAElement_Triangle6_membrane::get_normal_out( const int &qua,
    const Vector_3 &sur_pt, const Vector_3 &int_pt, double &len ) const
{
  // Construct a vector starting from the interior point to the triangle
  const Vector_3 mm = sur_pt - int_pt;

  // Do inner product with the normal vector
  const double mdotn = VEC3_T::dot_product( mm, un[qua] );

  SYS_T::print_fatal_if( std::abs(mdotn) < 1.0e-10, "Warning: FEAElement_Triangle3_membrane::get_normal_out, the element might be ill-shaped.\n");

  len = detJac[qua];

  if(mdotn < 0) return (-1.0) * un[qua];
  else return un[qua];
}

// EOF
