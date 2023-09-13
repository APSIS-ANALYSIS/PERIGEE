#include "FEAElement_Triangle6_membrane.hpp"

FEAElement_Triangle6_membrane::FEAElement_Triangle6_membrane( const int &in_nqua )
: nLocBas( 6 ), numQuapts( in_nqua )
{
  R     = new double [nLocBas * numQuapts];
  dR_dx = new double [nLocBas * numQuapts];
  dR_dy = new double [nLocBas * numQuapts];
  
  un = new Vector_3 [numQuapts];

  Jac    = new double [8 * numQuapts];
  detJac = new double [numQuapts];

  Q.resize(numQuapts);
}

FEAElement_Triangle6_membrane::~FEAElement_Triangle6_membrane()
{
  delete []     R;     R = nullptr;
  delete [] dR_dx; dR_dx = nullptr;
  delete [] dR_dy; dR_dy = nullptr;

  delete [] un; un = nullptr;

  delete []    Jac;    Jac = nullptr;
  delete [] detJac; detJac = nullptr;
}

void FEAElement_Triangle6_membrane::print_info() const
{
  SYS_T::commPrint("Triangle6_membrane: ");
  SYS_T::commPrint("6-node triangle element in the local lamina coordinate system. \n ");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d. \n", get_Type());
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

  double Rr [nLocBas], Rs [nLocBas];

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
 
    Rr[0] = 4.0 * qua_r + 4.0 * qua_s - 3.0;
    Rr[1] = 4.0 * qua_r - 1.0;
    Rr[2] = 0.0;
    Rr[3] = 4.0 - 8.0 * qua_r - 4.0 * qua_s;
    Rr[4] = 4.0 * qua_s;
    Rr[5] = -4.0 * qua_s;
    
    Rs[0] = 4.0 * qua_r + 4.0 * qua_s - 3.0;
    Rs[1] = 0.0;
    Rs[2] = 4.0 * qua_s - 1.0;
    Rs[3] = -4.0 * qua_r;
    Rs[4] = 4.0 * qua_r;
    Rs[5] = 4.0 - 4.0 * qua_r - 8.0 * qua_s;
    
    Vector_3 dx_dr( 0.0, 0.0, 0.0 );
    Vector_3 dx_ds( 0.0, 0.0, 0.0 );

    for( int ii=0; ii<nLocBas; ++ii )
    {
      const Vector_3 temp_dx_dr( ctrl_x[ii] * Rr[ii], ctrl_y[ii] * Rr[ii], ctrl_z[ii] * Rr[ii] );
      dx_dr += temp_dx_dr;

      const Vector_3 temp_dx_ds( ctrl_x[ii] * Rs[ii], ctrl_y[ii] * Rs[ii], ctrl_z[ii] * Rs[ii] );
      dx_ds += temp_dx_ds;
    }

    // vec(un) = vec(dx_dr) x vec(dx_ds)
    un[qua] = cross_product( dx_dr, dx_ds );
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
    Vector_3 e_b = cross_product( un[qua], e_a );
    e_b.normalize();

    // e_l1 = sqrt(2)/2 * (e_a - e_b)
    // e_l2 = sqrt(2)/2 * (e_a + e_b)
    const Vector_3 e_l1 = std::sqrt(2.0) * 0.5 * ( e_a - e_b );
    const Vector_3 e_l2 = std::sqrt(2.0) * 0.5 * ( e_a + e_b );

    // Q = transpose([ e_l1, e_l2, un ])
    Q[qua] = Matrix_3x3( e_l1, e_l2, un[qua] );
    Q[qua].transpose();

    // Rotated lamina coordinates
    double ctrl_xl [nLocBas], ctrl_yl [nLocBas];

    for(int ii = 0; ii < nLocBas; ++ii)
    {
      double temp_val;
      Q[qua].VecMult( ctrl_x[ii], ctrl_y[ii], ctrl_z[ii], ctrl_xl[ii], ctrl_yl[ii], temp_val );
    }

    // Rotated lamina 2D Jacobian & inverse Jacobian components
    double dxl_dr = 0.0, dxl_ds = 0.0, dyl_dr = 0.0, dyl_ds = 0.0;
    for( int ii = 0; ii < nLocBas; ++ii )
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

    for(int ii=0; ii<nLocBas; ++ii)
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
  const int offset = quaindex * nLocBas;
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
  const int offset = quaindex * nLocBas;
  return { R[offset], R[offset + 1], R[offset + 2],
   R[offset + 3], R[offset + 4], R[offset + 5] };
}

void FEAElement_Triangle6_membrane::get_gradR( const int &quaindex,
    double * const &basis_x, double * const &basis_y ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle6_membrane::get_gradR function error.\n" );
  const int offset = quaindex * nLocBas;
  for( int ii=0; ii<nLocBas; ++ii )
  {
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
  }
}

void FEAElement_Triangle6_membrane::get_R_gradR( const int &quaindex, 
    double * const &basis, double * const &basis_x, double * const &basis_y ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle6_membrane::get_R_gradR function error.\n" );
  const int offset = quaindex * nLocBas;
  for( int ii=0; ii<nLocBas; ++ii )
  {
    basis[ii]   = R[offset + ii];
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
  }
}

std::vector<double> FEAElement_Triangle6_membrane::get_dR_dx( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle6_membrane::get_dR_dx function error.\n" );
  const int offset = quaindex * nLocBas;
  return { dR_dx[offset], dR_dx[offset+1], dR_dx[offset+2],
    dR_dx[offset+3], dR_dx[offset+4], dR_dx[offset+5] }; 
}

std::vector<double> FEAElement_Triangle6_membrane::get_dR_dy( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle6_membrane::get_dR_dy function error.\n" );
  const int offset = quaindex * nLocBas;
  return { dR_dy[offset], dR_dy[offset+1], dR_dy[offset+2],
    dR_dy[offset+3], dR_dy[offset+4], dR_dy[offset+5] }; 
}

Vector_3 FEAElement_Triangle6_membrane::get_2d_normal_out( const int &qua,
    double &area ) const
{
  area = detJac[qua];
  return un[qua];
}

Matrix_3x3 FEAElement_Triangle6_membrane::get_rotationMatrix( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle6_membrane::get_rotationMatrix function error.\n" );
  return Q[quaindex];
}

Vector_3 FEAElement_Triangle6_membrane::get_normal_out( const int &qua,
    const double &sur_pt_x, const double &sur_pt_y, const double &sur_pt_z,
    const double &intpt_x, const double &intpt_y, const double &intpt_z,
    double &len ) const
{
  // Construct a vector starting from the interior point to the triangle
  const Vector_3 mm( sur_pt_x - intpt_x, sur_pt_y - intpt_y, sur_pt_z - intpt_z );

  // Do inner product with the normal vector
  const double mdotn = dot_product( mm, un[qua] );

  SYS_T::print_fatal_if( std::abs(mdotn) < 1.0e-10, "Warning: FEAElement_Triangle3_membrane::get_normal_out, the element might be ill-shaped.\n");

  len = detJac[qua];

  if(mdotn < 0) return (-1.0) * un[qua];
  else return un[qua];
}

// EOF
