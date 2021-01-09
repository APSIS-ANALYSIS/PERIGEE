#include "FEAElement_Triangle6_membrane.hpp"

FEAElement_Triangle6_membrane::FEAElement_Triangle6_membrane( const int &in_nqua )
: nLocBas( 6 ), numQuapts( in_nqua )
{
  R = new double [nLocBas * numQuapts];
  
  unx = new double [numQuapts];
  uny = new double [numQuapts];
  unz = new double [numQuapts];

  detJac = new double [numQuapts];

  e_r.resize(  numQuapts );
  e_s.resize(  numQuapts );
  e_l1.resize( numQuapts );
  e_l2.resize( numQuapts );

  for( int qua = 0; qua < numQuapts; ++qua )
  {
    e_r[qua].resize(3);
    e_s[qua].resize(3);
    e_l1[qua].resize(3);
    e_l2[qua].resize(3);
  }
}


FEAElement_Triangle6_membrane::~FEAElement_Triangle6_membrane()
{
  delete [] R; R = nullptr;

  delete [] unx; unx = nullptr; 
  delete [] uny; uny = nullptr; 
  delete [] unz; unz = nullptr;

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
  assert(quad->get_dim() == 3);

  double qua_r, qua_s, qua_t;

  double Rr [nLocBas];
  double Rs [nLocBas];

  double dx_dr[numQuapts], dx_ds[numQuapts];
  double dy_dr[numQuapts], dy_ds[numQuapts];
  double dz_dr[numQuapts], dz_ds[numQuapts];

  std::vector< std::vector<double> > e_a(numQuapts);
  std::vector< std::vector<double> > e_b(numQuapts);

  for( int qua = 0; qua < numQuapts; ++qua )
  {
    e_a[qua].resize(3);
    e_b[qua].resize(3);
  }

  for( int qua = 0; qua < numQuapts; ++qua )
  {
    qua_r = quad -> get_qp( qua, 0 );
    qua_s = quad -> get_qp( qua, 1 );
    qua_t = quad -> get_qp( qua, 2 );

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
    
    dx_dr[qua] = 0.0; dx_ds[qua] = 0.0;
    dy_dr[qua] = 0.0; dy_ds[qua] = 0.0;
    dz_dr[qua] = 0.0; dz_ds[qua] = 0.0;

    for( int ii=0; ii<nLocBas; ++ii )
    {
      dx_dr[qua] += ctrl_x[ii] * Rr[ii];
      dx_ds[qua] += ctrl_x[ii] * Rs[ii];
      
      dy_dr[qua] += ctrl_y[ii] * Rr[ii];
      dy_ds[qua] += ctrl_y[ii] * Rs[ii];
      
      dz_dr[qua] += ctrl_z[ii] * Rr[ii];
      dz_ds[qua] += ctrl_z[ii] * Rs[ii];
    }

    MATH_T::cross3d( dx_dr[qua], dy_dr[qua], dz_dr[qua], 
        dx_ds[qua], dy_ds[qua], dz_ds[qua],
        unx[qua], uny[qua], unz[qua] );
  
    detJac[qua] = MATH_T::normalize3d( unx[qua], uny[qua], unz[qua] );

    // ======= Global-to-local rotation matrix =======
    e_r[qua][0] = dx_dr[qua] / MATH_T::norm2( dx_dr[qua], dy_dr[qua], dz_dr[qua] );
    e_r[qua][1] = dy_dr[qua] / MATH_T::norm2( dx_dr[qua], dy_dr[qua], dz_dr[qua] );
    e_r[qua][2] = dz_dr[qua] / MATH_T::norm2( dx_dr[qua], dy_dr[qua], dz_dr[qua] );

    e_s[qua][0] = dx_ds[qua] / MATH_T::norm2( dx_ds[qua], dy_ds[qua], dz_ds[qua] );
    e_s[qua][1] = dy_ds[qua] / MATH_T::norm2( dx_ds[qua], dy_ds[qua], dz_ds[qua] );
    e_s[qua][2] = dz_ds[qua] / MATH_T::norm2( dx_ds[qua], dy_ds[qua], dz_ds[qua] );

    // e_a = 0.5*(e_r + e_s) / || 0.5*(e_r + e_s) ||
    for( unsigned int ii = 0; ii < e_r[qua].size(); ++ii )
    {
      e_a[qua][ii] = 0.5 * ( e_r[qua][ii] + e_s[qua][ii] );
    }
    MATH_T::normalize3d( e_a[qua][0], e_a[qua][1], e_a[qua][2] );

    // e_b = vec(un) x e_a / || vec(un) x e_a ||
    MATH_T::cross3d(unx[qua], uny[qua], unz[qua], e_a[qua][0], e_a[qua][1], e_a[qua][2],
        e_b[qua][0], e_b[qua][1], e_b[qua][2]);
    MATH_T::normalize3d( e_b[qua][0], e_b[qua][1], e_b[qua][2] );

    // e_l1 = sqrt(2)/2 * (e_a - e_b)
    // e_l2 = sqrt(2)/2 * (e_a + e_b)
    for( unsigned int ii = 0; ii < e_a[qua].size(); ++ii )
    {
      e_l1[qua][ii] = std::sqrt(2.0) / 2.0 * ( e_a[qua][ii] - e_b[qua][ii] );
      e_l2[qua][ii] = std::sqrt(2.0) / 2.0 * ( e_a[qua][ii] + e_b[qua][ii] );
    }

    // Q = transpose([ e_l1, e_l2, un ])
    Q[qua] = Matrix_3x3(e_l1[qua][0], e_l1[qua][1], e_l1[qua][2],
                        e_l2[qua][0], e_l2[qua][1], e_l2[qua][2],
                            unx[qua],     uny[qua],     unz[qua] );

  }
}


void FEAElement_Triangle6_membrane::get_R( 
    const int &quaindex, double * const &basis ) const
{
  assert(quaindex>=0 && quaindex < numQuapts);
  const int offset = quaindex * nLocBas;
  basis[0] = R[offset];
  basis[1] = R[offset+1];
  basis[2] = R[offset+2];
  basis[3] = R[offset+3];
  basis[4] = R[offset+4];
  basis[5] = R[offset+5];
}


void FEAElement_Triangle6_membrane::get_gradR( const int &quaindex,
    double * const &basis_x, double * const &basis_y, double * const &basis_z ) const
{
  // TODO
}


void FEAElement_Triangle6_membrane::get_R_gradR( const int &quaindex, 
    double * const &basis, double * const &basis_x, double * const &basis_y,
    double * const &basis_z ) const
{
  // TODO
}


void FEAElement_Triangle6_membrane::get_2d_normal_out( const int &qua,
        double &nx, double &ny, double &nz, double &len ) const
{
  nx = unx[qua]; ny = uny[qua]; nz = unz[qua];
  len = detJac[qua];
}


void FEAElement_Triangle6_membrane::get_rotationMatrix( const int &quaindex,
    Matrix_3x3 &rot_mat ) const
{
  rot_mat = Q[quaindex];
}


void FEAElement_Triangle6_membrane::get_normal_out( const int &qua,
    const double &sur_pt_x, const double &sur_pt_y, const double &sur_pt_z,
    const double &intpt_x, const double &intpt_y, const double &intpt_z,
    double &nx, double &ny, double &nz, double &len ) const
{
  // Construct a vector starting from the interior point to the triangle
  const double mx = sur_pt_x - intpt_x;
  const double my = sur_pt_y - intpt_y;
  const double mz = sur_pt_z - intpt_z;

  // Do inner product with the normal vector
  const double mdotn = mx * unx[qua] + my * uny[qua] + mz * unz[qua];

  SYS_T::print_fatal_if( std::abs(mdotn) < 1.0e-10, "Warning: FEAElement_Triangle3_membrane::get_normal_out, the element might be ill-shaped.\n");

  if(mdotn < 0)
  {
    nx = (-1.0) * unx[qua];
    ny = (-1.0) * uny[qua];
    nz = (-1.0) * unz[qua];
  }
  else
  {
    nx = unx[qua];
    ny = uny[qua];
    nz = unz[qua];
  }

  len = detJac[qua];
}

// EOF
