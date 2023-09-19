#include "FEAElement_Quad4_3D_der0.hpp"
#include "Vector_3.hpp"

FEAElement_Quad4_3D_der0::FEAElement_Quad4_3D_der0( const int &in_nqua )
: numQuapts( in_nqua )
{
  R = new double [4*numQuapts];
  
  unx = new double [numQuapts];
  uny = new double [numQuapts];
  unz = new double [numQuapts];

  detJac = new double [numQuapts];
}

FEAElement_Quad4_3D_der0::~FEAElement_Quad4_3D_der0()
{
  delete [] R;         R = nullptr;
  delete [] unx;     unx = nullptr; 
  delete [] uny;     uny = nullptr;
  delete [] unz;     unz = nullptr;
  delete [] detJac; detJac = nullptr;
}

void FEAElement_Quad4_3D_der0::print_info() const
{
  SYS_T::commPrint("Quad4_3D_der0: ");
  SYS_T::commPrint("4-node quad element with no derivative evaluated. \n ");
  SYS_T::commPrint("elemType: %d. \n", get_Type());
  SYS_T::commPrint("Note: This element is designed for natural BC integrals. \n ");
}

double FEAElement_Quad4_3D_der0::get_memory_usage() const
{
  const double dsize = 8 * numQuapts;
  return dsize * 8.0 + 4.0;
}

void FEAElement_Quad4_3D_der0::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z )
{
  ASSERT(quad->get_dim() == 2, "FEAElement_Quad4_3D_der0::buildBasis function error.\n" );

  for( int qua = 0; qua < numQuapts; ++qua )
  {
    const double qua_r = quad -> get_qp( qua, 0 );
    const double qua_s = quad -> get_qp( qua, 1 );

    const int offset = 4 * qua;

    R[offset + 0] = (1.0 - qua_r) * (1.0 - qua_s);
    R[offset + 1] = qua_r * (1.0 - qua_s);
    R[offset + 2] = qua_r * qua_s;
    R[offset + 3] = (1.0 - qua_r) * qua_s;
 
    const double Rr[4] { qua_s - 1.0, 1.0 - qua_s, qua_s, -qua_s };
    
    const double Rs[4] { qua_r - 1.0, -qua_r, qua_r, 1.0 - qua_r };
    
    double dx_dr = 0.0, dx_ds = 0.0, dy_dr = 0.0, dy_ds = 0.0, dz_dr = 0.0, dz_ds = 0.0;

    for( int ii=0; ii<4; ++ii )
    {
      dx_dr += ctrl_x[ii] * Rr[ii];
      dx_ds += ctrl_x[ii] * Rs[ii];
      
      dy_dr += ctrl_y[ii] * Rr[ii];
      dy_ds += ctrl_y[ii] * Rs[ii];
      
      dz_dr += ctrl_z[ii] * Rr[ii];
      dz_ds += ctrl_z[ii] * Rs[ii];
    }
    
    Vector_3 vec1(dx_dr, dy_dr, dz_dr);
    Vector_3 vec2(dx_ds, dy_ds, dz_ds);
    Vector_3 un = cross_product(vec1, vec2);
  
    detJac[qua] = un.normalize();
    unx[qua] = un.x();
    uny[qua] = un.y();
    unz[qua] = un.z();
  }
}

void FEAElement_Quad4_3D_der0::get_R( 
    const int &quaindex, double * const &basis ) const
{
  ASSERT(quaindex>=0 && quaindex < numQuapts, "FEAElement_Quad4_3D_der0::get_R function error.\n" );
  const int offset = quaindex * 4;
  basis[0] = R[offset];
  basis[1] = R[offset+1];
  basis[2] = R[offset+2];
  basis[3] = R[offset+3];
}

std::vector<double> FEAElement_Quad4_3D_der0::get_R( 
    const int &quaindex ) const
{
  ASSERT(quaindex>=0 && quaindex < numQuapts, "FEAElement_Quad4_3D_der0::get_R function error.\n" );
  const int offset = quaindex * 4;
  return { R[offset], R[offset+1], R[offset+2], R[offset+3] };
}

Vector_3 FEAElement_Quad4_3D_der0::get_2d_normal_out( const int &qua,
    double &area ) const
{
  ASSERT(qua >= 0 && qua < numQuapts, "FEAElement_Quad4_3D_der0::get_2d_normal_out function error.\n" );
  area = detJac[qua];
  return Vector_3( unx[qua], uny[qua], unz[qua] );
}

void FEAElement_Quad4_3D_der0::get_normal_out( const int &qua,
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

  SYS_T::print_fatal_if( std::abs(mdotn) < 1.0e-10, "Warning: FEAElement_Triangle3_3D_der0::get_normal_out, the element might be ill-shaped.\n");

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
