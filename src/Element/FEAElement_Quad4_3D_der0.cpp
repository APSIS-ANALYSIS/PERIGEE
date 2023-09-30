#include "FEAElement_Quad4_3D_der0.hpp"

FEAElement_Quad4_3D_der0::FEAElement_Quad4_3D_der0( const int &in_nqua )
: numQuapts( in_nqua )
{
  R = new double [4*numQuapts];
  detJac = new double [numQuapts];
  un.resize( numQuapts );
}

FEAElement_Quad4_3D_der0::~FEAElement_Quad4_3D_der0()
{
  delete [] R;         R = nullptr;
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
    
    Vector_3 dx_dr( 0.0, 0.0, 0.0 );
    Vector_3 dx_ds( 0.0, 0.0, 0.0 );

    for( int ii=0; ii<4; ++ii )
    {
      dx_dr += Vector_3( ctrl_x[ii] * Rr[ii], ctrl_y[ii] * Rr[ii], ctrl_z[ii] * Rr[ii] );
      dx_ds += Vector_3( ctrl_x[ii] * Rs[ii], ctrl_y[ii] * Rs[ii], ctrl_z[ii] * Rs[ii] );
    }

    un[qua] = Vec3::cross_product( dx_dr, dx_ds );
    detJac[qua] = un[qua].normalize();
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
  return un[qua]; 
}

Vector_3 FEAElement_Quad4_3D_der0::get_normal_out( const int &qua,
    const Vector_3 &sur_pt, const Vector_3 &int_pt, double &len ) const
{
  // Construct a vector starting from the interior point to the triangle
  const Vector_3 mm = sur_pt - int_pt;

  // Do inner product with the normal vector
  const double mdotn = Vec3::dot_product( mm, un[qua] );

  SYS_T::print_fatal_if( std::abs(mdotn) < 1.0e-10, "Warning: FEAElement_Quad4_3D_der0::get_normal_out, the element might be ill-shaped.\n");

  len = detJac[qua];

  if(mdotn < 0) return (-1.0)*un[qua];
  else return un[qua];
}

// EOF
