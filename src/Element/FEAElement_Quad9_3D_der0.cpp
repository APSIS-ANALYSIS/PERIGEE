#include "FEAElement_Quad9_3D_der0.hpp"

FEAElement_Quad9_3D_der0::FEAElement_Quad9_3D_der0( const int &in_nqua )
: numQuapts( in_nqua )
{
  R.resize(nLocBas * numQuapts, 0.0);
  detJac.resize(numQuapts, 0.0);
  un.resize(numQuapts);
}

void FEAElement_Quad9_3D_der0::print_info() const
{
  SYS_T::commPrint("Quad9_3D_der0: ");
  SYS_T::commPrint("Nine-node quad element with no derivative evaluated.\n");
  SYS_T::commPrint("Note: This element is designed for natural BC integrals.\n");
}

void FEAElement_Quad9_3D_der0::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z )
{
  ASSERT(quad->get_dim() == 2, "FEAElement_Quad9_3D_der0::buildBasis function error.\n" );

  for( int qua = 0; qua < numQuapts; ++qua )
  {
    const double qua_r = quad -> get_qp( qua, 0 );
    const double qua_s = quad -> get_qp( qua, 1 );
    
    const int offset = nLocBas * qua;
    
    const double Nr[3] = { (2.0 * qua_r - 1.0) * (qua_r - 1.0),
        - 4.0 * qua_r * (qua_r - 1.0), qua_r * (2.0 * qua_r - 1.0) };
    const double Ns[3] = { (2.0 * qua_s - 1.0) * (qua_s - 1.0),
        - 4.0 * qua_s * (qua_s - 1.0), qua_s * (2.0 * qua_s - 1.0) };
    
    R[offset    ] = Nr[0] * Ns[0];
    R[offset + 1] = Nr[2] * Ns[0];
    R[offset + 2] = Nr[2] * Ns[2];
    R[offset + 3] = Nr[0] * Ns[2];
    R[offset + 4] = Nr[1] * Ns[0];
    R[offset + 5] = Nr[2] * Ns[1];
    R[offset + 6] = Nr[1] * Ns[2];
    R[offset + 7] = Nr[0] * Ns[1];
    R[offset + 8] = Nr[1] * Ns[1];
 
    const double dNr[3] = { 4.0 * qua_r - 3.0, 
        - 8.0 * qua_r + 4.0, 4.0 * qua_r - 1.0 };
    const double dNs[3] = { 4.0 * qua_s - 3.0, 
        - 8.0 * qua_s + 4.0, 4.0 * qua_s - 1.0 };

    const double Rr[9] { 
    dNr[0] * Ns[0], dNr[2] * Ns[0], dNr[2] * Ns[2],
    dNr[0] * Ns[2], dNr[1] * Ns[0], dNr[2] * Ns[1],
    dNr[1] * Ns[2], dNr[0] * Ns[1], dNr[1] * Ns[1] };
    const double Rs[9] { 
    Nr[0] * dNs[0], Nr[2] * dNs[0], Nr[2] * dNs[2],
    Nr[0] * dNs[2], Nr[1] * dNs[0], Nr[2] * dNs[1],
    Nr[1] * dNs[2], Nr[0] * dNs[1], Nr[1] * dNs[1] };
    
    Vector_3 dx_dr( 0.0, 0.0, 0.0 );
    Vector_3 dx_ds( 0.0, 0.0, 0.0 );

    for( int ii=0; ii<nLocBas; ++ii )
    {
      dx_dr += Vector_3( ctrl_x[ii] * Rr[ii], ctrl_y[ii] * Rr[ii], ctrl_z[ii] * Rr[ii] );
      dx_ds += Vector_3( ctrl_x[ii] * Rs[ii], ctrl_y[ii] * Rs[ii], ctrl_z[ii] * Rs[ii] );
    }

    un[qua] = Vec3::cross_product( dx_dr, dx_ds );
    detJac[qua] = un[qua].normalize();
  }
}

void FEAElement_Quad9_3D_der0::get_R( 
    const int &quaindex, double * const &basis ) const
{
  ASSERT(quaindex>=0 && quaindex < numQuapts, "FEAElement_Quad9_3D_der0::get_R function error.\n" );
  const int offset = quaindex * nLocBas;
  for(int ii=0; ii<nLocBas; ++ii) basis[ii] = R[offset+ii];
}

std::vector<double> FEAElement_Quad9_3D_der0::get_R( 
    const int &quaindex ) const
{
  ASSERT(quaindex>=0 && quaindex < numQuapts, "FEAElement_Quad9_3D_der0::get_R function error.\n" );
  const int offset = quaindex * nLocBas;
  std::vector<double> vec(R.begin() + offset, R.begin() + offset + 9);
  return vec;
}

Vector_3 FEAElement_Quad9_3D_der0::get_2d_normal_out( const int &qua,
    double &area ) const
{
  ASSERT(qua >= 0 && qua < numQuapts, "FEAElement_Quad9_3D_der0::get_2d_normal_out function error.\n" );
  area = detJac[qua];
  return un[qua]; 
}

Vector_3 FEAElement_Quad9_3D_der0::get_normal_out( const int &qua,
    const Vector_3 &sur_pt, const Vector_3 &int_pt, double &len ) const
{
  // Construct a vector starting from the interior point to the triangle
  const Vector_3 mm = sur_pt - int_pt;

  // Do inner product with the normal vector
  const double mdotn = Vec3::dot_product( mm, un[qua] );

  SYS_T::print_fatal_if( std::abs(mdotn) < 1.0e-10, "Warning: FEAElement_Quad9_3D_der0::get_normal_out, the element might be ill-shaped.\n");

  len = detJac[qua];

  if(mdotn < 0) return (-1.0)*un[qua];
  else return un[qua];
}

// EOF
