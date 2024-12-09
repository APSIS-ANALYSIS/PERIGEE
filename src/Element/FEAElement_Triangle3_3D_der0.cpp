#include "FEAElement_Triangle3_3D_der0.hpp"

FEAElement_Triangle3_3D_der0::FEAElement_Triangle3_3D_der0( 
    const int &in_nqua ) : numQuapts( in_nqua )
{
  R = new double [ 3 * numQuapts ];
}

FEAElement_Triangle3_3D_der0::~FEAElement_Triangle3_3D_der0()
{
  delete [] R; R = nullptr;
}

void FEAElement_Triangle3_3D_der0::print_info() const
{
  SYS_T::commPrint("Triangle3_3D_der0: ");
  SYS_T::commPrint("3-node triangle element with no derivative evaluated. \n ");
  SYS_T::commPrint("Note: This element is designed for natural BC integrals. \n ");
}

double FEAElement_Triangle3_3D_der0::get_memory_usage() const
{
  double double_size = 3 * numQuapts + 10.0;
  double int_size = 2;
  return double_size * 8.0 + int_size * 4.0;
}

void FEAElement_Triangle3_3D_der0::buildBasis( const IQuadPts * const &quad,
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

  const Vector_3 dx_dr( - ctrl_x[0] + ctrl_x[1],
    - ctrl_y[0] + ctrl_y[1],
    - ctrl_z[0] + ctrl_z[1]);

  const Vector_3 dx_ds( - ctrl_x[0] + ctrl_x[2],
    - ctrl_y[0] + ctrl_y[2],
    - ctrl_z[0] + ctrl_z[2]);

  // vec(un) = vec(dx_dr) x vec(dx_ds)
  un = Vec3::cross_product( dx_dr, dx_ds );

  // area = || vec(un) ||
  detJac = un.normalize();
}

void FEAElement_Triangle3_3D_der0::get_R( const int &quaindex, 
    double * const &basis ) const
{
  ASSERT(quaindex>=0 && quaindex < numQuapts, "FEAElement_Triangle3_3D_der0::get_R function error.\n" );
  const int offset = quaindex * 3;
  basis[0] = R[offset];
  basis[1] = R[offset+1];
  basis[2] = R[offset+2];
}

std::vector<double> FEAElement_Triangle3_3D_der0::get_R( const int &quaindex ) const
{
  ASSERT(quaindex>=0 && quaindex < numQuapts, "FEAElement_Triangle3_3D_der0::get_R function error.\n" );
  const int offset = quaindex * 3;
  return { R[offset], R[offset+1], R[offset+2] };
}

Vector_3 FEAElement_Triangle3_3D_der0::get_2d_normal_out( const int &quaindex,
    double &area ) const
{
  ASSERT(quaindex>=0 && quaindex < numQuapts, "FEAElement_Triangle3_3D_der0::get_2d_normal_out function error.\n" );
  area = detJac;
  return un;
}

Vector_3 FEAElement_Triangle3_3D_der0::get_normal_out( const int &quaindex,
    const Vector_3 &sur_pt, const Vector_3 &int_pt, double &area ) const
{
  // Construct a vector from the interior point to the triangle first node
  const Vector_3 mm = sur_pt - int_pt;

  // Dot product of the defined vector with the calculated normal vector
  const double mdotn = Vec3::dot_product( mm, un );

  SYS_T::print_fatal_if( std::abs(mdotn) < 1.0e-10, "Warning: FEAElement_Triangle3_3D_der0::get_normal_out, the element might be ill-shaped.\n");

  area = detJac;

  // If dot product is negative, adjust the normal vector
  if(mdotn < 0) return (-1.0)*un;
  else return un;
}

// EOF
