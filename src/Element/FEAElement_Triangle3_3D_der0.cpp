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
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d. \n", get_Type());
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

  const double dx_dr = ctrl_x[0] * (-1.0) + ctrl_x[1];
  const double dy_dr = ctrl_y[0] * (-1.0) + ctrl_y[1];
  const double dz_dr = ctrl_z[0] * (-1.0) + ctrl_z[1];

  const double dx_ds = ctrl_x[0] * (-1.0) + ctrl_x[2];
  const double dy_ds = ctrl_y[0] * (-1.0) + ctrl_y[2];
  const double dz_ds = ctrl_z[0] * (-1.0) + ctrl_z[2];

  // vec(un) = vec(dx_dr) x vec(dx_ds)
  MATH_T::cross3d(dx_dr, dy_dr, dz_dr, dx_ds, dy_ds, dz_ds, unx, uny, unz);
  
  // area = || vec(un) ||
  detJac = MATH_T::normalize3d( unx, uny, unz );
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
  return Vector_3( unx, uny, unz );
}

void FEAElement_Triangle3_3D_der0::get_normal_out( const int &quaindex,
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

  SYS_T::print_fatal_if( std::abs(mdotn) < 1.0e-10, "Warning: FEAElement_Triangle3_3D_der0::get_normal_out, the element might be ill-shaped.\n");

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
