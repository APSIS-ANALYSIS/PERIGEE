#include "FEAElement_Line3_3D_der0.hpp"

FEAElement_Line3_3D_der0::FEAElement_Line3_3D_der0( const int &in_nqua ) : numQuapts( in_nqua )
{
  R      = new double [3*numQuapts];
  dx_dr  = new double [numQuapts];
  dy_dr  = new double [numQuapts];
  dz_dr  = new double [numQuapts];
  detJac = new double [numQuapts];
}

FEAElement_Line3_3D_der0::~FEAElement_Line3_3D_der0()
{
  delete [] R;      R      = nullptr;
  delete [] dx_dr;  dx_dr  = nullptr;
  delete [] dy_dr;  dy_dr  = nullptr;
  delete [] dz_dr;  dz_dr  = nullptr;
  delete [] detJac; detJac = nullptr;
}

void FEAElement_Line3_3D_der0::print_info() const
{
  SYS_T::commPrint("Line3_3D_der0: ");
  SYS_T::commPrint("P2 line element in 3D with no derivative evaluation. \n");
  SYS_T::commPrint("elemType: %d \n", get_Type());
  SYS_T::commPrint("Note: This element is designed for natural BC integrals. \n ");
}

double FEAElement_Line3_3D_der0::get_memory_usage() const
{
  double double_size = 7.0 * numQuapts;
  double int_size = 2.0;
  return double_size * 8.0 + int_size * 4.0;
}

void FEAElement_Line3_3D_der0::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x, const double * const &ctrl_y,
    const double * const &ctrl_z )
{
  for( int qua = 0; qua < numQuapts; ++qua )
  {
    const double qua_r = quad -> get_qp(qua);
    R[qua*3+0] = (2.0 * qua_r - 1.0) * (qua_r - 1.0);
    R[qua*3+1] = qua_r * (2.0 * qua_r - 1.0);
    R[qua*3+2] = 4.0 * (1.0 - qua_r) * qua_r;

    const double N0r = 4.0 * qua_r - 3.0;
    const double N1r = 4.0 * qua_r - 1.0;
    const double N2r = -8.0 * qua_r + 4.0;
    dx_dr[qua] = ctrl_x[0] * N0r + ctrl_x[1] * N1r + ctrl_x[2] * N2r;
    dy_dr[qua] = ctrl_y[0] * N0r + ctrl_y[1] * N1r + ctrl_y[2] * N2r;
    dz_dr[qua] = ctrl_z[0] * N0r + ctrl_z[1] * N1r + ctrl_z[2] * N2r;
    detJac[qua] = std::sqrt( dx_dr[qua] * dx_dr[qua]
        + dy_dr[qua] * dy_dr[qua] + dz_dr[qua] * dz_dr[qua] );
  }
}

void FEAElement_Line3_3D_der0::get_R( const int &quaindex, 
    double * const &basis ) const
{
  basis[0] = R[quaindex*3];
  basis[1] = R[quaindex*3+1];
  basis[2] = R[quaindex*3+2];
}

void FEAElement_Line3_3D_der0::get_1d_normal_out( const int &quaindex,
    const double * const &ctrl_x, const double * const &ctrl_y,
    const double * const &ctrl_z,
    const double &intpt_x, const double &intpt_y, const double &intpt_z,
    double &nx, double &ny, double &nz, double &len ) const
{
  const int offset = quaindex * 3;
  double tan_root_x = ctrl_x[0] * R[offset] + ctrl_x[1] * R[offset+1]
    + ctrl_x[2] * R[offset+2];
  double tan_root_y = ctrl_y[0] * R[offset] + ctrl_y[1] * R[offset+1]
    + ctrl_y[2] * R[offset+2];
  double tan_root_z = ctrl_z[0] * R[offset] + ctrl_z[1] * R[offset+1]
    + ctrl_z[2] * R[offset+2];

  FE_T::get_n_from_t( dx_dr[quaindex], dy_dr[quaindex], dz_dr[quaindex], 
      tan_root_x, tan_root_y, tan_root_z, intpt_x, intpt_y, intpt_z, 
      nx, ny, nz );

  len = detJac[quaindex];
}

// EOF
