#include "FEAElement_Line2_3D_der0.hpp"

FEAElement_Line2_3D_der0::FEAElement_Line2_3D_der0( 
    const int &in_nqua ) : numQuapts( in_nqua )
{
  R = new double [2 * numQuapts];
}

FEAElement_Line2_3D_der0::~FEAElement_Line2_3D_der0()
{
  delete [] R; R = nullptr;
}

void FEAElement_Line2_3D_der0::print_info() const
{
  SYS_T::commPrint("Line2_3D_der0: ");
  SYS_T::commPrint("P1 line element in 3D with no derivative evaluation. \n");
  SYS_T::commPrint("elemType: %d \n", get_Type());
  SYS_T::commPrint("Note: This element is designed for natural BC integrals. \n ");
}

double FEAElement_Line2_3D_der0::get_memory_usage() const
{
  double double_size = 2 * numQuapts + 4;
  double int_size = 2;
  return double_size * 8.0 + int_size * 4.0;
}

void FEAElement_Line2_3D_der0::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x, const double * const &ctrl_y,
    const double * const &ctrl_z )
{
  for( int qua = 0; qua < numQuapts; ++qua )
  {
    const double qua_r = quad -> get_qp( qua );

    R[qua*2 + 0] = 1 - qua_r;
    R[qua*2 + 1] = qua_r;
  }
  
  dx_dr = Vector_3( ctrl_x[1] - ctrl_x[0], ctrl_y[1] - ctrl_y[0], ctrl_z[1] - ctrl_z[0] );

  detJac = dx_dr.norm2();
}

void FEAElement_Line2_3D_der0::get_R( const int &quaindex, 
    double * const &basis ) const
{
  basis[0] = R[quaindex*2];
  basis[1] = R[quaindex*2+1];
}

Vector_3 FEAElement_Line2_3D_der0::get_normal_out( const int &quaindex,
    const std::vector<Vector_3> &ctrl_pt,
    const Vector_3 &int_pt, double &len ) const
{
  const int offset = quaindex * 2;
  const Vector_3 tan_root = R[offset] * ctrl_pt[0] + R[offset+1] * ctrl_pt[1];

  const Vector_3 un = FE_T::get_n_from_t( dx_dr, tan_root, int_pt );

  len = detJac;
  return un;
}

// EOF
