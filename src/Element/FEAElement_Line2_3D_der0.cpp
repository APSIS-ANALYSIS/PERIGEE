#include "FEAElement_Line2_3D_der0.hpp"

FEAElement_Line2_3D_der0::FEAElement_Line2_3D_der0( 
    const int &in_nqua )
: nLocBas( 2 ), numQuapts( in_nqua )
{
  R = new double [2 * numQuapts];
}


FEAElement_Line2_3D_der0::~FEAElement_Line2_3D_der0()
{
  delete [] R; R = NULL;
}


void FEAElement_Line2_3D_der0::print() const
{
  SYS_T::commPrint("Line2_3D_der0: ");
  SYS_T::commPrint("P1 line element in 3D with no derivative evaluation. \n");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d \n", get_Type());
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
  
  dx_dr = ctrl_x[1] - ctrl_x[0];
  dy_dr = ctrl_y[1] - ctrl_y[0];
  dz_dr = ctrl_z[1] - ctrl_z[0];

  detJac = std::sqrt(dx_dr*dx_dr + dy_dr*dy_dr + dz_dr*dz_dr);
}


void FEAElement_Line2_3D_der0::get_R( const int &quaindex, 
    double * const &basis ) const
{
  basis[0] = R[quaindex*2];
  basis[1] = R[quaindex*2+1];
}


void FEAElement_Line2_3D_der0::get_normal_out( const int &quaindex,
    const double * const &ctrl_x, const double * const &ctrl_y,
    const double * const &ctrl_z,
    const double &intpt_x, const double &intpt_y, const double &intpt_z,
    double &nx, double &ny, double &nz, double &len ) const
{
  const int offset = quaindex * 2;
  double tan_root_x = ctrl_x[0] * R[offset] + ctrl_x[1] * R[offset+1];
  double tan_root_y = ctrl_y[0] * R[offset] + ctrl_y[1] * R[offset+1];
  double tan_root_z = ctrl_z[0] * R[offset] + ctrl_z[1] * R[offset+1];

  MATH_T::get_n_from_t( dx_dr, dy_dr, dz_dr, tan_root_x, tan_root_y,
      tan_root_z, intpt_x, intpt_y, intpt_z, nx, ny, nz );

  len = detJac;
}


// EOF
