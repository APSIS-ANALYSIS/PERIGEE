#include "FEAElement_Triangle6.hpp"

FEAElement_Triangle6::FEAElement_Triangle6( const int &in_nqua )
: numQuapts( in_nqua )
{
  R = new double [6*numQuapts];
  dR_dx = new double [6*numQuapts];
  dR_dy = new double [6*numQuapts];
  
  Jac = new double [9*numQuapts];

  dR_dr = new double [6*numQuapts];
  dR_ds = new double [6*numQuapts];
}

FEAElement_Triangle6::~FEAElement_Triangle6()
{
  delete [] R;     R = nullptr;
  delete [] dR_dx; dR_dx = nullptr;
  delete [] dR_dy; dR_dy = nullptr;
  delete [] Jac;   Jac = nullptr;
  delete [] dR_dr; dR_dr = nullptr;
  delete [] dR_ds; dR_ds = nullptr;
}

void FEAElement_Triangle6::print() const
{
  SYS_T::commPrint("Tri6: ");
  SYS_T::commPrint("6-node triangle element with up to 2nd derivatives. \n");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d \n", get_Type());
  SYS_T::commPrint("Note: Jacobian and inverse Jacobian are evaluated. \n");
}

double FEAElement_Triangle6::get_memory_usage() const
{
  double double_size = 39 * numQuapts;
  double int_size = 1;
  return double_size * 8.0 + int_size * 4.0;
}

void FEAElement_Triangle6::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x,
    const double * const &ctrl_y )
{
  assert(quad->get_dim() == 3);

  for( int qua = 0; qua < numQuapts; ++qua )
  {
    const double qua_r = quad -> get_qp( qua, 0 );
    const double qua_s = quad -> get_qp( qua, 1 );
    const double qua_t = 1.0 - qua_r - qua_s;

    const int offset = 6 * qua;

    R[offset + 0] = qua_t * (2.0 * qua_t - 1.0);
    R[offset + 1] = qua_r * (2.0 * qua_r - 1.0);
    R[offset + 2] = qua_s * (2.0 * qua_s - 1.0);
    R[offset + 3] = 4.0 * qua_r * qua_t;
    R[offset + 4] = 4.0 * qua_r * qua_s;
    R[offset + 5] = 4.0 * qua_s * qua_t;

    dR_dr[offset+0] = 4.0 * qua_r + 4.0 * qua_s - 3.0;
    dR_dr[offset+1] = 4.0 * qua_r - 1.0;
    dR_dr[offset+2] = 0.0;
    dR_dr[offset+3] = 4.0 - 8.0 * qua_r - 4.0 * qua_s;
    dR_dr[offset+4] = 4.0 * qua_s;
    dR_dr[offset+5] = -4.0 * qua_s;
    
    dR_ds[offset+0] = 4.0 * qua_r + 4.0 * qua_s - 3.0;
    dR_ds[offset+1] = 0.0;
    dR_ds[offset+2] = 4.0 * qua_s - 1.0;
    dR_ds[offset+3] = -4.0 * qua_r;
    dR_ds[offset+4] = 4.0 * qua_r;
    dR_ds[offset+5] = 4.0 - 4.0 * qua_r - 8.0 * qua_s;  
    
    double dx_dr = 0.0, dx_ds = 0.0, dy_dr = 0.0, dy_ds = 0.0;
    for(int ii=0; ii<6; ++ii)
    {
      dx_dr += ctrl_x[ii] * dR_dr[offset+ii];
      dx_ds += ctrl_x[ii] * dR_ds[offset+ii];
      dy_dr += ctrl_y[ii] * dR_dr[offset+ii];
      dy_ds += ctrl_y[ii] * dR_ds[offset+ii];
    }
    
    Jac[4*qua+0] = dx_dr;
    Jac[4*qua+1] = dx_ds;
    Jac[4*qua+2] = dy_dr;
    Jac[4*qua+3] = dy_ds;

    Jac[8*numQuapts + qua] = dx_dr * dy_ds - dx_ds * dy_dr;

    const double inv_detJac = 1.0 / Jac[8*numQuapts + qua];

    const double dr_dx = dy_ds * inv_detJac;
    const double dr_dy = (-1.0) * dx_ds * inv_detJac;
    const double ds_dx = (-1.0) * dy_dr * inv_detJac;
    const double ds_dy = dx_dr * inv_detJac;

    Jac[4*numQuapts + 4*qua + 0] = dr_dx;
    Jac[4*numQuapts + 4*qua + 1] = dr_dy;
    Jac[4*numQuapts + 4*qua + 2] = ds_dx;
    Jac[4*numQuapts + 4*qua + 3] = ds_dy;
  
    for(int ii=0; ii<6; ++ii)
    {
      dR_dx[offset+ii] = dR_dr[offset+ii] * dr_dx + dR_ds[offset+ii] * ds_dx;
      dR_dy[offset+ii] = dR_dr[offset+ii] * dr_dy + dR_ds[offset+ii] * ds_dy;
    }
  }
}

double FEAElement_Triangle6::get_h( const double * const &ctrl_x,
    const double * const &ctrl_y ) const
{
  const double a = 2.0 * ( ctrl_x[1] - ctrl_x[0] );
  const double b = 2.0 * ( ctrl_y[1] - ctrl_y[0] );
  const double c = 2.0 * ( ctrl_x[2] - ctrl_x[0] );
  const double d = 2.0 * ( ctrl_y[2] - ctrl_y[0] );

  const double m = ctrl_x[1] * ctrl_x[1] + ctrl_y[1] * ctrl_y[1]
    - ctrl_x[0] * ctrl_x[0] - ctrl_y[0] * ctrl_y[0];

  const double n = ctrl_x[2] * ctrl_x[2] + ctrl_y[2] * ctrl_y[2]
    - ctrl_x[0] * ctrl_x[0] - ctrl_y[0] * ctrl_y[0];


  const double xc = (d * m - b * n) / (a*d-b*c);
  const double yc = (a * n - c * m) / (a*d-b*c);

  const double radius = std::sqrt( (xc-ctrl_x[0])*(xc-ctrl_x[0])
      + (yc - ctrl_y[0])*(yc - ctrl_y[0]) );

  return 2.0 * radius;
}

void FEAElement_Triangle6::get_R( const int &quaindex, 
    double * const &basis ) const
{
  const int offset = quaindex * 6;
  basis[0] = R[offset];   basis[1] = R[offset+1];
  basis[2] = R[offset+2]; basis[3] = R[offset+3];
  basis[4] = R[offset+4]; basis[5] = R[offset+5];
}

std::vector<double> FEAElement_Triangle6::get_R( const int &quaindex ) const
{
  const int offset = quaindex * 6;
  return { R[offset], R[offset+1], R[offset+2], R[offset+3], R[offset+4], R[offset+5] };
}

void FEAElement_Triangle6::get_gradR( const int &quaindex, 
    double * const &basis_x, double * const &basis_y ) const
{
  const int offset = quaindex * 6;
  for(int ii=0; ii<6; ++ii)
  {
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
  }
}

void FEAElement_Triangle6::get_R_gradR( const int &quaindex, 
    double * const &basis, double * const &basis_x, 
    double * const &basis_y ) const
{
  const int offset = quaindex * 6;
  for(int ii=0; ii<6; ++ii)
  {
    basis[ii]   = R[offset + ii];
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
  }
}

std::vector<double> FEAElement_Triangle6::get_dR_dx( const int &quaindex ) const
{
  const int offset = quaindex * 6;
  return { dR_dx[offset], dR_dx[offset+1], dR_dx[offset+2], dR_dx[offset+3], dR_dx[offset+4], dR_dx[offset+5] };
}

std::vector<double> FEAElement_Triangle6::get_dR_dy( const int &quaindex ) const
{
  const int offset = quaindex * 6;
  return { dR_dy[offset], dR_dy[offset+1], dR_dy[offset+2], dR_dy[offset+3], dR_dy[offset+4], dR_dy[offset+5] };
}

void FEAElement_Triangle6::get_2D_R_dR_d2R( const int &quaindex,
    double * const &basis,
    double * const &basis_x, double * const &basis_y,
    double * const &basis_xx, double * const &basis_yy,
    double * const &basis_xy ) const
{
  SYS_T::print_fatal("Error: FEAElement_Triangle6::get_2D_R_dR_d2R is not implemented. \n");
}

void FEAElement_Triangle6::get_Jacobian(const int &quaindex,
    double * const &jac_value) const
{
  jac_value[0] = Jac[4*quaindex];
  jac_value[1] = Jac[4*quaindex+1];
  jac_value[2] = Jac[4*quaindex+2];
  jac_value[3] = Jac[4*quaindex+3];
}

void FEAElement_Triangle6::get_invJacobian(const int &quaindex,
    double * const &jac_value) const
{
  const int offset = 4 * numQuapts + 4 * quaindex;
  jac_value[0] = Jac[offset];
  jac_value[1] = Jac[offset+1];
  jac_value[2] = Jac[offset+2];
  jac_value[3] = Jac[offset+3];
}

// EOF
