#include "FEAElement_Quad9.hpp"

FEAElement_Quad9::FEAElement_Quad9( const int &in_nqua ) : numQuapts( in_nqua )
{
  R = new double [9*numQuapts];
  dR_dx = new double [9*numQuapts];
  dR_dy = new double [9*numQuapts];
  
  Jac = new double [9*numQuapts];

  d2R_dxx = new double [9 * numQuapts];
  d2R_dyy = new double [9 * numQuapts];
  d2R_dxy = new double [9 * numQuapts];
}

FEAElement_Quad9::~FEAElement_Quad9()
{
  delete [] R;             R = nullptr;
  delete [] dR_dx;     dR_dx = nullptr;
  delete [] dR_dy;     dR_dy = nullptr;
  delete [] Jac;         Jac = nullptr;
  delete [] d2R_dxx; d2R_dxx = nullptr;
  delete [] d2R_dyy; d2R_dyy = nullptr;
  delete [] d2R_dxy; d2R_dxy = nullptr;
}

void FEAElement_Quad9::print_info() const
{
  SYS_T::commPrint("Quad9: ");
  SYS_T::commPrint("9-node quadrilateral element with up to 2nd derivatives. \n");
  SYS_T::commPrint("elemType: %d \n", get_Type());
  SYS_T::commPrint("Note: Jacobian and inverse Jacobian are evaluated. \n");
}

double FEAElement_Quad9::get_memory_usage() const
{
  double double_size = 63 * numQuapts;
  double int_size = 1;
  return double_size * 8.0 + int_size * 4.0;
}

void FEAElement_Quad9::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x,
    const double * const &ctrl_y )
{
  ASSERT(quad->get_dim() == 2, "FEAElement_Quad9::buildBasis function error.\n" );

  for( int qua = 0; qua < numQuapts; ++qua )
  {
    const double qua_r = quad -> get_qp( qua, 0 );
    const double qua_s = quad -> get_qp( qua, 1 );

    const int offset = 9 * qua;

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

    const double dR_dr[9] { 
    dNr[0] * Ns[0], dNr[2] * Ns[0], dNr[2] * Ns[2],
    dNr[0] * Ns[2], dNr[1] * Ns[0], dNr[2] * Ns[1],
    dNr[1] * Ns[2], dNr[0] * Ns[1], dNr[1] * Ns[1] };
    
    const double dR_ds[9] { 
    Nr[0] * dNs[0], Nr[2] * dNs[0], Nr[2] * dNs[2],
    Nr[0] * dNs[2], Nr[1] * dNs[0], Nr[2] * dNs[1],
    Nr[1] * dNs[2], Nr[0] * dNs[1], Nr[1] * dNs[1] };

    const double d2R_drs[9] {
    dNr[0] * dNs[0], dNr[2] * dNs[0], dNr[2] * dNs[2],
    dNr[0] * dNs[2], dNr[1] * dNs[0], dNr[2] * dNs[1],
    dNr[1] * dNs[2], dNr[0] * dNs[1], dNr[1] * dNs[1] };

    const double d2Nr[3] { 4.0, -8.0, 4.0 };
    const double d2Ns[3] { 4.0, -8.0, 4.0 };

    const double d2R_drr[9] {
    d2Nr[0] * Ns[0], d2Nr[2] * Ns[0], d2Nr[2] * Ns[2],
    d2Nr[0] * Ns[2], d2Nr[1] * Ns[0], d2Nr[2] * Ns[1],
    d2Nr[1] * Ns[2], d2Nr[0] * Ns[1], d2Nr[1] * Ns[1] };
    
    const double d2R_dss[9] {
    Nr[0] * d2Ns[0], Nr[2] * d2Ns[0], Nr[2] * d2Ns[2],
    Nr[0] * d2Ns[2], Nr[1] * d2Ns[0], Nr[2] * d2Ns[1],
    Nr[1] * d2Ns[2], Nr[0] * d2Ns[1], Nr[1] * d2Ns[1] };
    
    double dx_dr = 0.0, dx_ds = 0.0, dy_dr = 0.0, dy_ds = 0.0;
    double xrr = 0.0, xss = 0.0, yrr = 0.0, yss = 0.0;
    double xrs = 0.0, yrs = 0.0;

    for(int ii=0; ii<9; ++ii)
    {
      dx_dr += ctrl_x[ii] * dR_dr[ii];
      dx_ds += ctrl_x[ii] * dR_ds[ii];
      dy_dr += ctrl_y[ii] * dR_dr[ii];
      dy_ds += ctrl_y[ii] * dR_ds[ii];

      xrr += ctrl_x[ii] * d2R_drr[ii];
      xss += ctrl_x[ii] * d2R_dss[ii];
      yrr += ctrl_y[ii] * d2R_drr[ii];
      yss += ctrl_y[ii] * d2R_dss[ii];

      xrs += ctrl_x[ii] * d2R_drs[ii];
      yrs += ctrl_y[ii] * d2R_drs[ii];
    }
    
    Jac[4*qua+0] = dx_dr;
    Jac[4*qua+1] = dx_ds;
    Jac[4*qua+2] = dy_dr;
    Jac[4*qua+3] = dy_ds;

    Jac[8*numQuapts + qua] = dx_dr * dy_ds - dx_ds * dy_dr;

    const double inv_detJac = 1.0 / Jac[8*numQuapts + qua];

    const double dr_dx = dy_ds * inv_detJac;
    const double dr_dy = - dx_ds * inv_detJac;
    const double ds_dx = - dy_dr * inv_detJac;
    const double ds_dy = dx_dr * inv_detJac;

    Jac[4*numQuapts + 4*qua + 0] = dr_dx;
    Jac[4*numQuapts + 4*qua + 1] = dr_dy;
    Jac[4*numQuapts + 4*qua + 2] = ds_dx;
    Jac[4*numQuapts + 4*qua + 3] = ds_dy;
  
    for(int ii=0; ii<9; ++ii)
    {
      dR_dx[offset+ii] = dR_dr[ii] * dr_dx + dR_ds[ii] * ds_dx;
      dR_dy[offset+ii] = dR_dr[ii] * dr_dy + dR_ds[ii] * ds_dy;
    }

    // Setup the 3x3 matrix
    const double a11 = dx_dr * dx_dr;
    const double a12 = dy_dr * dy_dr;
    const double a13 = 2 * dx_dr * dy_dr;
    const double a21 = dx_dr * dx_ds;
    const double a22 = dy_dr * dy_ds;
    const double a23 = dx_dr * dy_ds + dx_ds * dy_dr;
    const double a31 = dx_ds * dx_ds;
    const double a32 = dy_ds * dy_ds;
    const double a33 = 2 * dx_ds * dy_ds;

    FE_T::Matrix_double_3by3_Array LHS(a11, a12, a13, a21, a22, a23, a31, a32, a33);

    // LU factorization
    LHS.LU_fac();

    for(int ii=0; ii<9; ++ii)
    {
      const std::array<double, 3> RHS {{ d2R_drr[ii] - dR_dx[offset+ii] * xrr - dR_dy[offset+ii] * yrr,
        d2R_drs[ii] - dR_dx[offset+ii] * xrs - dR_dy[offset+ii] * yrs,
        d2R_dss[ii] - dR_dx[offset+ii] * xss - dR_dy[offset+ii] * yss }};

      const auto sol = LHS.LU_solve(RHS);

      d2R_dxx[offset+ii] = sol[0];
      d2R_dyy[offset+ii] = sol[1];
      d2R_dxy[offset+ii] = sol[2];
    }
  }
}

double FEAElement_Quad9::get_h( const double * const &ctrl_x,
    const double * const &ctrl_y ) const
{
  const double diag[2] { std::pow((ctrl_x[0] - ctrl_x[2]), 2.0) 
      + std::pow((ctrl_y[0] - ctrl_y[2]), 2.0),
      std::pow((ctrl_x[1] - ctrl_x[3]), 2.0)
      + std::pow((ctrl_y[1] - ctrl_y[3]), 2.0) };
    
  const double d = (diag[0] > diag[1] ? diag[0] : diag[1]);

  return std::sqrt(d);
}

void FEAElement_Quad9::get_R( const int &quaindex, 
    double * const &basis ) const
{
  const int offset = quaindex * 9;
  for(int ii=0; ii<9; ++ii) basis[ii] = R[offset+ii];
}

std::vector<double> FEAElement_Quad9::get_R( const int &quaindex ) const
{
  const int offset = quaindex * 9;
  std::vector<double> vec(R + offset, R + offset + 9);
  return vec;
}

void FEAElement_Quad9::get_gradR( const int &quaindex, 
    double * const &basis_x, double * const &basis_y ) const
{
  const int offset = quaindex * 9;
  for(int ii=0; ii<9; ++ii)
  {
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
  }
}

void FEAElement_Quad9::get_R_gradR( const int &quaindex, 
    double * const &basis, double * const &basis_x, 
    double * const &basis_y ) const
{
  const int offset = quaindex * 9;
  for(int ii=0; ii<9; ++ii)
  {
    basis[ii]   = R[offset + ii];
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
  }
}

std::vector<double> FEAElement_Quad9::get_dR_dx( const int &quaindex ) const
{
  const int offset = quaindex * 9;
  std::vector<double> vec(dR_dx + offset, dR_dx + offset + 9);
  return vec;
}

std::vector<double> FEAElement_Quad9::get_dR_dy( const int &quaindex ) const
{
  const int offset = quaindex * 9;
  std::vector<double> vec(dR_dy + offset, dR_dy + offset + 9);
  return vec;
}

void FEAElement_Quad9::get_2D_R_dR_d2R( const int &quaindex,
    double * const &basis,
    double * const &basis_x, double * const &basis_y,
    double * const &basis_xx, double * const &basis_yy,
    double * const &basis_xy ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Quad::get_2D_R_dR_d2R function error.\n" );
  const int offset = quaindex * 9;
  for( int ii=0; ii<9; ++ii )
  {
    basis[ii]   = R[offset + ii];
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
    basis_xx[ii] = d2R_dxx[offset + ii];
    basis_yy[ii] = d2R_dyy[offset + ii];
    basis_xy[ii] = d2R_dxy[offset + ii];
  }
}

std::vector<double> FEAElement_Quad9::get_d2R_dxx( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Quad9::get_d2R_dxx function error.\n" );
  const int offset = quaindex * 9;
  std::vector<double> vec(d2R_dxx + offset, d2R_dxx + offset + 9);
  return vec;
}

std::vector<double> FEAElement_Quad9::get_d2R_dyy( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Quad9::get_d2R_dyy function error.\n" );
  const int offset = quaindex * 9;
  std::vector<double> vec(d2R_dyy + offset, d2R_dyy + offset + 9);
  return vec;
}

std::vector<double> FEAElement_Quad9::get_d2R_dxy( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Quad9::get_d2R_dxy function error.\n" );
  const int offset = quaindex * 9;
  std::vector<double> vec(d2R_dxy + offset, d2R_dxy + offset + 9);
  return vec;
}

void FEAElement_Quad9::get_Jacobian(const int &quaindex,
    double * const &jac_value) const
{
  jac_value[0] = Jac[4*quaindex];
  jac_value[1] = Jac[4*quaindex+1];
  jac_value[2] = Jac[4*quaindex+2];
  jac_value[3] = Jac[4*quaindex+3];
}

void FEAElement_Quad9::get_invJacobian(const int &quaindex,
    double * const &jac_value) const
{
  const int offset = 4 * numQuapts + 4 * quaindex;
  jac_value[0] = Jac[offset];
  jac_value[1] = Jac[offset+1];
  jac_value[2] = Jac[offset+2];
  jac_value[3] = Jac[offset+3];
}

// EOF
