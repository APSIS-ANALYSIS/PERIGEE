#include "FE_Tools.hpp"

double FE_T::L2Proj_DGP0( const double * const &f,
    const double * const &gwts, const int &nqp )
{
  double sum_top = 0.0, sum_bot = 0.0;
  for(int ii=0; ii<nqp; ++ii)
  {
    sum_top += f[ii] * gwts[ii];
    sum_bot += gwts[ii];
  }
  return sum_top / sum_bot;
}

void FE_T::L2Proj_DGP1_2D( const double * const &f,
    const double * const &gwts,
    const double * const &qp_x,
    const double * const &qp_y,
    const int &nqp,
    double &coeff_0, double &coeff_x, double &coeff_y )
{
  double a [9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double b [3] = {0.0, 0.0, 0.0};

  for(int ii=0; ii<nqp; ++ii)
  {
    a[0] += gwts[ii];
    a[1] += gwts[ii] * qp_x[ii];
    a[2] += gwts[ii] * qp_y[ii];
    a[4] += gwts[ii] * qp_x[ii] * qp_x[ii];
    a[5] += gwts[ii] * qp_x[ii] * qp_y[ii];
    a[8] += gwts[ii] * qp_y[ii] * qp_y[ii];

    b[0] += gwts[ii] * f[ii];
    b[1] += gwts[ii] * f[ii] * qp_x[ii];
    b[2] += gwts[ii] * f[ii] * qp_y[ii];
  }
  a[3] = a[1]; a[6] = a[2]; a[7] = a[5];

  Matrix_double_3by3_Array AA(a[0], a[1], a[2], a[3],
      a[4], a[5], a[6], a[7], a[8]);

  AA.LU_fac();

  AA.LU_solve(b[0], b[1], b[2], coeff_0, coeff_x, coeff_y);
}

void FE_T::L2Proj_DGP1_3D( const double * const &f,
    const double * const &gwts,
    const double * const &qp_x,
    const double * const &qp_y,
    const double * const &qp_z,
    const int &nqp,
    double &coeff_0, double &coeff_x, double &coeff_y, double &coeff_z )
{
  std::array<double,16> a {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  std::array<double,4> b {{0.0, 0.0, 0.0, 0.0}};

  for(int ii=0; ii<nqp; ++ii)
  {
    a[0] += gwts[ii];
    a[1] += gwts[ii] * qp_x[ii];
    a[2] += gwts[ii] * qp_y[ii];
    a[3] += gwts[ii] * qp_z[ii];
    a[5] += gwts[ii] * qp_x[ii] * qp_x[ii];
    a[6] += gwts[ii] * qp_x[ii] * qp_y[ii];
    a[7] += gwts[ii] * qp_x[ii] * qp_z[ii];
    a[10] += gwts[ii] * qp_y[ii] * qp_y[ii];
    a[11] += gwts[ii] * qp_y[ii] * qp_z[ii];
    a[15] += gwts[ii] * qp_z[ii] * qp_z[ii];

    b[0] += gwts[ii] * f[ii];
    b[1] += gwts[ii] * f[ii] * qp_x[ii];
    b[2] += gwts[ii] * f[ii] * qp_y[ii];
    b[3] += gwts[ii] * f[ii] * qp_z[ii];
  }
  a[4] = a[1]; a[8] = a[2]; a[12] = a[3];
  a[9] = a[6]; a[13] = a[7]; a[14] = a[11];

  MATH_T::Matrix_Dense<4> AA(a);

  AA.LU_fac();

  auto out = AA.LU_solve(b);

  coeff_0 = out[0]; coeff_x = out[1]; coeff_y = out[2]; coeff_z = out[3];
}

void FE_T::get_n_from_t( 
    const double &tx, const double &ty, const double &tz,
    const double &p0_x, const double &p0_y, const double &p0_z,
    const double &p1_x, const double &p1_y, const double &p1_z,
    double &nx, double &ny, double &nz )
{
  const double mx = p0_x - p1_x;
  const double my = p0_y - p1_y;
  const double mz = p0_z - p1_z;

  const double mdt = mx * tx + my * ty + mz * tz;
  const double tdt = tx * tx + ty * ty + tz * tz;
  const double fac = mdt / tdt;

  nx = mx - fac * tx;
  ny = my - fac * ty;
  nz = mz - fac * tz;

  const double len = std::sqrt(nx*nx + ny*ny + nz*nz);
  nx = nx / len;
  ny = ny / len;
  nz = nz / len;
}


Vector_3 FE_T::get_n_from_t( const Vector_3 &tan, const Vector_3 &p0, const Vector_3 &p1 )
{
  const Vector_3 mm = p0 - p1;
  const double mdt = VEC3_T::dot_product( mm, tan );
  const double tdt = VEC3_T::dot_product( tan, tan );
  const double fac = mdt / tdt;

  const Vector_3 nn = fac * tan;

  return VEC3_T::normalize(nn);
}

void FE_T::get_tet_sphere_info( const double &x0, const double &x1,
    const double &x2, const double &x3, const double &y0,
    const double &y1, const double &y2, const double &y3,
    const double &z0, const double &z1, const double &z2,
    const double &z3, double &x, double &y, double &z, double &r )
{
  Matrix_double_3by3_Array AA(
      2.0 * (x1-x0), 2.0 * (y1-y0), 2.0 * (z1-z0),
      2.0 * (x2-x0), 2.0 * (y2-y0), 2.0 * (z2-z0),
      2.0 * (x3-x0), 2.0 * (y3-y0), 2.0 * (z3-z0) );

  AA.LU_fac();

  const double xyz2 = x0*x0 + y0*y0 + z0*z0;

  AA.LU_solve( x1*x1 + y1*y1 + z1*z1 - xyz2,
      x2*x2 + y2*y2 + z2*z2 - xyz2,
      x3*x3 + y3*y3 + z3*z3 - xyz2,
      x, y, z );

  r = std::sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0) );
}

Vector_3 FE_T::get_tet_sphere_info( const Vector_3 &pt0,
    const Vector_3 &pt1, const Vector_3 &pt2, const Vector_3 &pt3, 
    double &radius ) 
{
  Matrix_double_3by3_Array AA(
      2.0 * (pt1.x()-pt0.x()), 2.0 * (pt1.y()-pt0.y()), 2.0 * (pt1.z()-pt0.z()),
      2.0 * (pt2.x()-pt0.x()), 2.0 * (pt2.y()-pt0.y()), 2.0 * (pt2.z()-pt0.z()),
      2.0 * (pt3.x()-pt0.x()), 2.0 * (pt3.y()-pt0.y()), 2.0 * (pt3.z()-pt0.z()) );

  AA.LU_fac();

  const double xyz2 = pt0.dot_product( pt0 );

  const Vector_3 centre = AA.LU_solve( Vector_3( pt1.dot_product(pt1) - xyz2,
      pt2.dot_product(pt2) - xyz2, pt3.dot_product(pt3) - xyz2 ) );

  radius = ( centre - pt0 ).norm2();

  return centre;
}

// EOF
