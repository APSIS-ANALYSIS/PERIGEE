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

  FE_T::Matrix_double_3by3_Array AA(a[0], a[1], a[2], a[3],
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

Vector_3 FE_T::get_n_from_t( const Vector_3 &tan, const Vector_3 &p0, const Vector_3 &p1 )
{
  const Vector_3 mm = p0 - p1;
  const double mdt = Vec3::dot_product( mm, tan );
  const double tdt = Vec3::dot_product( tan, tan );
  const double fac = mdt / tdt;

  const Vector_3 nn = mm - fac * tan;

  return Vec3::normalize(nn);
}

void FE_T::get_tet_sphere_info(
    const double &x0, const double &x1, const double &x2, const double &x3, 
    const double &y0, const double &y1, const double &y2, const double &y3,
    const double &z0, const double &z1, const double &z2, const double &z3,
    double &xx, double &yy, double &zz, double &rr )
{
  FE_T::Matrix_double_3by3_Array AA(
      2.0 * (x1-x0), 2.0 * (y1-y0), 2.0 * (z1-z0),
      2.0 * (x2-x0), 2.0 * (y2-y0), 2.0 * (z2-z0),
      2.0 * (x3-x0), 2.0 * (y3-y0), 2.0 * (z3-z0) );

  AA.LU_fac();

  const double xyz2 = x0*x0 + y0*y0 + z0*z0;

  AA.LU_solve( x1*x1 + y1*y1 + z1*z1 - xyz2,
      x2*x2 + y2*y2 + z2*z2 - xyz2,
      x3*x3 + y3*y3 + z3*z3 - xyz2,
      xx, yy, zz );

  rr = std::sqrt( (xx-x0)*(xx-x0) + (yy-y0)*(yy-y0) + (zz-z0)*(zz-z0) );
}

double FE_T::get_tet_sphere_radius (
    const double &x0, const double &x1, const double &x2, const double &x3,
    const double &y0, const double &y1, const double &y2, const double &y3,
    const double &z0, const double &z1, const double &z2, const double &z3 )
{
  FE_T::Matrix_double_3by3_Array AA(
      2.0 * (x1-x0), 2.0 * (y1-y0), 2.0 * (z1-z0),
      2.0 * (x2-x0), 2.0 * (y2-y0), 2.0 * (z2-z0),
      2.0 * (x3-x0), 2.0 * (y3-y0), 2.0 * (z3-z0) );

  AA.LU_fac();

  const double xyz2 = x0*x0 + y0*y0 + z0*z0;
  double xx {0.0}, yy {0.0}, zz {0.0};

  AA.LU_solve( x1*x1 + y1*y1 + z1*z1 - xyz2,
      x2*x2 + y2*y2 + z2*z2 - xyz2,
      x3*x3 + y3*y3 + z3*z3 - xyz2,
      xx, yy, zz );

  return std::sqrt( (xx-x0)*(xx-x0) + (yy-y0)*(yy-y0) + (zz-z0)*(zz-z0) );
}

Vector_3 FE_T::get_tet_sphere_info( const Vector_3 &pt0,
    const Vector_3 &pt1, const Vector_3 &pt2, const Vector_3 &pt3, 
    double &radius ) 
{
  FE_T::Matrix_double_3by3_Array AA(
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

double FE_T::get_circumradius( const std::array<Vector_3, 4> &pts )
{
  FE_T::Matrix_double_3by3_Array AA(
      2.0 * (pts[1].x()-pts[0].x()), 2.0 * (pts[1].y()-pts[0].y()), 2.0 * (pts[1].z()-pts[0].z()),
      2.0 * (pts[2].x()-pts[0].x()), 2.0 * (pts[2].y()-pts[0].y()), 2.0 * (pts[2].z()-pts[0].z()),
      2.0 * (pts[3].x()-pts[0].x()), 2.0 * (pts[3].y()-pts[0].y()), 2.0 * (pts[3].z()-pts[0].z()) );

  AA.LU_fac();

  const double xyz2 = pts[0].dot_product( pts[0] );

  const Vector_3 centre = AA.LU_solve( Vector_3( pts[1].dot_product(pts[1]) - xyz2,
        pts[2].dot_product(pts[2]) - xyz2, pts[3].dot_product(pts[3]) - xyz2 ) );

  return ( centre - pts[0] ).norm2();
}

bool FE_T::search_closest_point( const Vector_3 &target_xyz,
    FEAElement * const &elements,
    const double * const &electrl_x,
    const double * const &electrl_y,
    const double * const &electrl_z,
    IQuadPts * const &closest_point )
{
  // initial value
  elements->buildBasis(closest_point, electrl_x, electrl_y, electrl_z);

  const int nLocBas {elements->get_nLocBas()};

  auto R = elements->get_R(0); // only use the first point of the user-defined IQuadPts

  Vector_3 point_xyz(0.0, 0.0, 0.0);

  // physical xyz of the initial point
  for(int ii=0; ii<nLocBas; ++ii)
  {
    point_xyz(0) += R[ii] * electrl_x[ii];
    point_xyz(1) += R[ii] * electrl_y[ii];
    point_xyz(2) += R[ii] * electrl_z[ii];
  }

  // initial distance
  const double init_dist = Vec3::dist(point_xyz, target_xyz);

  // SYS_T::commPrint("      init_dist: %e\n", init_dist);
  if (init_dist < 1e-9) return true;  // lucky enouugh 
  if (init_dist > 8) return false;

  // initialize the nonlinear iteration
  double curr_dist = init_dist;
  double old_dist = curr_dist;
  int iter_counter = 0;

  const double eps = 1.0e-7;

  while(iter_counter < 20)
  {
    // SYS_T::commPrint("      iter = %d\n",iter_counter);
    // only use the first point of the user-defined IQuadPts
    const Vector_3 dx_dr   = elements->get_dx_dr(0, electrl_x, electrl_y, electrl_z);
    const Vector_3 dx_ds   = elements->get_dx_ds(0, electrl_x, electrl_y, electrl_z);
    const Vector_3 d2x_drr = elements->get_d2x_drr(0, electrl_x, electrl_y, electrl_z);
    const Vector_3 d2x_dss = elements->get_d2x_dss(0, electrl_x, electrl_y, electrl_z);
    const Vector_3 d2x_drs = elements->get_d2x_drs(0, electrl_x, electrl_y, electrl_z);

    // Distance funtion: (x_opp - x_tar)^2 + (y_opp - y_tar)^2 + (z_opp - z_tar)^2,
    // x_opp = x_opp(r, s), y_opp = y_opp(r, s), z_opp = z_opp(r, s)
    // To minimize Dist, we need d(Dist)/dr = 0 and d(Dist)/ds = 0.

    const Vector_3 dist_xyz = point_xyz - target_xyz;

    // d(Dist/dr) = 2 * d(x_opp)/dr * (x_opp - x_tar)
    //            + 2 * d(y_opp)/dr * (y_opp - y_tar)
    //            + 2 * d(z_opp)/dr * (z_opp - z_tar)
    // Let Res_r (Residual of 'r' component) = 0.5 * d(Dist/dr).
    const double Res_r = Vec3::dot_product(dx_dr, dist_xyz);

    // Newton-Raphson iteration:
    // d(Res_r)/dr * dr + d(Res_r)/ds * ds = - Res_r
    // Let dRes_r_dr = d(Res_r)/dr = d2(x_opp)/drr * (x_opp - x_tar) + (d(x_opp)/dr)^2
    //                             + d2(y_opp)/drr * (y_opp - y_tar) + (d(y_opp)/dr)^2
    //                             + d2(z_opp)/drr * (z_opp - x_tar) + (d(z_opp)/dr)^2.
    const double dRes_r_dr = Vec3::dot_product( d2x_drr, dist_xyz ) + Vec3::dot_product(dx_dr, dx_dr);

    // Let dRes_r_ds = d(Res_r)/ds = d2(x_opp)/drs * (x_opp - x_tar) + d(x_opp)/dr * d(x_opp)/ds
    //                             + d2(y_opp)/drs * (y_opp - y_tar) + d(y_opp)/dr * d(y_opp)/ds
    //                             + d2(z_opp)/drs * (z_opp - x_tar) + d(z_opp)/dr * d(z_opp)/ds.
    const double dRes_r_ds = Vec3::dot_product( d2x_drs, dist_xyz ) + Vec3::dot_product(dx_dr, dx_ds);

    // Same for Res_s.
    const double Res_s = Vec3::dot_product( dx_ds, dist_xyz );

    const double dRes_s_ds = Vec3::dot_product( d2x_dss, dist_xyz ) + Vec3::dot_product( dx_ds, dx_ds );

    // We use the fact that d(Res_s)/dr = d(Res_r)/ds

    // Solve [dRes_r_dr  dRes_r_ds  {dr  = {-Res_r
    //        dRes_s_dr  dRes_s_ds]  ds}    -Res_s}

    // Ax = b, A = [a b   inv(A) = 1/(ad - bc) * [d -b      x = inv(A) * b
    //              c d]                         -c  a]
    const double ad_bc = dRes_r_dr * dRes_s_ds - dRes_r_ds * dRes_r_ds;
    const double dxi   = (-dRes_s_ds * Res_r + dRes_r_ds * Res_s) / ad_bc;
    const double deta  = (dRes_r_ds * Res_r - dRes_r_dr * Res_s) / ad_bc;

    // SYS_T::commPrint("      dr = %e, ds = %e\n", dxi, deta);

    // Update the xi value
    const double new_xi  = closest_point->get_qp(0, 0) + dxi;
    const double new_eta = closest_point->get_qp(0, 1) + deta;

    closest_point->set_qp( new_xi, new_eta );
    // SYS_T::commPrint("      current [r,s,t] = [%e, %e, %e]\n", closest_point->get_qp(0, 0), closest_point->get_qp(0, 1), closest_point->get_qp(0, 2));

    // Update basis function and physical xyz
    elements->buildBasis(closest_point, electrl_x, electrl_y, electrl_z);

    R = elements->get_R(0);

    point_xyz = Vector_3(0.0, 0.0, 0.0);
    for(int ii=0; ii<nLocBas; ++ii)
    {
      point_xyz(0) += R[ii] * electrl_x[ii];
      point_xyz(1) += R[ii] * electrl_y[ii];
      point_xyz(2) += R[ii] * electrl_z[ii];
    }

    curr_dist = Vec3::dist(point_xyz, target_xyz);
    // SYS_T::commPrint("      curr_dist: %e\n", curr_dist);
    // SYS_T::commPrint("      old_dist: %e\n", old_dist);
    if(std::abs(curr_dist - old_dist) < eps)
    {
      // SYS_T::commPrint("        Cannot optimize further.\n");
      break;
    }

    iter_counter += 1;
    old_dist = curr_dist;
  }

  if( closest_point->check_qp_bound() == false )
    return false;
  else
    return true;
}

namespace FE_T
{
  Matrix_double_3by3_Array::Matrix_double_3by3_Array()
  {
    mat[0] = 1.0; mat[1] = 0.0; mat[2] = 0.0;
    mat[3] = 0.0; mat[4] = 1.0; mat[5] = 0.0;
    mat[6] = 0.0; mat[7] = 0.0; mat[8] = 1.0;

    pp[0] = 0; pp[1] = 1; pp[2] = 2;
  }

  Matrix_double_3by3_Array::Matrix_double_3by3_Array( 
      const double &a11, const double &a12, const double &a13,
      const double &a21, const double &a22, const double &a23,
      const double &a31, const double &a32, const double &a33 )
  {
    mat[0] = a11;  mat[1] = a12;  mat[2] = a13;
    mat[3] = a21;  mat[4] = a22;  mat[5] = a23;
    mat[6] = a31;  mat[7] = a32;  mat[8] = a33;

    pp[0] = 0; pp[1] = 1; pp[2] = 2; 
  }

  Matrix_double_3by3_Array& Matrix_double_3by3_Array::operator= (
      const Matrix_double_3by3_Array &input )
  {
    if(this != &input)
    {
      for(int ii=0; ii<9; ++ii) mat[ii] = input.mat[ii];
      pp[0] = input.pp[0];
      pp[1] = input.pp[1];
      pp[2] = input.pp[2];
    }
    return *this;
  }

  void Matrix_double_3by3_Array::gen_id()
  {
    mat[0] = 1.0;  mat[1] = 0.0;  mat[2] = 0.0;
    mat[3] = 0.0;  mat[4] = 1.0;  mat[5] = 0.0;
    mat[6] = 0.0;  mat[7] = 0.0;  mat[8] = 1.0;

    pp[0] = 0; pp[1] = 1; pp[2] = 2;
  }

  void Matrix_double_3by3_Array::gen_rand(const double &min, const double &max)
  {
    std::random_device rd;
    std::mt19937_64 gen( rd() );
    std::uniform_real_distribution<double> dis(min, max);
    for(int ii=0; ii<9; ++ii) mat[ii] = dis(gen); 

    pp[0] = 0; pp[1] = 1; pp[2] = 2;
  }

  void Matrix_double_3by3_Array::gen_hilb()
  {
    for(int ii=0; ii<3; ++ii)
      for(int jj=0; jj<3; ++jj)
        mat[ii*3+jj] = 1.0 / (ii + jj + 1.0);

    pp[0] = 0; pp[1] = 1; pp[2] = 2;
  }

  void Matrix_double_3by3_Array::LU_fac()
  {
    // 1st row
    double max_value = std::abs(mat[0]);
    int max_index = 0;
    bool pivot_flag = false;

    // Pivoting the first column
    if( max_value < std::abs(mat[3]) )
    {
      max_value = std::abs(mat[3]);
      max_index = 1;
      pivot_flag = true;
    }

    if( max_value < std::abs(mat[6]) )
    {
      max_value = std::abs(mat[6]);
      max_index = 2;
      pivot_flag = true;
    }

    if( pivot_flag )
    {
      int int_temp = pp[0];
      pp[0] = pp[max_index];
      pp[max_index] = int_temp;

      double temp = mat[0];
      mat[0] = mat[max_index*3];
      mat[max_index*3] = temp;

      temp = mat[1];
      mat[1] = mat[max_index*3+1];
      mat[max_index*3+1] = temp;

      temp = mat[2];
      mat[2] = mat[max_index*3+2];
      mat[max_index*3+2] = temp;
    }

    mat[3] = mat[3] / mat[0];
    mat[4] = mat[4] - mat[3] * mat[1];
    mat[5] = mat[5] - mat[3] * mat[2];
    mat[6] = mat[6] / mat[0];
    mat[7] = mat[7] - mat[6] * mat[1];
    mat[8] = mat[8] - mat[6] * mat[2];


    // Pivoting the second column
    max_value = std::abs(mat[4]);
    max_index = 1;
    pivot_flag = false;

    if(max_value< std::abs(mat[7]))
    {
      max_value = std::abs(mat[7]);
      max_index = 2;
      pivot_flag = true;
    }

    if(pivot_flag)
    {
      int int_temp = pp[1];
      pp[1] = pp[2];
      pp[2] = int_temp;

      double temp = mat[3];
      mat[3] = mat[6];
      mat[6] = temp;

      temp = mat[4];
      mat[4] = mat[7];
      mat[7] = temp;

      temp = mat[5];
      mat[5] = mat[8];
      mat[8] = temp;
    }


    mat[7] = mat[7] / mat[4];
    mat[8] = mat[8] - mat[7] * mat[5];
  }

  Vector_3 Matrix_double_3by3_Array::LU_solve( const Vector_3 &bb ) const
  {
    Vector_3 xx( bb(pp[0]), bb(pp[1]), bb(pp[2]) );

    xx(1) = xx(1) - mat[3] * xx(0);
    xx(2) = xx(2) - mat[6] * xx(0) - mat[7] * xx(1);

    xx(2) =  xx(2) / mat[8];
    xx(1) = (xx(1) - mat[5] * xx(2)) / mat[4];
    xx(0) = (xx(0) - mat[2] * xx(2)  - mat[1] * xx(1)) / mat[0];

    return xx;
  }

  std::array<double, 3> Matrix_double_3by3_Array::LU_solve( const std::array<double, 3> &bb ) const
  {
    std::array<double, 3> xx {{ bb[pp[0]], bb[pp[1]], bb[pp[2]] }};

    xx[1] = xx[1] - mat[3] * xx[0];
    xx[2] = xx[2] - mat[6] * xx[0] - mat[7] * xx[1];

    xx[2] =  xx[2] / mat[8];
    xx[1] = (xx[1] - mat[5] * xx[2]) / mat[4];
    xx[0] = (xx[0] - mat[2] * xx[2]  - mat[1] * xx[1]) / mat[0];

    return xx;
  }

  void Matrix_double_3by3_Array::LU_solve(const double &b1, const double &b2, const double &b3,
      double &x1, double &x2, double &x3) const
  {
    const double bb[3] = {b1, b2, b3};

    x1 = bb[pp[0]];
    x2 = bb[pp[1]];
    x3 = bb[pp[2]];

    x2 = x2 - mat[3] * x1;
    x3 = x3 - mat[6] * x1 - mat[7] * x2;

    x3 = x3 / mat[8];
    x2 = (x2 - mat[5] * x3) / mat[4];
    x1 = (x1 - mat[2] * x3 - mat[1] * x2) / mat[0];
  }

  void Matrix_double_3by3_Array::transpose()
  {
    double temp = mat[1]; mat[1] = mat[3]; mat[3] = temp;
    temp = mat[2]; mat[2] = mat[6]; mat[6] = temp;
    temp = mat[5]; mat[5] = mat[7]; mat[7] = temp; 
  }

  void Matrix_double_3by3_Array::inverse()
  {
    const double invdetA = 1.0 / det();

    const double temp[9] { invdetA * (mat[4] * mat[8] - mat[5] * mat[7]),
      invdetA * (mat[2] * mat[7] - mat[1] * mat[8]),
      invdetA * (mat[1] * mat[5] - mat[2] * mat[4]),
      invdetA * (mat[5] * mat[6] - mat[3] * mat[8]),
      invdetA * (mat[0] * mat[8] - mat[2] * mat[6]),
      invdetA * (mat[2] * mat[3] - mat[0] * mat[5]),
      invdetA * (mat[3] * mat[7] - mat[4] * mat[6]),
      invdetA * (mat[1] * mat[6] - mat[0] * mat[7]),
      invdetA * (mat[0] * mat[4] - mat[1] * mat[3]) };

    for(int ii=0; ii<9; ++ii) mat[ii] = temp[ii];
  }

  double Matrix_double_3by3_Array::det() const
  {
    return mat[0] * mat[4] * mat[8] + mat[1] * mat[5] * mat[6] 
      + mat[2] * mat[3] * mat[7] - mat[2] * mat[4] * mat[6] 
      - mat[0] * mat[5] * mat[7] - mat[1] * mat[3] * mat[8];
  }

  void Matrix_double_3by3_Array::VecMult(const double * const &xx, double * const &yy) const
  {
    yy[0] = mat[0] * xx[0] + mat[1] * xx[1] + mat[2] * xx[2];
    yy[1] = mat[3] * xx[0] + mat[4] * xx[1] + mat[5] * xx[2];
    yy[2] = mat[6] * xx[0] + mat[7] * xx[1] + mat[8] * xx[2];
  }

  void Matrix_double_3by3_Array::MatMult( const Matrix_double_3by3_Array &mleft,
      const Matrix_double_3by3_Array &mright )
  {
    mat[0] = mleft(0) * mright(0) + mleft(1) * mright(3) + mleft(2) * mright(6);
    mat[1] = mleft(0) * mright(1) + mleft(1) * mright(4) + mleft(2) * mright(7);
    mat[2] = mleft(0) * mright(2) + mleft(1) * mright(5) + mleft(2) * mright(8);

    mat[3] = mleft(3) * mright(0) + mleft(4) * mright(3) + mleft(5) * mright(6);
    mat[4] = mleft(3) * mright(1) + mleft(4) * mright(4) + mleft(5) * mright(7);
    mat[5] = mleft(3) * mright(2) + mleft(4) * mright(5) + mleft(5) * mright(8);

    mat[6] = mleft(6) * mright(0) + mleft(7) * mright(3) + mleft(8) * mright(6);
    mat[7] = mleft(6) * mright(1) + mleft(7) * mright(4) + mleft(8) * mright(7);
    mat[8] = mleft(6) * mright(2) + mleft(7) * mright(5) + mleft(8) * mright(8);
  }

  void Matrix_double_3by3_Array::print() const
  {
    std::cout<<std::setprecision(9)<<mat[0]<<'\t'<<mat[1]<<'\t'<<mat[2]<<std::endl;
    std::cout<<std::setprecision(9)<<mat[3]<<'\t'<<mat[4]<<'\t'<<mat[5]<<std::endl;
    std::cout<<std::setprecision(9)<<mat[6]<<'\t'<<mat[7]<<'\t'<<mat[8]<<std::endl;
  }

  void Matrix_double_3by3_Array::print_full() const
  {
    std::cout<<"Matrix: \n";
    print();
    std::cout<<"pivoting flag: \t";
    std::cout<<std::setprecision(9)<<pp[0]<<'\t'<<pp[1]<<'\t'<<pp[2]<<'\n';
  }

  Matrix_double_6by6_Array::Matrix_double_6by6_Array(const double &aa, const double &bb,
      const double &cc, const double &dd, const double &ee, const double &ff,
      const double &gg, const double &hh, const double &ii)
  {
    Mat[0] = aa*aa;     Mat[1] = dd*dd;     Mat[2] = gg*gg;
    Mat[3] = 2.0*aa*dd; Mat[4] = 2.0*aa*gg; Mat[5] = 2.0*dd*gg;

    Mat[6] = aa*bb;       Mat[7]  = dd*ee;       Mat[8]  = gg*hh;
    Mat[9] = aa*ee+bb*dd; Mat[10] = aa*hh+bb*gg; Mat[11] = dd*hh+ee*gg;

    Mat[12] = aa*cc;       Mat[13] = dd*ff;       Mat[14] = gg*ii;
    Mat[15] = aa*ff+cc*dd; Mat[16] = aa*ii+cc*gg; Mat[17] = dd*ii+ff*gg;

    Mat[18] = bb*bb;     Mat[19] = ee*ee;     Mat[20] = hh*hh;
    Mat[21] = 2.0*bb*ee; Mat[22] = 2.0*bb*hh; Mat[23] = 2.0*ee*hh;

    Mat[24] = bb*cc;       Mat[25] = ee*ff;       Mat[26] = hh*ii;
    Mat[27] = bb*ff+cc*ee; Mat[28] = bb*ii+cc*hh; Mat[29] = ee*ii+ff*hh;

    Mat[30] = cc*cc;     Mat[31] = ff*ff;     Mat[32] = ii*ii;
    Mat[33] = 2.0*cc*ff; Mat[34] = 2.0*cc*ii; Mat[35] = 2.0*ff*ii;

    pp[0] = 0; pp[1] = 1; pp[2] = 2; pp[3] = 3; pp[4] = 4; pp[5] = 5;
  }

  void Matrix_double_6by6_Array::LU_fac()
  {
    for(int kk=0; kk<5; ++kk)
    {
      double max_value = std::abs( Mat[6*kk+kk] );
      int max_index = kk;
      bool pivot_flag = false;
      // find the column pivoting
      for(int ii=kk+1; ii<6; ++ii)
      {
        if( max_value < std::abs( Mat[6*ii+kk] ) )
        {
          max_value = std::abs( Mat[6*ii+kk] );
          max_index = ii;
          pivot_flag = true;
        }
      }

      if(pivot_flag)
      {
        int int_temp = pp[kk];
        pp[kk] = pp[max_index];
        pp[max_index] = int_temp;

        for(int ii=0; ii<6; ++ii)
        {
          double temp = Mat[6*kk+ii];
          Mat[6*kk+ii] = Mat[6*max_index+ii];
          Mat[6*max_index+ii] = temp;
        }
      } 

      for(int ii=kk+1; ii<6; ++ii)
      {
        Mat[6*ii+kk] = Mat[6*ii+kk] / Mat[6*kk+kk];
        for(int jj=kk+1; jj<6; ++jj)
          Mat[6*ii+jj] = Mat[6*ii+jj] - Mat[6*ii+kk] * Mat[6*kk+jj];
      }
    }
  }

  std::array<double, 6> Matrix_double_6by6_Array::LU_solve( const std::array<double, 6> &rhs ) const
  {
    std::array<double, 6> xx {{ rhs[pp[0]], rhs[pp[1]], rhs[pp[2]], rhs[pp[3]], rhs[pp[4]], rhs[pp[5]] }};

    xx[1] = xx[1] - Mat[6]  * xx[0];
    xx[2] = xx[2] - Mat[12] * xx[0] - Mat[13] * xx[1];
    xx[3] = xx[3] - Mat[18] * xx[0] - Mat[19] * xx[1] - Mat[20] * xx[2];
    xx[4] = xx[4] - Mat[24] * xx[0] - Mat[25] * xx[1] - Mat[26] * xx[2] - Mat[27] * xx[3];
    xx[5] = xx[5] - Mat[30] * xx[0] - Mat[31] * xx[1] - Mat[32] * xx[2] - Mat[33] * xx[3] - Mat[34] * xx[4];

    xx[5] =  xx[5] / Mat[35];
    xx[4] = (xx[4] - Mat[29] * xx[5]) / Mat[28];
    xx[3] = (xx[3] - Mat[23] * xx[5]  - Mat[22] * xx[4]) / Mat[21];
    xx[2] = (xx[2] - Mat[17] * xx[5]  - Mat[16] * xx[4]  - Mat[15] * xx[3]) / Mat[14];
    xx[1] = (xx[1] - Mat[11] * xx[5]  - Mat[10] * xx[4]  - Mat[9]  * xx[3]  - Mat[8] * xx[2]) / Mat[7];
    xx[0] = (xx[0] - Mat[5]  * xx[5]  - Mat[4]  * xx[4]  - Mat[3]  * xx[3]  - Mat[2] * xx[2]  - Mat[1] * xx[1]) / Mat[0];

    return xx;
  }

  void Matrix_double_6by6_Array::print() const
  {
    for(int ii=0; ii<6; ++ii)
    {
      for(int jj=0; jj<6; ++jj)
        std::cout<<Mat[6*ii+jj]<<'\t';
      std::cout<<'\n';
    }
    std::cout<<'\n';
    for(int ii=0; ii<6; ++ii)
      std::cout<<pp[ii]<<'\t';

    std::cout<<'\n'<<std::endl;
  }

  QuadPts_on_face::QuadPts_on_face(const FEType &vol_elemType, const int &face_id, 
      const IQuadPts * const lower_quad_rule)
    : dim(lower_quad_rule->get_dim() + 1), num_pts( lower_quad_rule -> get_num_quadPts() )
  {
    if(vol_elemType == FEType::Tet4 || vol_elemType == FEType::Tet10) // Tet element
    {
      //                     t
      //                     ^
      //                     |
      //                     3
      //                    /| `.
      //                   / |    `.
      //                  /  |       `.
      //                 /   |          `.
      //                /    |             `.     
      //               /     |                `.
      //              /      |                   `.   
      //             /      ,0 - - - - - - - - - - -`2 - - -> s
      //            /     ,'  (u)               ,  "
      //           /    ,'                ,  "
      //          /   ,'            ,  "
      //         /  ,'        ,  "
      //        / ,'    ,  "
      //       /,',  "
      //      1'
      //    ,'
      // r *

      SYS_T::print_fatal_if( lower_quad_rule -> get_dim() != 3,
        "Error: FE_T::QuadPts_on_face, wrong surface quadrature rule.\n" );

      qp.assign( 4 * lower_quad_rule->get_num_quadPts(), 0.0 );
      
      switch(face_id)
      {
        case 0: // u = 0 : node1 = node0', node2 = node1', node3 = node2'       //      t                                     s'
          for(int ii {0}; ii < lower_quad_rule->get_num_quadPts(); ++ii)        //      ^                                     ^
          {                                                                     //      3                                     2'
            qp[4*ii + 0] = lower_quad_rule->get_qp(ii, 2);  // r = t'           //      |  `.                     map         |  `.  
            qp[4*ii + 1] = lower_quad_rule->get_qp(ii, 0);  // s = r'           //      |     `.                 <----        |     `.
            qp[4*ii + 2] = lower_quad_rule->get_qp(ii, 1);  // t = s'           //      | front  `.                           |        `.
          }                                                                     //      |           `.                        |           `.
          break;                                                                //  (r) 1 - - - - - - 2 - - -> s         (t') 0'- - - - - - 1'- - -> r'
        
        case 1: // r = 0 : node0 = node0', node3 = node1', node2 = node2'       //      s                                     s'
          for(int ii {0}; ii < lower_quad_rule->get_num_quadPts(); ++ii)        //      ^                                     ^
          {                                                                     //      2                                     2'
            qp[4*ii + 1] = lower_quad_rule->get_qp(ii, 1);  // s = s'           //      |  `.                     map         |  `.
            qp[4*ii + 2] = lower_quad_rule->get_qp(ii, 0);  // t = r'           //      |     `.                 <----        |     `.
            qp[4*ii + 3] = lower_quad_rule->get_qp(ii, 2);  // u = t'           //      | back   `.                           |        `.
          }                                                                     //      |           `.                        |           `.
          break;                                                                //  (u) 0 - - - - - - 3 - - -> t         (t') 0'- - - - - - 1'- - -> r'

        case 2: // s = 0 : node0 = node0', node1 = node1', node3 = node2'       //      t                                     s'
          for(int ii {0}; ii < lower_quad_rule->get_num_quadPts(); ++ii)        //      ^                                     ^
          {                                                                     //      3                                     2'
            qp[4*ii + 0] = lower_quad_rule->get_qp(ii, 0);  // r = r'           //      |  `.                     map         |  `.
            qp[4*ii + 2] = lower_quad_rule->get_qp(ii, 1);  // t = s'           //      |     `.                 <----        |     `.
            qp[4*ii + 3] = lower_quad_rule->get_qp(ii, 2);  // u = t'           //      | left   `.                           |        `.
          }                                                                     //      |           `.                        |           `.
          break;                                                                //  (u) 0 - - - - - - 1 - - -> r         (t') 0'- - - - - - 1'- - -> r'

        case 3: // t = 0 : node0 = node0', node2 = node1', node1 = node2'       //      r                                     s'
          for(int ii {0}; ii < lower_quad_rule->get_num_quadPts(); ++ii)        //      ^                                     ^
          {                                                                     //      1                                     2'
            qp[4*ii + 0] = lower_quad_rule->get_qp(ii, 1);  // r = s'           //      |  `.                     map         |  `.
            qp[4*ii + 1] = lower_quad_rule->get_qp(ii, 0);  // s = r'           //      |     `.                 <----        |     `.
            qp[4*ii + 3] = lower_quad_rule->get_qp(ii, 2);  // u = t'           //      | bottom `.                           |        `.
          }                                                                     //      |           `.                        |           `.
          break;                                                                //  (u) 0 - - - - - - 2 - - -> s         (t') 0'- - - - - - 1'- - -> r'

        default:
          SYS_T::print_fatal("Error: FE_T::QuadPts_on_face, wrong face id input.\n");
          break;
      }
    }
    else if (vol_elemType == FEType::Hex8 || vol_elemType == FEType::Hex27) // Hex element
    {
      //                    t
      //                    ^
      //                    |
      //                    4------------------7
      //                   /.                 /|
      //                  / .                / |
      //                 /  .               /  |
      //                /   .              /   |
      //               /    .             /    |
      //              /     .            /     |
      //             5------------------6      |
      //             |      0...........|......3-------> s
      //             |     .            |     /
      //             |    .             |    /
      //             |   .              |   /
      //             |  .               |  /
      //             | .                | /
      //             |.                 |/
      //             1------------------2
      //            /
      //           *
      //           r
      SYS_T::print_fatal_if( lower_quad_rule -> get_dim() != 2,
        "Error: FE_T::QuadPts_on_face, wrong surface quadrature rule.\n" );

      qp.assign( 3 * lower_quad_rule->get_num_quadPts(), 0.0 );
      
      switch(face_id)
      {
        case 0: // t = 0 : node0 = node0', node3 = node1', node2 = node2', node1 = node3' //    r                          s'
          for(int ii {0}; ii < lower_quad_rule->get_num_quadPts(); ++ii)                  //    ^                          ^
          {                                                                               //    |                          |
            qp[3*ii + 0] = lower_quad_rule->get_qp(ii, 1);  // r = s'                     //    1 -------- 2               3'-------- 2'
            qp[3*ii + 1] = lower_quad_rule->get_qp(ii, 0);  // s = r'                     //    |          |        map    |          |
            qp[3*ii + 2] = 0.0;                             // t = 0                      //    |  bottom  |       <----   |          |
          }                                                                               //    |          |               |          |
          break;                                                                          //    0 -------- 3 - -> s        0'-------- 1'- -> r'

        case 1: // t = 1 : node4 = node0', node5 = node1', node6 = node2', node7 = node3' //    s                          s'
          for(int ii {0}; ii < lower_quad_rule->get_num_quadPts(); ++ii)                  //    ^                          ^
          {                                                                               //    |                          |
            qp[3*ii + 0] = lower_quad_rule->get_qp(ii, 0);  // r = r'                     //    7 -------- 6               3'-------- 2'
            qp[3*ii + 1] = lower_quad_rule->get_qp(ii, 1);  // s = s'                     //    |          |        map    |          |
            qp[3*ii + 2] = 1.0;                             // t = 1                      //    |   top    |       <----   |          |
          }                                                                               //    |          |               |          |
          break;                                                                          //    4 -------- 5 - -> r        0'-------- 1'- -> r'
        
        case 2: // s = 0 : node0 = node0', node1 = node1', node5 = node2', node4 = node3' //    t                          s'
          for(int ii {0}; ii < lower_quad_rule->get_num_quadPts(); ++ii)                  //    ^                          ^
          {                                                                               //    |                          |
            qp[3*ii + 0] = lower_quad_rule->get_qp(ii, 0);  // r = r'                     //    4 -------- 5               3'-------- 2'
            qp[3*ii + 1] = 0.0;                             // s = 0                      //    |          |        map    |          |
            qp[3*ii + 2] = lower_quad_rule->get_qp(ii, 1);  // t = s'                     //    |   left   |       <----   |          |
          }                                                                               //    |          |               |          |
          break;                                                                          //    0 -------- 1 - -> r        0'-------- 1'- -> r'

        case 3: // r = 1 : node1 = node0', node2 = node1', node6 = node2', node5 = node3' //    t                          s'
          for(int ii {0}; ii < lower_quad_rule->get_num_quadPts(); ++ii)                  //    ^                          ^
          {                                                                               //    |                          |
            qp[3*ii + 0] = 1.0;                             // r = 1                      //    5 -------- 6               3'-------- 2'
            qp[3*ii + 1] = lower_quad_rule->get_qp(ii, 0);  // s = r'                     //    |          |        map    |          |
            qp[3*ii + 2] = lower_quad_rule->get_qp(ii, 1);  // t = s'                     //    |   front  |       <----   |          |
          }                                                                               //    |          |               |          |
          break;                                                                          //    1 -------- 2 - -> s        0'-------- 1'- -> r'

        case 4: // s = 1 : node3 = node0', node7 = node1', node6 = node2', node2 = node3' //    r                          s'
          for(int ii {0}; ii < lower_quad_rule->get_num_quadPts(); ++ii)                  //    ^                          ^
          {                                                                               //    |                          |
            qp[3*ii + 0] = lower_quad_rule->get_qp(ii, 1);  // r = s'                     //    2 -------- 6               3'-------- 2'
            qp[3*ii + 1] = 1.0;                             // s = 1                      //    |          |        map    |          |
            qp[3*ii + 2] = lower_quad_rule->get_qp(ii, 0);  // t = r'                     //    |   right  |       <----   |          |
          }                                                                               //    |          |               |          |
          break;                                                                          //    3 -------- 7 - -> t        0'-------- 1'- -> r'
        
        case 5: // r = 0 : node0 = node0', node4 = node1', node7 = node2', node3 = node3' //    s                          s'
          for(int ii {0}; ii < lower_quad_rule->get_num_quadPts(); ++ii)                  //    ^                          ^
          {                                                                               //    |                          |
            qp[3*ii + 0] = 0.0;                             // r = 0                      //    3 -------- 7               3'-------- 2'
            qp[3*ii + 1] = lower_quad_rule->get_qp(ii, 1);  // s = s'                     //    |          |        map    |          |
            qp[3*ii + 2] = lower_quad_rule->get_qp(ii, 0);  // t = r'                     //    |   back   |       <----   |          |
          }                                                                               //    |          |               |          |
          break;                                                                          //    0 -------- 4 - -> t        0'-------- 1'- -> r'
        
        default:
          SYS_T::print_fatal("Error: FE_T::QuadPts_on_face, wrong face id input.\n");
          break;
      }
    }
    else
      SYS_T::print_fatal("Error: FE_T::QuadPts_on_face, unknown element type.\n");
  }

  void QuadPts_on_face::print_info() const
  {
    SYS_T::commPrint("====== Gauss Points on Face =======\n");
    SYS_T::commPrint("Number of points = %d\n", num_pts);
    SYS_T::commPrint("qp.size() = %d\n", qp.size());
    for(int ii=0; ii<num_pts; ++ii)
      SYS_T::commPrint("  %.15f %.15f %.15f\n",
          qp[3*ii], qp[3*ii+1], qp[3*ii+2]);
    SYS_T::commPrint("========================================\n");
  }

}

// EOF
