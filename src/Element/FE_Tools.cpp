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
  double a [16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double b [4]  = {0.0, 0.0, 0.0, 0.0};

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

  MATH_T::Matrix_dense AA(4);

  AA.set_values( a );

  AA.LU_fac();

  double out[4] = {0.0, 0.0, 0.0, 0.0};

  AA.LU_solve(b, out);

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

namespace FE_T
{
  Matrix_double_3by3_Array::Matrix_double_3by3_Array()
  {
    mat[0] = 1.0; mat[1] = 0.0; mat[2] = 0.0;
    mat[3] = 0.0; mat[4] = 1.0; mat[5] = 0.0;
    mat[6] = 0.0; mat[7] = 0.0; mat[8] = 1.0;

    p[0] = 0; p[1] = 1; p[2] = 2;

    invm0 = 1.0; invm1 = 1.0; invm2 = 1.0;
  }

  Matrix_double_3by3_Array::Matrix_double_3by3_Array( 
      const double &a11, const double &a12, const double &a13,
      const double &a21, const double &a22, const double &a23,
      const double &a31, const double &a32, const double &a33 )
  {
    mat[0] = a11;  mat[1] = a12;  mat[2] = a13;
    mat[3] = a21;  mat[4] = a22;  mat[5] = a23;
    mat[6] = a31;  mat[7] = a32;  mat[8] = a33;

    p[0] = 0; p[1] = 1; p[2] = 2;

    invm0 = 1.0; invm1 = 1.0; invm2 = 1.0;
  }

  Matrix_double_3by3_Array::Matrix_double_3by3_Array( 
      const Matrix_double_3by3_Array &other )
  {
    for(int ii=0; ii<9; ++ii) mat[ii] = other.mat[ii];

    p[0] = other.p[0]; 
    p[1] = other.p[1]; 
    p[2] = other.p[2];

    invm0 = other.invm0;
    invm1 = other.invm1;
    invm2 = other.invm2; 
  }

  Matrix_double_3by3_Array::~Matrix_double_3by3_Array()
  {}

  Matrix_double_3by3_Array& Matrix_double_3by3_Array::operator= (
      const Matrix_double_3by3_Array &input )
  {
    if(this != &input)
    {
      for(int ii=0; ii<9; ++ii) mat[ii] = input.mat[ii];
      p[0] = input.p[0];
      p[1] = input.p[1];
      p[2] = input.p[2];
      invm0 = input.invm0;
      invm1 = input.invm1;
      invm2 = input.invm2;
    }
    return *this;
  }

  void Matrix_double_3by3_Array::gen_id()
  {
    mat[0] = 1.0;  mat[1] = 0.0;  mat[2] = 0.0;
    mat[3] = 0.0;  mat[4] = 1.0;  mat[5] = 0.0;
    mat[6] = 0.0;  mat[7] = 0.0;  mat[8] = 1.0;

    p[0] = 0; p[1] = 1; p[2] = 2;

    invm0 = 1.0; invm1 = 1.0; invm2 = 1.0;
  }

  void Matrix_double_3by3_Array::gen_rand(const double &min, const double &max)
  {
    std::random_device rd;
    std::mt19937_64 gen( rd() );
    std::uniform_real_distribution<double> dis(min, max);
    for(int ii=0; ii<9; ++ii) mat[ii] = dis(gen); 

    p[0] = 0; p[1] = 1; p[2] = 2;

    invm0 = 1.0; invm1 = 1.0; invm2 = 1.0;
  }

  void Matrix_double_3by3_Array::gen_hilb()
  {
    for(int ii=0; ii<3; ++ii)
      for(int jj=0; jj<3; ++jj)
        mat[ii*3+jj] = 1.0 / (ii + jj + 1.0);

    p[0] = 0; p[1] = 1; p[2] = 2;

    invm0 = 1.0; invm1 = 1.0; invm2 = 1.0;
  }

  void Matrix_double_3by3_Array::LU_fac()
  {
    double temp;
    int int_temp;
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
      int_temp = p[0];
      p[0] = p[max_index];
      p[max_index] = int_temp;

      temp = mat[0];
      mat[0] = mat[max_index*3];
      mat[max_index*3] = temp;

      temp = mat[1];
      mat[1] = mat[max_index*3+1];
      mat[max_index*3+1] = temp;

      temp = mat[2];
      mat[2] = mat[max_index*3+2];
      mat[max_index*3+2] = temp;
    }

    double invA0 = 1.0 / mat[0];
    mat[3] = mat[3] * invA0;
    mat[4] = mat[4] - mat[3] * mat[1];
    mat[5] = mat[5] - mat[3] * mat[2];
    mat[6] = mat[6] * invA0;
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
      int_temp = p[1];
      p[1] = p[2];
      p[2] = int_temp;

      temp = mat[3];
      mat[3] = mat[6];
      mat[6] = temp;

      temp = mat[4];
      mat[4] = mat[7];
      mat[7] = temp;

      temp = mat[5];
      mat[5] = mat[8];
      mat[8] = temp;
    }

    double invA1 = 1.0 / mat[4];
    mat[7] = mat[7] * invA1;
    mat[8] = mat[8] - mat[7] * mat[5];

    invm0 = 1.0 / mat[0]; invm1 = 1.0 / mat[4]; invm2 = 1.0 / mat[8];
  }

  Vector_3 Matrix_double_3by3_Array::LU_solve( const Vector_3 &b ) const
  {
    Vector_3 x( b(p[0]), b(p[1]), b(p[2]) );

    x(1) = x(1) - mat[3] * x(0);
    x(2) = x(2) - mat[6] * x(0) - mat[7] * x(1);

    x(2) = x(2) * invm2;
    x(1) = (x(1) - mat[5] * x(2)) * invm1;
    x(0) = (x(0) - mat[2] * x(2) - mat[1] * x(1)) * invm0;

    return x;
  }

  std::array<double, 3> Matrix_double_3by3_Array::LU_solve( const std::array<double, 3> &b ) const
  {
    std::array<double, 3> x {{ b[p[0]], b[p[1]], b[p[2]] }};

    x[1] = x[1] - mat[3] * x[0];
    x[2] = x[2] - mat[6] * x[0] - mat[7] * x[1];

    x[2] = x[2] * invm2;
    x[1] = (x[1] - mat[5] * x[2]) * invm1;
    x[0] = (x[0] - mat[2] * x[2] - mat[1] * x[1]) * invm0;

    return x;
  }

  void Matrix_double_3by3_Array::LU_solve(const double &b1, const double &b2, const double &b3,
      double &x1, double &x2, double &x3) const
  {
    const double b[3] = {b1, b2, b3};

    x1 = b[p[0]];
    x2 = b[p[1]];
    x3 = b[p[2]];

    x2 = x2 - mat[3] * x1;
    x3 = x3 - mat[6] * x1 - mat[7] * x2;

    x3 = x3 * invm2;
    x2 = (x2 - mat[5] * x3) * invm1;
    x1 = (x1 - mat[2] * x3 - mat[1] * x2) * invm0;
  }

  void Matrix_double_3by3_Array::transpose()
  {
    double temp;
    temp = mat[1]; mat[1] = mat[3]; mat[3] = temp;
    temp = mat[2]; mat[2] = mat[6]; mat[6] = temp;
    temp = mat[5]; mat[5] = mat[7]; mat[7] = temp; 
  }

  void Matrix_double_3by3_Array::inverse()
  {
    const double invdetA = 1.0 / det();

    double temp[9];

    temp[0] = invdetA * (mat[4] * mat[8] - mat[5] * mat[7]);
    temp[1] = invdetA * (mat[2] * mat[7] - mat[1] * mat[8]);
    temp[2] = invdetA * (mat[1] * mat[5] - mat[2] * mat[4]);
    temp[3] = invdetA * (mat[5] * mat[6] - mat[3] * mat[8]);
    temp[4] = invdetA * (mat[0] * mat[8] - mat[2] * mat[6]);
    temp[5] = invdetA * (mat[2] * mat[3] - mat[0] * mat[5]);
    temp[6] = invdetA * (mat[3] * mat[7] - mat[4] * mat[6]);
    temp[7] = invdetA * (mat[1] * mat[6] - mat[0] * mat[7]);
    temp[8] = invdetA * (mat[0] * mat[4] - mat[1] * mat[3]);

    for(int ii=0; ii<9; ++ii) mat[ii] = temp[ii];
  }

  double Matrix_double_3by3_Array::det() const
  {
    return mat[0] * mat[4] * mat[8] + mat[1] * mat[5] * mat[6] 
      + mat[2] * mat[3] * mat[7] - mat[2] * mat[4] * mat[6] 
      - mat[0] * mat[5] * mat[7] - mat[1] * mat[3] * mat[8];
  }

  void Matrix_double_3by3_Array::VecMult(const double * const &x, double * const &y) const
  {
    y[0] = mat[0] * x[0] + mat[1] * x[1] + mat[2] * x[2];
    y[1] = mat[3] * x[0] + mat[4] * x[1] + mat[5] * x[2];
    y[2] = mat[6] * x[0] + mat[7] * x[1] + mat[8] * x[2];
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
    std::cout<<std::setprecision(9)<<p[0]<<'\t'<<p[1]<<'\t'<<p[2]<<'\n';
    std::cout<<"invm0 = "<<invm0<<'\t';
    std::cout<<"invm1 = "<<invm1<<'\t';
    std::cout<<"invm2 = "<<invm2<<'\n';
  }

}

// EOF