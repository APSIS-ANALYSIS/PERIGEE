#include "Matrix_3x3.hpp"

Matrix_3x3::Matrix_3x3()
{
  mat[0] = 1.0; mat[1] = 0.0; mat[2] = 0.0;
  mat[3] = 0.0; mat[4] = 1.0; mat[5] = 0.0;
  mat[6] = 0.0; mat[7] = 0.0; mat[8] = 1.0;
}


Matrix_3x3::Matrix_3x3( const Matrix_3x3 &source )
{
  for(int ii=0; ii<9; ++ii) mat[ii] = source(ii);
}

Matrix_3x3::Matrix_3x3( 
    const double &a11, const double &a12, const double &a13,
    const double &a21, const double &a22, const double &a23,
    const double &a31, const double &a32, const double &a33 )
{
  mat[0] = a11; mat[1] = a12; mat[2] = a13;
  mat[3] = a21; mat[4] = a22; mat[5] = a23;
  mat[6] = a31; mat[7] = a32; mat[8] = a33;
}


Matrix_3x3::~Matrix_3x3()
{}


bool Matrix_3x3::is_identical( const Matrix_3x3 source ) const
{
  for(int ii=0; ii<9; ++ii) 
    if(source(ii) != mat[ii]) return false;
  return true;
}


void Matrix_3x3::copy( const Matrix_3x3 &source )
{
  for(int ii=0; ii<9; ++ii) mat[ii] = source(ii);
}


void Matrix_3x3::copy( double source[9] )
{
  for(int ii=0; ii<9; ++ii) mat[ii] = source[ii];
}


Matrix_3x3& Matrix_3x3::operator= (const Matrix_3x3 &source)
{
  if(this == &source) return *this;

  for(int ii=0; ii<9; ++ii) mat[ii] = source(ii);
  return *this; 
}


Matrix_3x3 operator+(const Matrix_3x3 &left, const Matrix_3x3 &right)
{
  Matrix_3x3 result;
  for(int ii=0; ii<9; ++ii) result.mat[ii] = left.mat[ii] + right.mat[ii];
  
  return result;
}


Matrix_3x3 operator-(const Matrix_3x3 &left, const Matrix_3x3 &right)
{
  Matrix_3x3 result;
  for(int ii=0; ii<9; ++ii) result.mat[ii] = left.mat[ii] - right.mat[ii];
  
  return result;
}

Matrix_3x3& Matrix_3x3::operator+= (const Matrix_3x3 &source)
{
  for(int ii=0; ii<9; ++ii) mat[ii] += source(ii);
  return *this;
}


Matrix_3x3& Matrix_3x3::operator-= (const Matrix_3x3 &source)
{
  for(int ii=0; ii<9; ++ii) mat[ii] -= source(ii);
  return *this;
}


Matrix_3x3& Matrix_3x3::operator*= (const double &val)
{
  for(int ii=0; ii<9; ++ii) mat[ii] *= val;
  return *this;
}

void Matrix_3x3::gen_zero()
{
  for(int ii=0; ii<9; ++ii) mat[ii] = 0.0; 
}


void Matrix_3x3::gen_id()
{
  mat[0] = 1.0; mat[1] = 0.0; mat[2] = 0.0;
  mat[3] = 0.0; mat[4] = 1.0; mat[5] = 0.0;
  mat[6] = 0.0; mat[7] = 0.0; mat[8] = 1.0;
}


void Matrix_3x3::gen_rand()
{
  srand(time(NULL));

  for(int ii=0; ii<9; ++ii)
  {
    double value = rand() % 10000;

    mat[ii] = value * 1.0e-3 - 5.0;
  }
}


void Matrix_3x3::gen_hilb()
{
  for(int ii=0; ii<3; ++ii)
    for(int jj=0; jj<3; ++jj)
      mat[ii*3+jj] = 1.0 / (ii + jj + 1.0);
}


void Matrix_3x3::gen_outprod( const double * const &a, const double * const &b )
{
  mat[0] = a[0] * b[0]; mat[1] = a[0] * b[1]; mat[2] = a[0] * b[2];
  mat[3] = a[1] * b[0]; mat[4] = a[1] * b[1]; mat[5] = a[1] * b[2];
  mat[6] = a[2] * b[0]; mat[7] = a[2] * b[1]; mat[8] = a[2] * b[2];
}


void Matrix_3x3::gen_outprod( const double * const &a )
{
  mat[0] = a[0] * a[0]; mat[1] = a[0] * a[1]; mat[2] = a[0] * a[2];
  mat[3] = a[1] * a[0]; mat[4] = a[1] * a[1]; mat[5] = a[1] * a[2];
  mat[6] = a[2] * a[0]; mat[7] = a[2] * a[1]; mat[8] = a[2] * a[2];
}


void Matrix_3x3::transpose()
{
  double temp;
  temp = mat[1]; mat[1] = mat[3]; mat[3] = temp;
  temp = mat[2]; mat[2] = mat[6]; mat[6] = temp;
  temp = mat[5]; mat[5] = mat[7]; mat[7] = temp;
}


void Matrix_3x3::inverse()
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


void Matrix_3x3::scale( const double &val )
{
  for(int ii=0; ii<9; ++ii) mat[ii] = mat[ii] * val;
}


void Matrix_3x3::AXPY( const double &val, const Matrix_3x3 &source )
{
  for(int ii=0; ii<9; ++ii) mat[ii] = mat[ii] + val * source(ii);
}


void Matrix_3x3::PY( const Matrix_3x3 &source )
{
  for(int ii=0; ii<9; ++ii) mat[ii] += source(ii);
}


double Matrix_3x3::det() const
{
  return mat[0] * mat[4] * mat[8] + mat[1] * mat[5] * mat[6]
    + mat[2] * mat[3] * mat[7] - mat[2] * mat[4] * mat[6]
    - mat[0] * mat[5] * mat[7] - mat[1] * mat[3] * mat[8];
}


double Matrix_3x3::I2() const
{
  return 0.5 * ( I1() * I1() - mat[0]*mat[0] - mat[4] * mat[4]
     - mat[8] * mat[8] - 2.0 * ( mat[1]*mat[3] + mat[2] * mat[6] + mat[5] * mat[7] ) );
}

double Matrix_3x3::VecMatVec( const double * const &x,
    const double * const &y ) const
{
  return x[0] * ( mat[0] * y[0] + mat[1] * y[1] + mat[2] * y[2] )
    + x[1] * ( mat[3] * y[0] + mat[4] * y[1] + mat[5] * y[2] )
    + x[2] * ( mat[6] * y[0] + mat[7] * y[1] + mat[8] * y[2] );    
}


double Matrix_3x3::VecMatVec( const Vector_3 &x, const Vector_3 &y ) const
{
  return x(0) * ( mat[0] * y(0) + mat[1] * y(1) + mat[2] * y(2) )
    + x(1) * ( mat[3] * y(0) + mat[4] * y(1) + mat[5] * y(2) )
    + x(2) * ( mat[6] * y(0) + mat[7] * y(1) + mat[8] * y(2) );
}


void Matrix_3x3::VecMult( const double * const &x, double * const &y ) const
{
  y[0] = mat[0] * x[0] + mat[1] * x[1] + mat[2] * x[2];
  y[1] = mat[3] * x[0] + mat[4] * x[1] + mat[5] * x[2];
  y[2] = mat[6] * x[0] + mat[7] * x[1] + mat[8] * x[2];
}


void Matrix_3x3::VecMult( const Vector_3 &x, Vector_3 &y ) const
{
  y(0) = mat[0] * x(0) + mat[1] * x(1) + mat[2] * x(2);
  y(1) = mat[3] * x(0) + mat[4] * x(1) + mat[5] * x(2);
  y(2) = mat[6] * x(0) + mat[7] * x(1) + mat[8] * x(2);
}


void Matrix_3x3::VecMult( const double &x0, const double &x1, const double &x2, 
    double * const &y ) const
{
  y[0] = mat[0] * x0 + mat[1] * x1 + mat[2] * x2;
  y[1] = mat[3] * x0 + mat[4] * x1 + mat[5] * x2;
  y[2] = mat[6] * x0 + mat[7] * x1 + mat[8] * x2;
}


void Matrix_3x3::VecMultT( const double * const &x, double * const &y ) const
{
  y[0] = mat[0] * x[0] + mat[3] * x[1] + mat[6] * x[2];
  y[1] = mat[1] * x[0] + mat[4] * x[1] + mat[7] * x[2];
  y[2] = mat[2] * x[0] + mat[5] * x[1] + mat[8] * x[2];
}


void Matrix_3x3::VecMultT( const Vector_3 &x, Vector_3 &y ) const
{
  y(0) = mat[0] * x(0) + mat[3] * x(1) + mat[6] * x(2);
  y(1) = mat[1] * x(0) + mat[4] * x(1) + mat[7] * x(2);
  y(2) = mat[2] * x(0) + mat[5] * x(1) + mat[8] * x(2);
}


void Matrix_3x3::VecMultT(const double &x0, const double &x1, const double &x2,
    double * const &y ) const 
{
  y[0] = mat[0] * x0 + mat[3] * x1 + mat[6] * x2;
  y[1] = mat[1] * x0 + mat[4] * x1 + mat[7] * x2;
  y[2] = mat[2] * x0 + mat[5] * x1 + mat[8] * x2;
}


void Matrix_3x3::VecMult( double * const &x ) const
{
  double y[3] = {0.0, 0.0, 0.0};
  y[0] = mat[0] * x[0] + mat[1] * x[1] + mat[2] * x[2];
  y[1] = mat[3] * x[0] + mat[4] * x[1] + mat[5] * x[2];
  y[2] = mat[6] * x[0] + mat[7] * x[1] + mat[8] * x[2];
  x[0] = y[0]; x[1] = y[1]; x[2] = y[2];
}


void Matrix_3x3::VecMult( Vector_3 &x ) const
{
  double y[3] = {0.0, 0.0, 0.0};
  y[0] = mat[0] * x(0) + mat[1] * x(1) + mat[2] * x(2);
  y[1] = mat[3] * x(0) + mat[4] * x(1) + mat[5] * x(2);
  y[2] = mat[6] * x(0) + mat[7] * x(1) + mat[8] * x(2);
  x(0) = y[0]; x(1) = y[1]; x(2) = y[2];
}


void Matrix_3x3::MatMult( const Matrix_3x3 &mleft, const Matrix_3x3 &mright )
{
  double temp[9];

  temp[0] = mleft(0) * mright(0) + mleft(1) * mright(3) + mleft(2) * mright(6);
  temp[1] = mleft(0) * mright(1) + mleft(1) * mright(4) + mleft(2) * mright(7);
  temp[2] = mleft(0) * mright(2) + mleft(1) * mright(5) + mleft(2) * mright(8);

  temp[3] = mleft(3) * mright(0) + mleft(4) * mright(3) + mleft(5) * mright(6);
  temp[4] = mleft(3) * mright(1) + mleft(4) * mright(4) + mleft(5) * mright(7);
  temp[5] = mleft(3) * mright(2) + mleft(4) * mright(5) + mleft(5) * mright(8);

  temp[6] = mleft(6) * mright(0) + mleft(7) * mright(3) + mleft(8) * mright(6);
  temp[7] = mleft(6) * mright(1) + mleft(7) * mright(4) + mleft(8) * mright(7);
  temp[8] = mleft(6) * mright(2) + mleft(7) * mright(5) + mleft(8) * mright(8);

  for(int ii=0; ii<9; ++ii) mat[ii] = temp[ii];
}


void Matrix_3x3::MatMultTransposeLeft( const Matrix_3x3 &in )
{
  double temp[9];
  temp[0] = in(0) * in(0) + in(3) * in(3) + in(6) * in(6);
  temp[1] = in(0) * in(1) + in(3) * in(4) + in(6) * in(7);
  temp[2] = in(0) * in(2) + in(3) * in(5) + in(6) * in(8);
  temp[3] = temp[1];
  temp[4] = in(1) * in(1) + in(4) * in(4) + in(7) * in(7);
  temp[5] = in(1) * in(2) + in(4) * in(5) + in(7) * in(8);
  temp[6] = temp[2];
  temp[7] = temp[5];
  temp[8] = in(2) * in(2) + in(5) * in(5) + in(8) * in(8);

  for(int ii=0; ii<9; ++ii) mat[ii] = temp[ii];
}


void Matrix_3x3::MatMultTransposeRight( const Matrix_3x3 &in )
{
  double temp[9];

  temp[0] = in(0) * in(0) + in(1) * in(1) + in(2) * in(2);
  temp[1] = in(0) * in(3) + in(1) * in(4) + in(2) * in(5);
  temp[2] = in(0) * in(6) + in(1) * in(7) + in(2) * in(8);
  temp[3] = temp[1];
  temp[4] = in(3) * in(3) + in(4) * in(4) + in(5) * in(5);
  temp[5] = in(3) * in(6) + in(4) * in(7) + in(5) * in(8);
  temp[6] = temp[2];
  temp[7] = temp[5];
  temp[8] = in(6) * in(6) + in(7) * in(7) + in(8) * in(8);

  for(int ii=0; ii<9; ++ii) mat[ii] = temp[ii];
}


void Matrix_3x3::print() const
{
  std::cout<<std::setprecision(9)<<mat[0]<<'\t'<<mat[1]<<'\t'<<mat[2]<<std::endl;
  std::cout<<std::setprecision(9)<<mat[3]<<'\t'<<mat[4]<<'\t'<<mat[5]<<std::endl;
  std::cout<<std::setprecision(9)<<mat[6]<<'\t'<<mat[7]<<'\t'<<mat[8]<<std::endl;
}


double Matrix_3x3::MatContraction( const Matrix_3x3 &in ) const
{
  return mat[0] * in(0) + mat[1] * in(1) + mat[2] * in(2) + mat[3] * in(3)
    + mat[4] * in(4) + mat[5] * in(5) + mat[6] * in(6) + mat[7] * in(7)
    + mat[8] * in(8);
}


double Matrix_3x3::MatTContraction( const Matrix_3x3 &in ) const
{
  return mat[0] * in(0) + mat[1] * in(3) + mat[2] * in(6) + mat[3] * in(1)
    + mat[4] * in(4) + mat[5] * in(7) + mat[6] * in(2) + mat[7] * in(5)
    + mat[8] * in(8);
}

// EOF
