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

Matrix_3x3::Matrix_3x3(
    const Vector_3 &vec1, const Vector_3 &vec2, const Vector_3 &vec3 )
{
  mat[0] = vec1(0); mat[1] = vec2(0); mat[2] = vec3(0);
  mat[3] = vec1(1); mat[4] = vec2(1); mat[5] = vec3(1);
  mat[6] = vec1(2); mat[7] = vec2(2); mat[8] = vec3(2);
}

Matrix_3x3::~Matrix_3x3()
{}

bool Matrix_3x3::is_identical( const Matrix_3x3 &source, const double &tol ) const
{
  for(int ii=0; ii<9; ++ii) 
    if( std::abs( source(ii) - mat[ii]) > tol ) return false;
  return true;
}

void Matrix_3x3::copy( const Matrix_3x3 &source )
{
  for(int ii=0; ii<9; ++ii) mat[ii] = source(ii);
}

void Matrix_3x3::copy( const double source[9] )
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

Matrix_3x3 Matrix_3x3::operator- () const
{
  return Matrix_3x3( -mat[0], -mat[1], -mat[2], -mat[3], -mat[4],
      -mat[5], -mat[6], -mat[7], -mat[8] );
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

void Matrix_3x3::gen_outprod( const Vector_3 &va, const Vector_3 &vb ) 
{
  mat[0] = va(0) * vb(0); mat[1] = va(0) * vb(1); mat[2] = va(0) * vb(2);
  mat[3] = va(1) * vb(0); mat[4] = va(1) * vb(1); mat[5] = va(1) * vb(2);
  mat[6] = va(2) * vb(0); mat[7] = va(2) * vb(1); mat[8] = va(2) * vb(2);
}

void Matrix_3x3::gen_outprod( const Vector_3 &va )
{
  mat[0] = va(0) * va(0); mat[1] = va(0) * va(1); mat[2] = va(0) * va(2);
  mat[3] = va(1) * va(0); mat[4] = va(1) * va(1); mat[5] = va(1) * va(2);
  mat[6] = va(2) * va(0); mat[7] = va(2) * va(1); mat[8] = va(2) * va(2);
}

void Matrix_3x3::add_outprod( const double &val, const Vector_3 &va, const Vector_3 &vb ) 
{
  mat[0] += val * va(0) * vb(0); mat[1] += val * va(0) * vb(1); mat[2] += val * va(0) * vb(2);
  mat[3] += val * va(1) * vb(0); mat[4] += val * va(1) * vb(1); mat[5] += val * va(1) * vb(2);
  mat[6] += val * va(2) * vb(0); mat[7] += val * va(2) * vb(1); mat[8] += val * va(2) * vb(2);
}

void Matrix_3x3::transpose()
{
  double temp; // temperary variable for swapping off diagonal entries
  temp = mat[1]; mat[1] = mat[3]; mat[3] = temp;
  temp = mat[2]; mat[2] = mat[6]; mat[6] = temp;
  temp = mat[5]; mat[5] = mat[7]; mat[7] = temp;
}

void Matrix_3x3::inverse()
{
  const double invdetA = 1.0 / det();

  const double temp[9] = { 
    invdetA * (mat[4] * mat[8] - mat[5] * mat[7]),
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

void Matrix_3x3::scale( const double &val )
{
  for(int ii=0; ii<9; ++ii) mat[ii] *= val;
}

void Matrix_3x3::AXPY( const double &val, const Matrix_3x3 &source )
{
  for(int ii=0; ii<9; ++ii) mat[ii] += val * source(ii);
}

void Matrix_3x3::AXPI( const double &val )
{
  mat[0] += val; mat[4] += val; mat[8] += val;
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

double Matrix_3x3::VecMatVec( const Vector_3 &x, const Vector_3 &y ) const
{
  return x(0) * ( mat[0] * y(0) + mat[1] * y(1) + mat[2] * y(2) )
    + x(1) * ( mat[3] * y(0) + mat[4] * y(1) + mat[5] * y(2) )
    + x(2) * ( mat[6] * y(0) + mat[7] * y(1) + mat[8] * y(2) );
}

Vector_3 Matrix_3x3::VecMult( const Vector_3 &x ) const
{
  return Vector_3( mat[0] * x(0) + mat[1] * x(1) + mat[2] * x(2), 
      mat[3] * x(0) + mat[4] * x(1) + mat[5] * x(2), 
      mat[6] * x(0) + mat[7] * x(1) + mat[8] * x(2) );
}

void Matrix_3x3::VecMult( const double &x0, const double &x1, const double &x2, 
    double &y0, double &y1, double &y2 ) const
{
  y0 = mat[0] * x0 + mat[1] * x1 + mat[2] * x2;
  y1 = mat[3] * x0 + mat[4] * x1 + mat[5] * x2;
  y2 = mat[6] * x0 + mat[7] * x1 + mat[8] * x2;
}

Vector_3 Matrix_3x3::VecMultT( const Vector_3 &x ) const
{
  return Vector_3( mat[0] * x(0) + mat[3] * x(1) + mat[6] * x(2),
      mat[1] * x(0) + mat[4] * x(1) + mat[7] * x(2),
      mat[2] * x(0) + mat[5] * x(1) + mat[8] * x(2) );
}

void Matrix_3x3::VecMultT(const double &x0, const double &x1, const double &x2,
    double &y0, double &y1, double &y2 ) const 
{
  y0 = mat[0] * x0 + mat[3] * x1 + mat[6] * x2;
  y1 = mat[1] * x0 + mat[4] * x1 + mat[7] * x2;
  y2 = mat[2] * x0 + mat[5] * x1 + mat[8] * x2;
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

void Matrix_3x3::MatRot( const Matrix_3x3 &Q )
{
  double temp[9] = {0.0};
  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      temp[ii*3+jj] = Q(0,ii) * ( mat[0]*Q(0,jj) + mat[1]*Q(1,jj) + mat[2]*Q(2,jj) )
                    + Q(1,ii) * ( mat[3]*Q(0,jj) + mat[4]*Q(1,jj) + mat[5]*Q(2,jj) )
                    + Q(2,ii) * ( mat[6]*Q(0,jj) + mat[7]*Q(1,jj) + mat[8]*Q(2,jj) );
    }
  }
  
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

void Matrix_3x3::print_in_row() const
{
  std::cout<<std::setprecision(9)<<mat[0]<<'\t'<<mat[1]<<'\t'<<mat[2]<<'\t';
  std::cout<<std::setprecision(9)<<mat[3]<<'\t'<<mat[4]<<'\t'<<mat[5]<<'\t';
  std::cout<<std::setprecision(9)<<mat[6]<<'\t'<<mat[7]<<'\t'<<mat[8]<<std::endl;
}

void Matrix_3x3::print_Voigt() const
{
  std::cout<<std::setprecision(9)<<mat[0]<<'\t'<<mat[4]<<'\t'<<mat[8]<<'\t';
  std::cout<<std::setprecision(9)<<mat[5]<<'\t'<<mat[2]<<'\t'<<mat[1]<<std::endl;
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

void Matrix_3x3::find_eigen_vector( const double &eta, Vector_3 &v,
    Vector_3 &s1, Vector_3 &s2 ) const
{
  const double frac13_tr = tr() / 3.0; // value used to shift the eigenvalue
  
  Vector_3 a = (*this) * Vector_3(1.0, 0.0, 0.0); 
  Vector_3 b = (*this) * Vector_3(0.0, 1.0, 0.0); 
  Vector_3 c = (*this) * Vector_3(0.0, 0.0, 1.0); 
  
  a(0) -= ( eta + frac13_tr );
  b(1) -= ( eta + frac13_tr );
  c(2) -= ( eta + frac13_tr );

  const double len_a = a.norm2();
  const double len_b = b.norm2();
  const double len_c = c.norm2();

  if( len_a >= len_b && len_a >= len_c )
  {
    a.normalize(); // a is s1 now

    s1 = a;

    double val = dot_product(s1, b);
    b.AXPY(-1.0*val, s1);

    val = dot_product(s1, c);
    c.AXPY(-1.0*val, s1);

    if( b.norm2() >= c.norm2() )
    {
      b.normalize(); s2 = b; v = cross_product(s1,s2);
    }
    else
    {
      c.normalize(); s2 = c; v = cross_product(s1,s2);
    }
  }
  else if( len_b >= len_a && len_b >= len_c )
  {
    b.normalize(); // b is s1 now

    s1 = b;

    double val = dot_product(s1, a);
    a.AXPY(-1.0*val, s1);

    val = dot_product(s1, c);
    c.AXPY(-1.0*val, s1);

    if( a.norm2() >= c.norm2() )
    {
      a.normalize(); s2 = a; v = cross_product(s1, s2);
    }
    else
    {
      c.normalize(); s2 = c; v = cross_product(s1, s2);
    }
  }
  else
  {
    c.normalize(); // c is s1 now

    s1 = c;

    double val = dot_product(s1, a);
    a.AXPY(-1.0*val, s1);

    val = dot_product(s1, b);
    b.AXPY(-1.0*val, s1);

    if(a.norm2() >= b.norm2())
    {
      a.normalize(); s2 = a; v = cross_product(s1, s2);
    }
    else
    {
      b.normalize(); s2 = b; v = cross_product(s1, s2);
    }
  }
}

double Matrix_3x3::J2() const
{
  const double a = mat[0] * mat[0] + 2.0 * mat[1] * mat[3] + 2.0 * mat[2] * mat[6]
    + mat[4] * mat[4] + 2.0 * mat[5] * mat[7] + mat[8] * mat[8];

  const double b = mat[0] + mat[4] + mat[8];

  return 0.5 * a - b * b / 6.0;
}

double Matrix_3x3::J3() const
{
  const double a = ( mat[0] + mat[4] + mat[8] ) / 3.0;

  const double m0 = mat[0] - a;
  const double m4 = mat[4] - a;
  const double m8 = mat[8] - a;

  return m0 * m4 * m8 + mat[1] * mat[5] * mat[6]
    + mat[2] * mat[3] * mat[7] - mat[2] * m4 * mat[6]
    - m0 * mat[5] * mat[7] - mat[1] * mat[3] * m8;
}

int Matrix_3x3::eigen_decomp( double &eta1, double &eta2, double &eta3,
           Vector_3 &v1, Vector_3 &v2, Vector_3 &v3 ) const
{
  const double frac13 = 1.0 / 3.0;
  const double frac13_tr = tr() * frac13; // value used to shift the eigenvalue

  const double mJ2 = J2();
  
  const double val = 2.0 * std::pow( mJ2*frac13, 0.5 );

  const double PI = atan(1.0) * 4.0;

  // If J2 is zero, we directly get three identical eignevalues,
  // meaning the matrix is eta I.
  if( mJ2 <= 1.0e-14 )
  {
    eta1 = frac13_tr; eta2 = eta1; eta3 = eta1;
    v1.gen_e1(); v2.gen_e2(); v3.gen_e3();
    return 1;
  }
  else
  {
    const double mJ3 = J3();

    double val_cos = 0.5*mJ3*std::pow(3.0/mJ2, 1.5);

    if( val_cos > 1.0 ) val_cos = 1.0;
  
    if( val_cos < -1.0 ) val_cos = -1.0;
    
    // Find the angle pincipal value
    const double alpha = acos( val_cos ) * frac13;

    // Find the most distinct eigenvalue as eta1
    if( alpha <= PI*frac13*0.5 ) eta1 = val * cos(alpha); 
    else eta1 = val * cos(alpha + 2.0*PI*frac13);

    // v1 is determined, v2 and v3 are used to hold s1 s2 for now
    find_eigen_vector(eta1, v1, v2, v3);

    // Form the reduced matrix
    const double A22 = VecMatVec(v2,v2) - frac13_tr * dot_product(v2, v2);
    const double A23 = VecMatVec(v2,v3) - frac13_tr * dot_product(v2, v3);
    const double A33 = VecMatVec(v3,v3) - frac13_tr * dot_product(v3, v3);

    const double diff_2233 = A22 - A33;
    if( diff_2233 >= 0.0 )
      eta2 = 0.5 * ( (A22+A33) - std::sqrt( diff_2233 * diff_2233 + 4.0*A23*A23 ) );
    else
      eta2 = 0.5 * ( (A22+A33) + std::sqrt( diff_2233 * diff_2233 + 4.0*A23*A23 ) );

    eta3 = A22 + A33 - eta2;

    if( std::abs( eta3 - eta2 ) > 1.0e-10 )
    {
      // Now form u1 and u2
      // u1 = (A - 0.333 tr - eta2 ) s1
      // u2 = (A - 0.333 tr - eta2 ) s2
      Vector_3 temp(v2);
      v2 = (*this) * v2;
      v2.AXPY( -1.0*(frac13_tr + eta2), temp );

      temp.copy(v3);
      v3 = (*this) * v3; 
      v3.AXPY( -1.0*(frac13_tr + eta2), temp );

      if( v2.norm2() >= v3.norm2() )
      {
        v2.normalize();               // w1 
        v2 = cross_product( v1, v2 ); // v2 = w1 x v1
        v3 = cross_product( v1, v2 ); // v3 = v1 x v2
      }
      else
      {
        v3.normalize();               // w1 
        v2 = cross_product( v1, v3 ); // v2 = w1 x v1
        v3 = cross_product( v1, v2 ); // v3 = v1 x v2
      }

      // Shift back from deviatoric to the original matrix eigenvalues
      eta1 += frac13_tr; eta2 += frac13_tr; eta3 += frac13_tr;
      return 3;
    }
    else
    {
      // Shift back from deviatoric to the original matrix eigenvalues
      eta1 += frac13_tr; eta2 += frac13_tr; eta3 += frac13_tr;
      
      return 2;
    }
  }
}

Vector_3 operator*(const Matrix_3x3 &left, const Vector_3 &right)
{
  return Vector_3( left.xx() * right.x() + left.xy() * right.y() + left.xz() * right.z(), 
      left.yx() * right.x() + left.yy() * right.y() + left.yz() * right.z(),
      left.zx() * right.x() + left.zy() * right.y() + left.zz() * right.z() );
}

Matrix_3x3 operator*(const Matrix_3x3 &mleft, const Matrix_3x3 &mright)
{
  return Matrix_3x3( mleft(0) * mright(0) + mleft(1) * mright(3) + mleft(2) * mright(6),
   mleft(0) * mright(1) + mleft(1) * mright(4) + mleft(2) * mright(7),
   mleft(0) * mright(2) + mleft(1) * mright(5) + mleft(2) * mright(8),
   mleft(3) * mright(0) + mleft(4) * mright(3) + mleft(5) * mright(6),
   mleft(3) * mright(1) + mleft(4) * mright(4) + mleft(5) * mright(7),
   mleft(3) * mright(2) + mleft(4) * mright(5) + mleft(5) * mright(8),
   mleft(6) * mright(0) + mleft(7) * mright(3) + mleft(8) * mright(6),
   mleft(6) * mright(1) + mleft(7) * mright(4) + mleft(8) * mright(7),
   mleft(6) * mright(2) + mleft(7) * mright(5) + mleft(8) * mright(8) );
}

Matrix_3x3 operator*( const double &val, const Matrix_3x3 &input )
{
  return Matrix_3x3( val * input(0), val * input(1), val * input(2),
      val * input(3), val * input(4), val * input(5),
      val * input(6), val * input(7), val * input(8) );
}

Matrix_3x3 inverse( const Matrix_3x3 &input )
{
  const double invdet = 1.0 / input.det();

  return Matrix_3x3( invdet * (input(4) * input(8) - input(5) * input(7)),
    invdet * (input(2) * input(7) - input(1) * input(8)),
    invdet * (input(1) * input(5) - input(2) * input(4)),
    invdet * (input(5) * input(6) - input(3) * input(8)),
    invdet * (input(0) * input(8) - input(2) * input(6)),
    invdet * (input(2) * input(3) - input(0) * input(5)),
    invdet * (input(3) * input(7) - input(4) * input(6)),
    invdet * (input(1) * input(6) - input(0) * input(7)),
    invdet * (input(0) * input(4) - input(1) * input(3)) );
}

Matrix_3x3 cofactor( const Matrix_3x3 &input )
{
  return Matrix_3x3( (input(4) * input(8) - input(5) * input(7)),
      (input(5) * input(6) - input(3) * input(8)),
      (input(3) * input(7) - input(4) * input(6)),
      (input(2) * input(7) - input(1) * input(8)),
      (input(0) * input(8) - input(2) * input(6)),
      (input(1) * input(6) - input(0) * input(7)),
      (input(1) * input(5) - input(2) * input(4)),
      (input(2) * input(3) - input(0) * input(5)),
      (input(0) * input(4) - input(1) * input(3)) );
}

Matrix_3x3 transpose( const Matrix_3x3 &input )
{
  return Matrix_3x3( input(0), input(3), input(6),
      input(1), input(4), input(7),
      input(2), input(5), input(8) );
}

Matrix_3x3 gen_identity_matrix()
{
  return Matrix_3x3( 1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0 );
}

// EOF
