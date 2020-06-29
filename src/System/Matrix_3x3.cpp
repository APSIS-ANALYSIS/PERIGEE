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


void Matrix_3x3::gen_outprod( const Vector_3 &a, const Vector_3 &b ) 
{
  mat[0] = a(0) * b(0); mat[1] = a(0) * b(1); mat[2] = a(0) * b(2);
  mat[3] = a(1) * b(0); mat[4] = a(1) * b(1); mat[5] = a(1) * b(2);
  mat[6] = a(2) * b(0); mat[7] = a(2) * b(1); mat[8] = a(2) * b(2);
}


void Matrix_3x3::gen_outprod( const double * const &a )
{
  mat[0] = a[0] * a[0]; mat[1] = a[0] * a[1]; mat[2] = a[0] * a[2];
  mat[3] = a[1] * a[0]; mat[4] = a[1] * a[1]; mat[5] = a[1] * a[2];
  mat[6] = a[2] * a[0]; mat[7] = a[2] * a[1]; mat[8] = a[2] * a[2];
}


void Matrix_3x3::gen_outprod( const Vector_3 &a )
{
  mat[0] = a(0) * a(0); mat[1] = a(0) * a(1); mat[2] = a(0) * a(2);
  mat[3] = a(1) * a(0); mat[4] = a(1) * a(1); mat[5] = a(1) * a(2);
  mat[6] = a(2) * a(0); mat[7] = a(2) * a(1); mat[8] = a(2) * a(2);
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


void Matrix_3x3::find_eigen_vector( const double &eta, Vector_3 &v,
    Vector_3 &s1, Vector_3 &s2 ) const
{
  const double frac13_tr = tr() / 3.0; // value used to shift the eigenvalue
  
  Vector_3 a, b, c;

  a.gen_e1(); b.gen_e2(); c.gen_e3();

  VecMult(a); a(0) -= ( eta + frac13_tr );
  VecMult(b); b(1) -= ( eta + frac13_tr );
  VecMult(c); c(2) -= ( eta + frac13_tr );

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
    else eta1 = val * cos(alpha + 4.0*PI*frac13);

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

    if( eta3 != eta2)
    {
      // Now form u1 and u2
      // u1 = (A - 0.333 tr - eta2 ) s1
      // u2 = (A - 0.333 tr - eta2 ) s2
      Vector_3 temp; temp.copy(v2);
      VecMult(v2);
      v2.AXPY( -1.0*(frac13_tr + eta2), temp );

      temp.copy(v3);
      VecMult(v3); v3.AXPY( -1.0*(frac13_tr + eta2), temp );

      if( v2.norm2() >= v3.norm2() )
      {
        v2.normalize(); // w1 
        v2 = cross_product( v1, v2 ); // v2 = w1 x v1
        v3 = cross_product( v1, v2 ); // v3 = v1 x v2
      }
      else
      {
        v3.normalize(); // w1 
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

// EOF
