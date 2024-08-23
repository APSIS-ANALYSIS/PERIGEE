#include "Tensor2_3D.hpp"

bool Tensor2_3D::is_identical( const Tensor2_3D &source, const double &tol ) const
{
  for(int ii=0; ii<9; ++ii) 
    if( std::abs( source(ii) - mat[ii]) > tol ) return false;
  return true;
}

Tensor2_3D& Tensor2_3D::operator= (const Tensor2_3D &source)
{
  if(this != &source) mat = source.mat; // use std::array assignment operator

  return *this; 
}

Tensor2_3D operator+(const Tensor2_3D &left, const Tensor2_3D &right)
{
  Tensor2_3D result;
  for(int ii=0; ii<9; ++ii) result.mat[ii] = left.mat[ii] + right.mat[ii];
  
  return result;
}

Tensor2_3D operator-(const Tensor2_3D &left, const Tensor2_3D &right)
{
  Tensor2_3D result;
  for(int ii=0; ii<9; ++ii) result.mat[ii] = left.mat[ii] - right.mat[ii];
  
  return result;
}

Tensor2_3D& Tensor2_3D::operator+= (const Tensor2_3D &source)
{
  for(int ii=0; ii<9; ++ii) mat[ii] += source(ii);
  return *this;
}

Tensor2_3D& Tensor2_3D::operator-= (const Tensor2_3D &source)
{
  for(int ii=0; ii<9; ++ii) mat[ii] -= source(ii);
  return *this;
}

Tensor2_3D& Tensor2_3D::operator*= (const double &val)
{
  for(int ii=0; ii<9; ++ii) mat[ii] *= val;
  return *this;
}

Tensor2_3D Tensor2_3D::operator- () const
{
  return Tensor2_3D( -mat[0], -mat[1], -mat[2], -mat[3], -mat[4],
      -mat[5], -mat[6], -mat[7], -mat[8] );
}

void Tensor2_3D::gen_rand(const double &left, const double &right)
{
  std::random_device rd;
  std::mt19937_64 gen( rd() );
  std::uniform_real_distribution<double> dis(left, right);
  for(int ii=0; ii<9; ++ii) mat[ii] = dis(gen);
}

void Tensor2_3D::gen_hilb()
{
  for(int ii=0; ii<3; ++ii)
    for(int jj=0; jj<3; ++jj)
      mat[ii*3+jj] = 1.0 / (ii + jj + 1.0);
}

void Tensor2_3D::gen_outprod( const Vector_3 &va, const Vector_3 &vb ) 
{
  mat[0] = va(0) * vb(0); mat[1] = va(0) * vb(1); mat[2] = va(0) * vb(2);
  mat[3] = va(1) * vb(0); mat[4] = va(1) * vb(1); mat[5] = va(1) * vb(2);
  mat[6] = va(2) * vb(0); mat[7] = va(2) * vb(1); mat[8] = va(2) * vb(2);
}

void Tensor2_3D::gen_outprod( const Vector_3 &va )
{
  mat[0] = va(0) * va(0); mat[1] = va(0) * va(1); mat[2] = va(0) * va(2);
  mat[3] = va(1) * va(0); mat[4] = va(1) * va(1); mat[5] = va(1) * va(2);
  mat[6] = va(2) * va(0); mat[7] = va(2) * va(1); mat[8] = va(2) * va(2);
}

void Tensor2_3D::add_outprod( const double &val, const Vector_3 &va, const Vector_3 &vb ) 
{
  mat[0] += val * va(0) * vb(0); mat[1] += val * va(0) * vb(1); mat[2] += val * va(0) * vb(2);
  mat[3] += val * va(1) * vb(0); mat[4] += val * va(1) * vb(1); mat[5] += val * va(1) * vb(2);
  mat[6] += val * va(2) * vb(0); mat[7] += val * va(2) * vb(1); mat[8] += val * va(2) * vb(2);
}

void Tensor2_3D::transpose()
{
  double temp; // temperary variable for swapping off diagonal entries
  temp = mat[1]; mat[1] = mat[3]; mat[3] = temp;
  temp = mat[2]; mat[2] = mat[6]; mat[6] = temp;
  temp = mat[5]; mat[5] = mat[7]; mat[7] = temp;
}

void Tensor2_3D::inverse()
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

void Tensor2_3D::scale( const double &val )
{
  for(int ii=0; ii<9; ++ii) mat[ii] *= val;
}

void Tensor2_3D::AXPY( const double &val, const Tensor2_3D &source )
{
  for(int ii=0; ii<9; ++ii) mat[ii] += val * source(ii);
}

void Tensor2_3D::AXPI( const double &val )
{
  mat[0] += val; mat[4] += val; mat[8] += val;
}

double Tensor2_3D::det() const
{
  return mat[0] * mat[4] * mat[8] + mat[1] * mat[5] * mat[6]
    + mat[2] * mat[3] * mat[7] - mat[2] * mat[4] * mat[6]
    - mat[0] * mat[5] * mat[7] - mat[1] * mat[3] * mat[8];
}

double Tensor2_3D::I2() const
{
  return 0.5 * ( I1() * I1() - mat[0]*mat[0] - mat[4] * mat[4]
     - mat[8] * mat[8] - 2.0 * ( mat[1]*mat[3] + mat[2] * mat[6] + mat[5] * mat[7] ) );
}

double Tensor2_3D::VecMatVec( const Vector_3 &x, const Vector_3 &y ) const
{
  return x(0) * ( mat[0] * y(0) + mat[1] * y(1) + mat[2] * y(2) )
    + x(1) * ( mat[3] * y(0) + mat[4] * y(1) + mat[5] * y(2) )
    + x(2) * ( mat[6] * y(0) + mat[7] * y(1) + mat[8] * y(2) );
}

Vector_3 Tensor2_3D::VecMult( const Vector_3 &x ) const
{
  return Vector_3( mat[0] * x(0) + mat[1] * x(1) + mat[2] * x(2), 
      mat[3] * x(0) + mat[4] * x(1) + mat[5] * x(2), 
      mat[6] * x(0) + mat[7] * x(1) + mat[8] * x(2) );
}

void Tensor2_3D::VecMult( const double &x0, const double &x1, const double &x2, 
    double &y0, double &y1, double &y2 ) const
{
  y0 = mat[0] * x0 + mat[1] * x1 + mat[2] * x2;
  y1 = mat[3] * x0 + mat[4] * x1 + mat[5] * x2;
  y2 = mat[6] * x0 + mat[7] * x1 + mat[8] * x2;
}

Vector_3 Tensor2_3D::VecMultT( const Vector_3 &x ) const
{
  return Vector_3( mat[0] * x(0) + mat[3] * x(1) + mat[6] * x(2),
      mat[1] * x(0) + mat[4] * x(1) + mat[7] * x(2),
      mat[2] * x(0) + mat[5] * x(1) + mat[8] * x(2) );
}

void Tensor2_3D::VecMultT(const double &x0, const double &x1, const double &x2,
    double &y0, double &y1, double &y2 ) const 
{
  y0 = mat[0] * x0 + mat[3] * x1 + mat[6] * x2;
  y1 = mat[1] * x0 + mat[4] * x1 + mat[7] * x2;
  y2 = mat[2] * x0 + mat[5] * x1 + mat[8] * x2;
}

void Tensor2_3D::MatMult( const Tensor2_3D &mleft, const Tensor2_3D &mright )
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

void Tensor2_3D::MatRot( const Tensor2_3D &Q )
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

void Tensor2_3D::MatMultTransposeLeft( const Tensor2_3D &in )
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

void Tensor2_3D::MatMultTransposeRight( const Tensor2_3D &in )
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

void Tensor2_3D::print(std::ostream& os, const std::string& delimiter) const
{
  os << std::setprecision(12);
  os<<mat[0]<<delimiter<<mat[1]<<delimiter<<mat[2]<<'\n';
  os<<mat[3]<<delimiter<<mat[4]<<delimiter<<mat[5]<<'\n';
  os<<mat[6]<<delimiter<<mat[7]<<delimiter<<mat[8]<<std::endl;
}

void Tensor2_3D::print_in_row() const
{
  std::cout<<std::setprecision(9)<<mat[0]<<'\t'<<mat[1]<<'\t'<<mat[2]<<'\t';
  std::cout<<std::setprecision(9)<<mat[3]<<'\t'<<mat[4]<<'\t'<<mat[5]<<'\t';
  std::cout<<std::setprecision(9)<<mat[6]<<'\t'<<mat[7]<<'\t'<<mat[8]<<std::endl;
}

void Tensor2_3D::print_Voigt() const
{
  std::cout<<std::setprecision(9)<<mat[0]<<'\t'<<mat[4]<<'\t'<<mat[8]<<'\t';
  std::cout<<std::setprecision(9)<<mat[5]<<'\t'<<mat[2]<<'\t'<<mat[1]<<std::endl;
}

double Tensor2_3D::MatContraction( const Tensor2_3D &in ) const
{
  return mat[0] * in(0) + mat[1] * in(1) + mat[2] * in(2) + mat[3] * in(3)
    + mat[4] * in(4) + mat[5] * in(5) + mat[6] * in(6) + mat[7] * in(7)
    + mat[8] * in(8);
}

double Tensor2_3D::MatTContraction( const Tensor2_3D &in ) const
{
  return mat[0] * in(0) + mat[1] * in(3) + mat[2] * in(6) + mat[3] * in(1)
    + mat[4] * in(4) + mat[5] * in(7) + mat[6] * in(2) + mat[7] * in(5)
    + mat[8] * in(8);
}

void Tensor2_3D::find_eigen_vector( const double &eta, Vector_3 &v,
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
    s1 = Vec3::normalize(a); // s1 = a/|a|

    b -= Vec3::dot_product(s1, b) * s1;

    c -= Vec3::dot_product(s1, c) * s1;

    if( b.norm2() >= c.norm2() )
    {
      s2 = Vec3::normalize(b); v = Vec3::cross_product(s1, s2);
    }
    else
    {
      s2 = Vec3::normalize(c); v = Vec3::cross_product(s1, s2);
    }
  }
  else if( len_b >= len_a && len_b >= len_c )
  {
    s1 = Vec3::normalize(b); // s1 = b/|b|

    a -= Vec3::dot_product(s1, a) * s1;

    c -= Vec3::dot_product(s1, c) * s1;

    if( a.norm2() >= c.norm2() )
    {
      s2 = Vec3::normalize(a); v = Vec3::cross_product(s1, s2);
    }
    else
    {
      s2 = Vec3::normalize(c); v = Vec3::cross_product(s1, s2);
    }
  }
  else
  {
    s1 = Vec3::normalize(c); // s1 = c/|c|

    a -= Vec3::dot_product(s1, a) * s1;

    b -= Vec3::dot_product(s1, b) * s1;

    if(a.norm2() >= b.norm2())
    {
      s2 = Vec3::normalize(a); v = Vec3::cross_product(s1, s2);
    }
    else
    {
      s2 = Vec3::normalize(b); v = Vec3::cross_product(s1, s2);
    }
  }
}

double Tensor2_3D::J2() const
{
  const double a = mat[0] * mat[0] + 2.0 * mat[1] * mat[3] + 2.0 * mat[2] * mat[6]
    + mat[4] * mat[4] + 2.0 * mat[5] * mat[7] + mat[8] * mat[8];

  const double b = mat[0] + mat[4] + mat[8];

  return 0.5 * a - b * b / 6.0;
}

double Tensor2_3D::J3() const
{
  const double a = ( mat[0] + mat[4] + mat[8] ) / 3.0;

  const double m0 = mat[0] - a;
  const double m4 = mat[4] - a;
  const double m8 = mat[8] - a;

  return m0 * m4 * m8 + mat[1] * mat[5] * mat[6]
    + mat[2] * mat[3] * mat[7] - mat[2] * m4 * mat[6]
    - m0 * mat[5] * mat[7] - mat[1] * mat[3] * m8;
}

int Tensor2_3D::eigen_decomp( double &eta1, double &eta2, double &eta3,
    Vector_3 &v1, Vector_3 &v2, Vector_3 &v3 ) const
{
  constexpr double frac13 = 1.0 / 3.0;
  const double frac13_tr = tr() * frac13; // value used to shift the eigenvalue

  const double mJ2 = J2();
  
  const double val = 2.0 * std::pow( mJ2*frac13, 0.5 );

  const double PI = atan(1.0) * 4.0;

  // If J2 is zero, we directly get three identical eignevalues,
  // meaning the matrix is eta I.
  if( mJ2 <= 1.0e-14 )
  {
    eta1 = frac13_tr; eta2 = eta1; eta3 = eta1;
    v1 = Vec3::gen_e1(); v2 = Vec3::gen_e2(); v3 = Vec3::gen_e3();
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
    const double A22 = VecMatVec(v2,v2) - frac13_tr * Vec3::dot_product(v2, v2);
    const double A23 = VecMatVec(v2,v3) - frac13_tr * Vec3::dot_product(v2, v3);
    const double A33 = VecMatVec(v3,v3) - frac13_tr * Vec3::dot_product(v3, v3);

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
      v2 -= (frac13_tr + eta2) * temp;

      temp = v3;
      v3 = (*this) * v3; 
      v3 -= (frac13_tr + eta2) * temp;

      if( v2.norm2() >= v3.norm2() )
      {
        v2.normalize();                     // w1 
        v2 = Vec3::cross_product( v1, v2 ); // v2 = w1 x v1
        v3 = Vec3::cross_product( v1, v2 ); // v3 = v1 x v2
      }
      else
      {
        v3.normalize();                     // w1 
        v2 = Vec3::cross_product( v1, v3 ); // v2 = w1 x v1
        v3 = Vec3::cross_product( v1, v2 ); // v3 = v1 x v2
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

Vector_3 operator*(const Tensor2_3D &left, const Vector_3 &right)
{
  return Vector_3( left.xx() * right.x() + left.xy() * right.y() + left.xz() * right.z(), 
      left.yx() * right.x() + left.yy() * right.y() + left.yz() * right.z(),
      left.zx() * right.x() + left.zy() * right.y() + left.zz() * right.z() );
}

Tensor2_3D operator*(const Tensor2_3D &mleft, const Tensor2_3D &mright)
{
  return Tensor2_3D( mleft(0) * mright(0) + mleft(1) * mright(3) + mleft(2) * mright(6),
   mleft(0) * mright(1) + mleft(1) * mright(4) + mleft(2) * mright(7),
   mleft(0) * mright(2) + mleft(1) * mright(5) + mleft(2) * mright(8),
   mleft(3) * mright(0) + mleft(4) * mright(3) + mleft(5) * mright(6),
   mleft(3) * mright(1) + mleft(4) * mright(4) + mleft(5) * mright(7),
   mleft(3) * mright(2) + mleft(4) * mright(5) + mleft(5) * mright(8),
   mleft(6) * mright(0) + mleft(7) * mright(3) + mleft(8) * mright(6),
   mleft(6) * mright(1) + mleft(7) * mright(4) + mleft(8) * mright(7),
   mleft(6) * mright(2) + mleft(7) * mright(5) + mleft(8) * mright(8) );
}

Tensor2_3D operator*( const double &val, const Tensor2_3D &input )
{
  return Tensor2_3D( val * input(0), val * input(1), val * input(2),
      val * input(3), val * input(4), val * input(5),
      val * input(6), val * input(7), val * input(8) );
}

Tensor2_3D Ten2::inverse( const Tensor2_3D &input )
{
  const double invdet = 1.0 / input.det();

  return Tensor2_3D( invdet * (input(4) * input(8) - input(5) * input(7)),
    invdet * (input(2) * input(7) - input(1) * input(8)),
    invdet * (input(1) * input(5) - input(2) * input(4)),
    invdet * (input(5) * input(6) - input(3) * input(8)),
    invdet * (input(0) * input(8) - input(2) * input(6)),
    invdet * (input(2) * input(3) - input(0) * input(5)),
    invdet * (input(3) * input(7) - input(4) * input(6)),
    invdet * (input(1) * input(6) - input(0) * input(7)),
    invdet * (input(0) * input(4) - input(1) * input(3)) );
}

Tensor2_3D Ten2::cofactor( const Tensor2_3D &input )
{
  return Tensor2_3D( (input(4) * input(8) - input(5) * input(7)),
      (input(5) * input(6) - input(3) * input(8)),
      (input(3) * input(7) - input(4) * input(6)),
      (input(2) * input(7) - input(1) * input(8)),
      (input(0) * input(8) - input(2) * input(6)),
      (input(1) * input(6) - input(0) * input(7)),
      (input(1) * input(5) - input(2) * input(4)),
      (input(2) * input(3) - input(0) * input(5)),
      (input(0) * input(4) - input(1) * input(3)) );
}

Tensor2_3D Ten2::exp( const Tensor2_3D &input )
{
  double nn = 0.0;
  double nn_fac = 1.0;

  Tensor2_3D input_pow = Ten2::gen_id();
  Tensor2_3D input_exp = Ten2::gen_id();

  do
  {
    nn += 1.0;
    nn_fac *= nn;

    input_pow = input_pow * input;
    input_exp += ( 1.0/nn_fac ) * input_pow;

  }while( std::sqrt( input_pow.MatContraction(input_pow) ) / nn_fac >= 1.0e-16 );

  return input_exp;
}

// EOF
