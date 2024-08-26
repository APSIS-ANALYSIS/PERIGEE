#include "SymmTensor2_3D.hpp"

constexpr std::array<int, 9> SymmTensor2_3D::VoigtMap;

SymmTensor2_3D& SymmTensor2_3D::operator= (const SymmTensor2_3D &source)
{
  if (this != &source) mat = source.mat; // use std::array assignment operator
  
  return *this;
}

SymmTensor2_3D operator+( const SymmTensor2_3D &left, const SymmTensor2_3D &right )
{
  return SymmTensor2_3D( left.mat[0] + right.mat[0], left.mat[1] + right.mat[1], 
		         left.mat[2] + right.mat[2], left.mat[3] + right.mat[3], 
		         left.mat[4] + right.mat[4], left.mat[5] + right.mat[5] );
}

SymmTensor2_3D operator-( const SymmTensor2_3D &left, const SymmTensor2_3D &right )
{
  return SymmTensor2_3D( left.mat[0] - right.mat[0], left.mat[1] - right.mat[1], 
		         left.mat[2] - right.mat[2], left.mat[3] - right.mat[3], 
		         left.mat[4] - right.mat[4], left.mat[5] - right.mat[5] );
}

SymmTensor2_3D& SymmTensor2_3D::operator+=( const SymmTensor2_3D &source )
{
  for(int ii=0; ii<6; ++ii) mat[ii] += source(ii);
  return *this; 
}

 SymmTensor2_3D& SymmTensor2_3D::operator-=( const SymmTensor2_3D &source )
{
  for(int ii=0; ii<6; ++ii) mat[ii] -= source(ii);
  return *this; 
}

SymmTensor2_3D& SymmTensor2_3D::operator*=( const double &val )
{
  for(int ii=0; ii<6; ++ii) mat[ii] *= val;
  return *this;
}

SymmTensor2_3D SymmTensor2_3D::operator- () const
{
  return SymmTensor2_3D( -mat[0], -mat[1], -mat[2], -mat[3], -mat[4], -mat[5] );
}

bool SymmTensor2_3D::is_identical( const SymmTensor2_3D &source, const double &tol ) const
{
  for(int ii=0; ii<6; ++ii)
    if( std::abs( source(ii) - mat[ii]) > tol ) return false;
  return true;  
}

void SymmTensor2_3D::inverse()
{
  const double invdetA = 1.0 / det();

  const double temp[6] = {
    invdetA * (mat[1] * mat[2] - mat[3] * mat[3]),
    invdetA * (mat[0] * mat[2] - mat[4] * mat[4]),
    invdetA * (mat[0] * mat[1] - mat[5] * mat[5]),
    invdetA * (mat[4] * mat[5] - mat[0] * mat[3]),
    invdetA * (mat[5] * mat[3] - mat[4] * mat[1]),
    invdetA * (mat[4] * mat[3] - mat[5] * mat[2]) };

  for(int ii=0; ii<6; ++ii) mat[ii] = temp[ii];
}

void SymmTensor2_3D::AXPY( const double &val, const SymmTensor2_3D &source )
{
  for(int ii=0; ii<6; ++ii) mat[ii] += val * source(ii);
}

void SymmTensor2_3D::AXPI( const double &val )
{
  mat[0] += val; mat[1] += val; mat[2] += val;
}

double SymmTensor2_3D::det() const
{
  return mat[0] * mat[1] * mat[2] + mat[5] * mat[3] * mat[4]
    + mat[4] * mat[5] * mat[3] - mat[4] * mat[1] * mat[4]
    - mat[0] * mat[3] * mat[3] - mat[5] * mat[5] * mat[2];
}

double SymmTensor2_3D::I2() const
{
  return 0.5 * ( I1() * I1() - mat[0] * mat[0] - mat[1] * mat[1] 
     - mat[2] * mat[2] - 2.0 * ( mat[5] * mat[5] + mat[4] * mat[4] + mat[3] * mat[3] ) );
}

double SymmTensor2_3D::VecMatVec( const Vector_3 &x, const Vector_3 &y ) const
{
  return x(0) * ( mat[0] * y(0) + mat[5] * y(1) + mat[4] * y(2) )
    + x(1) * ( mat[5] * y(0) + mat[1] * y(1) + mat[3] * y(2) )
    + x(2) * ( mat[4] * y(0) + mat[3] * y(1) + mat[2] * y(2) );
}

Vector_3 SymmTensor2_3D::VecMult( const Vector_3 &x ) const
{
  return Vector_3( mat[0] * x(0) + mat[5] * x(1) + mat[4] * x(2),
      mat[5] * x(0) + mat[1] * x(1) + mat[3] * x(2),
      mat[4] * x(0) + mat[3] * x(1) + mat[2] * x(2) );
}

void SymmTensor2_3D::VecMult( const double &x0, const double &x1, const double &x2,
       double &y0, double &y1, double &y2 ) const
{
  y0 = mat[0] * x0 + mat[5] * x1 + mat[4] * x2;
  y1 = mat[5] * x0 + mat[1] * x1 + mat[3] * x2;
  y2 = mat[4] * x0 + mat[3] * x1 + mat[2] * x2;
}

void SymmTensor2_3D::MatRot( const Tensor2_3D &Q )
{
  double temp[9] = {0.0};
  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      temp[ii*3+jj] = Q(0,ii) * ( mat[0]*Q(0,jj) + mat[5]*Q(1,jj) + mat[4]*Q(2,jj) )
                    + Q(1,ii) * ( mat[5]*Q(0,jj) + mat[1]*Q(1,jj) + mat[3]*Q(2,jj) )
                    + Q(2,ii) * ( mat[4]*Q(0,jj) + mat[3]*Q(1,jj) + mat[2]*Q(2,jj) );
    }
  }
  
  mat[0] = temp[0]; mat[5] = temp[1]; mat[4] = temp[2];
  mat[1] = temp[4]; mat[3] = temp[5]; mat[2] = temp[8];
}

void SymmTensor2_3D::push_forward_stress( const Tensor2_3D &F )
{
  double temp[9] = {0.0};
  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=ii; jj<3; ++jj)
    {
      temp[ii*3+jj] = F(ii,0) * ( mat[0]*F(jj,0) + mat[5]*F(jj,1) + mat[4]*F(jj,2) )
                    + F(ii,1) * ( mat[5]*F(jj,0) + mat[1]*F(jj,1) + mat[3]*F(jj,2) )
                    + F(ii,2) * ( mat[4]*F(jj,0) + mat[3]*F(jj,1) + mat[2]*F(jj,2) );
    }
  }
  
  mat[0] = temp[0]; mat[5] = temp[1]; mat[4] = temp[2];
  mat[1] = temp[4]; mat[3] = temp[5]; mat[2] = temp[8];
}

double SymmTensor2_3D::MatContraction( const Tensor2_3D &source ) const
{
  return mat[0] * source(0) + mat[5] * source(1) + mat[4] * source(2) + mat[5] * source(3)
    + mat[1] * source(4) + mat[3] * source(5) + mat[4] * source(6) + mat[3] * source(7)
    + mat[2] * source(8);
}

double SymmTensor2_3D::MatContraction( const SymmTensor2_3D &source ) const
{
  return mat[0] * source(0) + mat[1] * source(1) + mat[2] * source(2)
    + 2.0 * ( mat[3] * source(3) + mat[4] * source(4) + mat[5] * source(5) );
}

void SymmTensor2_3D::print() const
{
  std::cout<<std::setprecision(9)<<mat[0]<<'\t'<<mat[5]<<'\t'<<mat[4]<<std::endl;
  std::cout<<std::setprecision(9)<<mat[5]<<'\t'<<mat[1]<<'\t'<<mat[3]<<std::endl;
  std::cout<<std::setprecision(9)<<mat[4]<<'\t'<<mat[3]<<'\t'<<mat[2]<<std::endl;
}

void SymmTensor2_3D::print_in_row() const
{
  std::cout<<std::setprecision(9)<<mat[0]<<'\t'<<mat[5]<<'\t'<<mat[4]<<'\t';
  std::cout<<std::setprecision(9)<<mat[5]<<'\t'<<mat[1]<<'\t'<<mat[3]<<'\t';
  std::cout<<std::setprecision(9)<<mat[4]<<'\t'<<mat[3]<<'\t'<<mat[2]<<std::endl;
}

void SymmTensor2_3D::print_Voigt() const
{
  std::cout<<std::setprecision(9)<<mat[0]<<'\t'<<mat[1]<<'\t'<<mat[2]<<'\t';
  std::cout<<std::setprecision(9)<<mat[3]<<'\t'<<mat[4]<<'\t'<<mat[5]<<std::endl;
}

int SymmTensor2_3D::eigen_decomp( double &eta1, double &eta2, double &eta3,
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
        v2.normalize();               // w1 
        v2 = Vec3::cross_product( v1, v2 ); // v2 = w1 x v1
        v3 = Vec3::cross_product( v1, v2 ); // v3 = v1 x v2
      }
      else
      {
        v3.normalize();               // w1 
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

void SymmTensor2_3D::find_eigen_vector( const double &eta, Vector_3 &v,
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

    b -= Vec3::dot_product(s1, b) * s1;

    c -= Vec3::dot_product(s1, c) * s1;

    if( b.norm2() >= c.norm2() )
    {
      b.normalize(); s2 = b; v = Vec3::cross_product(s1,s2);
    }
    else
    {
      c.normalize(); s2 = c; v = Vec3::cross_product(s1,s2);
    }
  }
  else if( len_b >= len_a && len_b >= len_c )
  {
    b.normalize(); // b is s1 now

    s1 = b;

    a -= Vec3::dot_product(s1, a) * s1;

    c -= Vec3::dot_product(s1, c) * s1;

    if( a.norm2() >= c.norm2() )
    {
      a.normalize(); s2 = a; v = Vec3::cross_product(s1, s2);
    }
    else
    {
      c.normalize(); s2 = c; v = Vec3::cross_product(s1, s2);
    }
  }
  else
  {
    c.normalize(); // c is s1 now

    s1 = c;

    a -= Vec3::dot_product(s1, a) * s1;

    b -= Vec3::dot_product(s1, b) * s1;

    if(a.norm2() >= b.norm2())
    {
      a.normalize(); s2 = a; v = Vec3::cross_product(s1, s2);
    }
    else
    {
      b.normalize(); s2 = b; v = Vec3::cross_product(s1, s2);
    }
  }
}

double SymmTensor2_3D::J2() const
{
  const double a = mat[0] * mat[0] + 2.0 * mat[5] * mat[5] + 2.0 * mat[4] * mat[4]
    + mat[1] * mat[1] + 2.0 * mat[3] * mat[3] + mat[2] * mat[2];

  const double b = mat[0] + mat[1] + mat[2];

  return 0.5 * a - b * b / 6.0;
}

double SymmTensor2_3D::J3() const
{
  const double a = ( mat[0] + mat[1] + mat[2] ) / 3.0;

  const double m0 = mat[0] - a;
  const double m4 = mat[1] - a;
  const double m8 = mat[2] - a;

  return m0 * m4 * m8 + mat[5] * mat[3] * mat[4]
    + mat[4] * mat[5] * mat[3] - mat[4] * m4 * mat[4]
    - m0 * mat[3] * mat[3] - mat[5] * mat[5] * m8;
}

Vector_3 operator*( const SymmTensor2_3D &left, const Vector_3 &right )
{
  return Vector_3( left.xx() * right.x() + left.xy() * right.y() + left.xz() * right.z(),
      left.yx() * right.x() + left.yy() * right.y() + left.yz() * right.z(),
      left.zx() * right.x() + left.zy() * right.y() + left.zz() * right.z() ); 
}

Tensor2_3D operator*( const SymmTensor2_3D &left, const Tensor2_3D &right )
{
  return Tensor2_3D( left(0) * right(0) + left(5) * right(3) + left(4) * right(6),
   left(0) * right(1) + left(5) * right(4) + left(4) * right(7),
   left(0) * right(2) + left(5) * right(5) + left(4) * right(8),
   left(5) * right(0) + left(1) * right(3) + left(3) * right(6),
   left(5) * right(1) + left(1) * right(4) + left(3) * right(7),
   left(5) * right(2) + left(1) * right(5) + left(3) * right(8),
   left(4) * right(0) + left(3) * right(3) + left(2) * right(6),
   left(4) * right(1) + left(3) * right(4) + left(2) * right(7),
   left(4) * right(2) + left(3) * right(5) + left(2) * right(8) );
}

Tensor2_3D operator*( const Tensor2_3D &left, const SymmTensor2_3D &right )
{
  return Tensor2_3D( left(0) * right(0) + left(1) * right(5) + left(2) * right(4),
   left(0) * right(5) + left(1) * right(1) + left(2) * right(3),
   left(0) * right(4) + left(1) * right(3) + left(2) * right(2),
   left(3) * right(0) + left(4) * right(5) + left(5) * right(4),
   left(3) * right(5) + left(4) * right(1) + left(5) * right(3),
   left(3) * right(4) + left(4) * right(3) + left(5) * right(2),
   left(6) * right(0) + left(7) * right(5) + left(8) * right(4),
   left(6) * right(5) + left(7) * right(1) + left(8) * right(3),
   left(6) * right(4) + left(7) * right(3) + left(8) * right(2) );
}

Tensor2_3D operator*( const SymmTensor2_3D &left, const SymmTensor2_3D &right )
{
  return Tensor2_3D( left(0) * right(0) + left(5) * right(5) + left(4) * right(4),
   left(0) * right(5) + left(5) * right(1) + left(4) * right(3),
   left(0) * right(4) + left(5) * right(3) + left(4) * right(2),
   left(5) * right(0) + left(1) * right(5) + left(3) * right(4),
   left(5) * right(5) + left(1) * right(1) + left(3) * right(3),
   left(5) * right(4) + left(1) * right(3) + left(3) * right(2),
   left(4) * right(0) + left(3) * right(5) + left(2) * right(4),
   left(4) * right(5) + left(3) * right(1) + left(2) * right(3),
   left(4) * right(4) + left(3) * right(3) + left(2) * right(2) );
}

SymmTensor2_3D operator*( const double &val, const SymmTensor2_3D &input )
{
  return SymmTensor2_3D( val * input(0), val * input(1), val * input(2), 
   val * input(3), val * input(4), val * input(5) );
}

SymmTensor2_3D STen2::gen_id()
{
  return SymmTensor2_3D(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
}

SymmTensor2_3D STen2::gen_zero()
{
  return SymmTensor2_3D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}

SymmTensor2_3D STen2::gen_rand(const double &left, const double &right)
{
  std::random_device rd;
  std::mt19937_64 gen( rd() );
  std::uniform_real_distribution<double> dis(left, right);
  return SymmTensor2_3D( dis(gen), dis(gen), dis(gen), dis(gen), dis(gen), dis(gen) );
}

SymmTensor2_3D STen2::gen_dyad( const Vector_3 &input )
{
  return SymmTensor2_3D( input(0) * input(0), input(1) * input(1),
      input(2) * input(2), input(1) * input(2), input(0) * input(2),
      input(0) * input(1) );
}

SymmTensor2_3D STen2::inverse( const SymmTensor2_3D &input )
{
  const double invdetA = 1.0 / input.det();

  return SymmTensor2_3D (
  invdetA * (input(1) * input(2) - input(3) * input(3)),
  invdetA * (input(0) * input(2) - input(4) * input(4)),
  invdetA * (input(0) * input(1) - input(5) * input(5)),
  invdetA * (input(4) * input(5) - input(0) * input(3)),
  invdetA * (input(5) * input(3) - input(4) * input(1)),
  invdetA * (input(4) * input(3) - input(5) * input(2)) );
}

SymmTensor2_3D STen2::gen_right_Cauchy_Green( const Tensor2_3D &input )
{
  return SymmTensor2_3D (
   input(0) * input(0) + input(3) * input(3) + input(6) * input(6),
   input(1) * input(1) + input(4) * input(4) + input(7) * input(7),
   input(2) * input(2) + input(5) * input(5) + input(8) * input(8),
   input(1) * input(2) + input(4) * input(5) + input(7) * input(8),
   input(0) * input(2) + input(3) * input(5) + input(6) * input(8),
   input(0) * input(1) + input(3) * input(4) + input(6) * input(7) );
}

SymmTensor2_3D STen2::gen_left_Cauchy_Green( const Tensor2_3D &input )
{
  return SymmTensor2_3D (
  input(0) * input(0) + input(1) * input(1) + input(2) * input(2),
  input(3) * input(3) + input(4) * input(4) + input(5) * input(5),
  input(6) * input(6) + input(7) * input(7) + input(8) * input(8),
  input(3) * input(6) + input(4) * input(7) + input(5) * input(8),
  input(0) * input(6) + input(1) * input(7) + input(2) * input(8),
  input(0) * input(3) + input(1) * input(4) + input(2) * input(5) );
}

SymmTensor2_3D STen2::gen_symm_part( const Tensor2_3D &input )
{
  return SymmTensor2_3D( input(0), input(4), input(8),
                        0.5 * ( input(5) + input(7) ),
                        0.5 * ( input(2) + input(6) ),
                        0.5 * ( input(1) + input(3) ) );
}

SymmTensor2_3D STen2::gen_DEV_part( const SymmTensor2_3D &input, const SymmTensor2_3D &CC )
{
  const double val = input.MatContraction(CC) / 3.0;
  const auto invC  = inverse(CC);
  return SymmTensor2_3D( input(0) - val * invC(0), 
                         input(1) - val * invC(1), 
                         input(2) - val * invC(2), 
                         input(3) - val * invC(3), 
                         input(4) - val * invC(4), 
                         input(5) - val * invC(5) );
}

// EOF
