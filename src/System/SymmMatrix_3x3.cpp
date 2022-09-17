#include "SymmMatrix_3x3.hpp"

SymmMatrix_3x3::SymmMatrix_3x3()
{
  mat[0] = 1.0; mat[1] = 1.0; mat[2] = 1.0;
  mat[3] = 0.0; mat[4] = 0.0; mat[5] = 0.0;
}

SymmMatrix_3x3::SymmMatrix_3x3( const SymmMatrix_3x3 &source )
{
  mat[0] = source(0); mat[1] = source(1); mat[2] = source(2);
  mat[3] = source(3); mat[4] = source(4); mat[5] = source(5);
}

SymmMatrix_3x3::SymmMatrix_3x3( const double &m0, const double &m1, 
    const double &m2, const double &m3, const double &m4, const double &m5 )
{
  mat[0] = m0; mat[1] = m1; mat[2] = m2;
  mat[3] = m3; mat[4] = m4; mat[5] = m5;
}

SymmMatrix_3x3::~SymmMatrix_3x3()
{}

SymmMatrix_3x3 operator+( const SymmMatrix_3x3 &left, const SymmMatrix_3x3 &right )
{
  SymmMatrix_3x3 result;
  for(int ii=0; ii<6; ii++) result.mat[ii] = left.mat[ii] + right.mat[ii];

  return result;
}

SymmMatrix_3x3 operator-( const SymmMatrix_3x3 &left, const SymmMatrix_3x3 &right )
{
  SymmMatrix_3x3 result;
  for(int ii=0; ii<6; ii++) result.mat[ii] = left.mat[ii] - right.mat[ii];

  return result;
}

SymmMatrix_3x3& SymmMatrix_3x3::operator+=( const SymmMatrix_3x3 &source )
{
  for(int ii=0; ii<6; ii++) mat[ii] += source(ii);
  return *this; 
}

 SymmMatrix_3x3& SymmMatrix_3x3::operator-=( const SymmMatrix_3x3 &source )
{
  for(int ii=0; ii<6; ii++) mat[ii] -= source(ii);
  return *this; 
}

SymmMatrix_3x3& SymmMatrix_3x3::operator*=( const double &val )
{
  for(int ii=0; ii<6; ii++) mat[ii] *= val;
  return *this;
}

bool SymmMatrix_3x3::is_identical( const SymmMatrix_3x3 source ) const
{
  for(int ii=0; ii<6; ii++)
    if(source(ii) != mat[ii]) return false;
  return true;  
}

void SymmMatrix_3x3::gen_zero()
{
  for(int ii=0; ii<6; ii++) mat[ii] = 0.0;
}

void SymmMatrix_3x3::gen_id()
{
  mat[0] = 1.0; mat[1] = 1.0; mat[2] = 1.0;
  mat[3] = 0.0; mat[4] = 0.0; mat[5] = 0.0; 
}

void SymmMatrix_3x3::gen_rand()
{
  srand(time(NULL));

  for(int ii=0; ii<6; ii++)
  {
    double value = rand() % 10000;

    mat[ii] = value * 1.0e-3 - 5.0;
  }
}

void SymmMatrix_3x3::inverse()
{
  const double invdetA = 1.0 / det();

  double temp[6];

  temp[0] = invdetA * (mat[1] * mat[2] - mat[3] * mat[3]);
  temp[5] = invdetA * (mat[4] * mat[3] - mat[5] * mat[2]);
  temp[4] = invdetA * (mat[5] * mat[3] - mat[4] * mat[1]);
  temp[1] = invdetA * (mat[0] * mat[2] - mat[4] * mat[4]);
  temp[3] = invdetA * (mat[4] * mat[5] - mat[0] * mat[3]);
  temp[2] = invdetA * (mat[0] * mat[1] - mat[5] * mat[5]);

  for(int ii=0; ii<6; ii++) mat[ii] = temp[ii];
}

void SymmMatrix_3x3::scale( const double &val )
{
  for(int ii=0; ii<6; ii++) mat[ii] = mat[ii] * val;
}

void SymmMatrix_3x3::AXPY( const double &val, const SymmMatrix_3x3 &source )
{
  for(int ii=0; ii<6; ii++) mat[ii] = mat[ii] + val * source(ii);
}

void SymmMatrix_3x3::AXPI( const double &val )
{
  mat[0] += val; mat[1] += val; mat[2] += val;
}

void SymmMatrix_3x3::PY( const SymmMatrix_3x3 &source )
{
  for(int ii=0; ii<6; ii++) mat[ii] += source(ii);
}

// EOF
