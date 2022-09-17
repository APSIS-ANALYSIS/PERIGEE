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



// EOF
