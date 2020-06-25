#include "Vector_3.hpp"

Vector_3::Vector_3()
{
  vec[0] = 0.0; vec[1] = 0.0; vec[2] = 0.0;
}


Vector_3::Vector_3( const Vector_3 &source )
{
  vec[0] = source(0); vec[1] = source(1); vec[2] = source(2);
}


Vector_3::Vector_3( const double &v0, const double &v1, const double &v2 )
{
  vec[0] = v0; vec[1] = v1; vec[2] = v2;
}


Vector_3::~Vector_3()
{}


void Vector_3::copy( const Vector_3 &source )
{
  vec[0] = source(0); vec[1] = source(1); vec[2] = source(2);
}


void Vector_3::copy( double source[3] )
{
  vec[0] = source[0]; vec[1] = source[1]; vec[2] = source[2];
}

Vector_3& Vector_3::operator= (const Vector_3 &source)
{
  if(this == &source) return *this;

  vec[0] = source(0); vec[1] = source(1); vec[2] = source(2);

  return *this;
}


void Vector_3::print() const
{
  std::cout<<std::setprecision(9)<<vec[0]<<'\t'<<vec[1]<<'\t'<<vec[2]<<std::endl;
}



// EOF
