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

Vector_3 operator+( const Vector_3 &left, const Vector_3 &right )
{
  Vector_3 result;
  result.vec[0] = left(0) + right(0);
  result.vec[1] = left(1) + right(1);
  result.vec[2] = left(2) + right(2);

  return result;
}

Vector_3 operator-( const Vector_3 &left, const Vector_3 &right )
{
  Vector_3 result;
  result.vec[0] = left(0) - right(0);
  result.vec[1] = left(1) - right(1);
  result.vec[2] = left(2) - right(2);

  return result;
}

Vector_3& Vector_3::operator+=( const Vector_3 &source )
{
  vec[0] += source(0);
  vec[1] += source(1);
  vec[2] += source(2);
 
  return *this;
}

Vector_3& Vector_3::operator-=( const Vector_3 &source )
{
  vec[0] -= source(0);
  vec[1] -= source(1);
  vec[2] -= source(2);

  return *this;
}

Vector_3& Vector_3::operator*=( const double &val )
{
  vec[0] *= val;
  vec[1] *= val;
  vec[2] *= val;

  return *this;
}

void Vector_3::print() const
{
  std::cout<<std::setprecision(12)<<vec[0]<<'\t'<<vec[1]<<'\t'<<vec[2]<<std::endl;
}

void Vector_3::gen_zero()
{
  vec[0] = 0.0;
  vec[1] = 0.0;
  vec[2] = 0.0;
}

void Vector_3::gen_val(const double &val)
{
  vec[0] = val;
  vec[1] = val;
  vec[2] = val;
}

void Vector_3::gen_rand()
{
  srand(time(NULL));

  for(int ii=0; ii<3; ++ii)
  {
    const double value = rand() % 10000;

    vec[ii] = value * 1.0e-3 - 5.0;
  }
}

void Vector_3::scale( const double &val )
{
  vec[0] *= val; vec[1] *= val; vec[2] *= val;
}

void Vector_3::AXPY( const double &val, const Vector_3 &source )
{
  vec[0] += val * source(0);
  vec[1] += val * source(1);
  vec[2] += val * source(2);
}

void Vector_3::normalize()
{
  const double len = norm2();
  scale(1.0/len);
}

double Vector_3::dot_product( const Vector_3 &source ) const
{
  return vec[0]*source(0) + vec[1]*source(1) + vec[2]*source(2);
}

double dot_product( const Vector_3 &a, const Vector_3 &b )
{
  return a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
}

Vector_3 cross_product( const Vector_3 &a, const Vector_3 &b )
{
  Vector_3 result;
  result.vec[0] = a(1) * b(2) - a(2) * b(1);
  result.vec[1] = a(2) * b(0) - a(0) * b(2);
  result.vec[2] = a(0) * b(1) - a(1) * b(0);

  return result;
}

int Vector_3::get_dominant_comp() const
{
  int dominant_comp = 0;
  double dominant_val = -INFINITY;

  for(int ii=0; ii<3; ++ii)
  {
    if( std::abs(vec[ii]) > dominant_val )
    {
      dominant_comp = ii;
      dominant_val  = std::abs(vec[ii]);
    }
  }

  return dominant_comp;
}

// EOF
