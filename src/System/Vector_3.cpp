#include "Vector_3.hpp"

Vector_3& Vector_3::operator= (const Vector_3 &source)
{
  if(this != &source) vec = source.vec;

  return *this;
}

Vector_3 operator+( const Vector_3 &left, const Vector_3 &right )
{
  return Vector_3( left.vec[0] + right.vec[0], left.vec[1] + right.vec[1], left.vec[2] + right.vec[2] );
}

Vector_3 operator-( const Vector_3 &left, const Vector_3 &right )
{
  return Vector_3( left.vec[0] - right.vec[0], left.vec[1] - right.vec[1], left.vec[2] - right.vec[2] );
}

Vector_3& Vector_3::operator+=( const Vector_3 &source )
{
  vec[0] += source.vec[0];
  vec[1] += source.vec[1];
  vec[2] += source.vec[2];
 
  return *this;
}

Vector_3& Vector_3::operator-=( const Vector_3 &source )
{
  vec[0] -= source.vec[0];
  vec[1] -= source.vec[1];
  vec[2] -= source.vec[2];

  return *this;
}

Vector_3& Vector_3::operator*=( const double &val )
{
  vec[0] *= val;
  vec[1] *= val;
  vec[2] *= val;

  return *this;
}

void Vector_3::print(std::ostream& os, const std::string& delimiter) const
{
  os<<std::setprecision(12);
  os<<vec[0]<<delimiter<<vec[1]<<delimiter<<vec[2]<<std::endl;
}

void Vector_3::gen_rand(const double &left, const double &right)
{
  std::random_device rd;
  std::mt19937_64 gen( rd() );
  std::uniform_real_distribution<double> dis(left, right);
  vec[0] = dis(gen);
  vec[1] = dis(gen);
  vec[2] = dis(gen);
}

double Vector_3::normalize()
{
  const double len = norm2();
  const double inv_len = 1.0 / len;
  vec[0] *= inv_len;
  vec[1] *= inv_len;
  vec[2] *= inv_len;
  
  return len;
}

double Vector_3::dot_product( const Vector_3 &source ) const
{
  return vec[0]*source.vec[0] + vec[1]*source.vec[1] + vec[2]*source.vec[2];
}

double Vec3::dist( const Vector_3 &a, const Vector_3 &b )
{
  const double dist_x = a.x() - b.x();
  const double dist_y = a.y() - b.y();
  const double dist_z = a.z() - b.z();
  return std::sqrt( dist_x*dist_x + dist_y*dist_y + dist_z*dist_z );
}

double Vec3::dot_product( const Vector_3 &a, const Vector_3 &b )
{
  return a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
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

Vector_3 Vec3::cross_product( const Vector_3 &a, const Vector_3 &b )
{
  return Vector_3( a(1) * b(2) - a(2) * b(1), 
      a(2) * b(0) - a(0) * b(2), a(0) * b(1) - a(1) * b(0) );
}

Vector_3 operator*( const double &val, const Vector_3 &source )
{
  return Vector_3( source.x() * val, source.y() * val, source.z() * val );
}

Vector_3 Vec3::normalize( const Vector_3 &val )
{
  const double len = val.norm2();
  return (1.0/len) * val;
}

// EOF
