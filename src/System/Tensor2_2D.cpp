#include "Tensor2_2D.hpp"

Tensor2_2D::Tensor2_2D()
{
  mat[0] = 1.0; mat[1] = 0.0; mat[2] = 0.0; mat[3] = 1.0;
}

Tensor2_2D::Tensor2_2D( const Tensor2_2D &source )
{
  mat[0] = source(0); mat[1] = source(1);
  mat[2] = source(2); mat[3] = source(3);
}

Tensor2_2D::Tensor2_2D( const double &a11, const double &a12,
    const double &a21, const double &a22 )
{
  mat[0] = a11; mat[1] = a12;
  mat[2] = a21; mat[3] = a22;
}

void Tensor2_2D::copy( const Tensor2_2D &source )
{
  mat[0] = source(0);
  mat[1] = source(1);
  mat[2] = source(2);
  mat[3] = source(3);
}

void Tensor2_2D::copy( double source[4] )
{
  mat[0] = source[0];
  mat[1] = source[1];
  mat[2] = source[2];
  mat[3] = source[3];
}

void Tensor2_2D::gen_zero()
{
  mat[0] = 0.0; mat[1] = 0.0; mat[2] = 0.0; mat[3] = 0.0;
}

void Tensor2_2D::gen_id()
{
  mat[0] = 1.0; mat[1] = 0.0; mat[2] = 0.0; mat[3] = 1.0;
}

void Tensor2_2D::gen_rand(const double &min, const double &max)
{
  std::random_device rd;
  std::mt19937_64 gen( rd() );
  std::uniform_real_distribution<double> dis(min, max);
  for(int ii=0; ii<4; ++ii) mat[ii] = dis(gen);
}

void Tensor2_2D::gen_hilb()
{
  for(int ii=0; ii<2; ++ii)
    for(int jj=0; jj<2; ++jj)
      mat[ii*2+jj] = 1.0 / (ii + jj + 1.0);
}

void Tensor2_2D::transpose()
{
  double temp = mat[1]; mat[1] = mat[2]; mat[2] = temp;
}

void Tensor2_2D::inverse()
{
  const double invdetA = 1.0 / det();
  double temp[4];
  temp[0] = invdetA * mat[3];
  temp[1] = -1.0 * invdetA * mat[1];
  temp[2] = -1.0 * invdetA * mat[2];
  temp[3] = invdetA * mat[0];

  mat[0] = temp[0]; mat[1] = temp[1]; mat[2] = temp[2]; mat[3] = temp[3];
}

void Tensor2_2D::scale( const double &val )
{
  for(int ii=0; ii<4; ++ii) mat[ii] *= val;
}

void Tensor2_2D::AXPY( const double &val, const Tensor2_2D &source )
{
  for(int ii=0; ii<4; ++ii) mat[ii] += val * source(ii);
}

void Tensor2_2D::PY( const Tensor2_2D &source )
{
  for(int ii=0; ii<4; ++ii) mat[ii] += source(ii);
}

double Tensor2_2D::det() const
{
  return mat[0] * mat[3] - mat[1] * mat[2];
}

double Tensor2_2D:: VecMatVec( const double * const &x,
    const double * const &y ) const
{
  return x[0] * (mat[0]*y[0] + mat[1]*y[1]) + x[1] * (mat[2]*y[0] + mat[3]*y[1]);
}

void Tensor2_2D::VecMult( const double * const &x, double * const &y ) const
{
  y[0] = mat[0] * x[0] + mat[1] * x[1];
  y[1] = mat[2] * x[0] + mat[3] * x[1];
}

void Tensor2_2D::VecMult( double * const &x ) const
{
  double a, b;
  a = mat[0] * x[0] + mat[1] * x[1];
  b = mat[2] * x[0] + mat[3] * x[1];
  x[0] = a;
  x[1] = b;
}

void Tensor2_2D::MatMult( const Tensor2_2D &mleft, const Tensor2_2D &mright )
{
  double temp[4];
  temp[0] = mleft(0) * mright(0) + mleft(1) * mright(2);
  temp[1] = mleft(0) * mright(1) + mleft(1) * mright(3);
  temp[2] = mleft(2) * mright(0) + mleft(3) * mright(2);
  temp[3] = mleft(2) * mright(1) + mleft(3) * mright(3);

  mat[0] = temp[0]; mat[1] = temp[1];
  mat[2] = temp[2]; mat[3] = temp[3];
}

void Tensor2_2D::MatMultTransposeLeft( const Tensor2_2D &source )
{
  double temp[4];
  temp[0] = source(0) * source(0) + source(2) * source(2);
  temp[1] = source(0) * source(1) + source(2) * source(3);
  temp[2] = temp[1];
  temp[3] = source(1) * source(1) + source(3) * source(3);

  mat[0] = temp[0]; mat[1] = temp[1];
  mat[2] = temp[2]; mat[3] = temp[3];
}

void Tensor2_2D::MatMultTransposeRight( const Tensor2_2D &source )
{
  double temp[4];
  temp[0] = source(0) * source(0) + source(1) * source(1);
  temp[1] = source(0) * source(2) + source(1) * source(3);
  temp[2] = temp[1];
  temp[3] = source(2) * source(2) + source(3) * source(3); 

  mat[0] = temp[0]; mat[1] = temp[1];
  mat[2] = temp[2]; mat[3] = temp[3];
}

double Tensor2_2D::MatContraction( const Tensor2_2D &source ) const
{
  return mat[0] * source(0) + mat[1] * source(1) + mat[2] * source(2) + mat[3] * source(3);
}

void Tensor2_2D::print() const
{
  std::cout<<std::setprecision(9)<<mat[0]<<'\t'<<mat[1]<<'\n';
  std::cout<<std::setprecision(9)<<mat[2]<<'\t'<<mat[3]<<'\n';
}

// EOF
