#include "Tensor4_2D.hpp"

Tensor4_2D::Tensor4_2D()
{
  for(int ii=0; ii<16; ++ii) ten[ii] = 0.0;

  for(int aa=0; aa<2; ++aa)
  {
    for(int bb=0; bb<2; ++bb) ten[8*aa+4*bb+2*aa+bb] = 1.0;
  }
}

Tensor4_2D::Tensor4_2D( const Tensor4_2D &source )
{
  for(int ii=0; ii<16; ++ii) ten[ii] = source(ii);
}

void Tensor4_2D::print() const
{
  std::cout<<"Tensor4_2D: \n";
  for(int kk=0; kk<2; ++kk)
  {
    for(int ll=0; ll<2; ++ll)
    {
      std::cout<<"k = "<<kk<<"\tl = "<<ll<<'\n';
      for(int ii=0; ii<2; ++ii)
      {
        for(int jj=0; jj<2; ++jj)
        {
          std::cout<<"i = "<<ii<<'\t'<<"j = "<<jj<<'\t'
            <<std::setprecision(6)<<ten[8*ii+4*jj+2*kk+ll]<<'\t';
        }
        std::cout<<'\n';
      }
      std::cout<<'\n';
    }
  }
}

void Tensor4_2D::gen_id()
{
  for(int ii=0; ii<16; ++ii) ten[ii] = 0.0;

  for(int aa=0; aa<2; ++aa)
  {
    for(int bb=0; bb<2; ++bb) ten[8*aa+4*bb+2*aa+bb] = 1.0;
  }
}

void Tensor4_2D::gen_symm_id()
{
  for(int ii=0; ii<16; ++ii) ten[ii] = 0.0;
  for(int aa=0; aa<2; ++aa)
  {
    for(int bb=0; bb<2; ++bb)
    {
      ten[ 8 * aa + 4 * bb + 2 * aa + bb ] += 0.5;
      ten[ 8 * aa + 4 * bb + 2 * bb + aa ] += 0.5;
    }
  }
}

void Tensor4_2D::gen_rand(const double &min, const double &max)
{
  std::random_device rd;
  std::mt19937_64 gen( rd() );
  std::uniform_real_distribution<double> dis(min, max);
  for(int ii=0; ii<16; ++ii) ten[ii] = dis(gen);
}

void Tensor4_2D::gen_zero()
{
  for(int ii=0; ii<16; ++ii) ten[ii] = 0.0;
}

void Tensor4_2D::scale( const double &val )
{
  for(int ii=0; ii<16; ++ii) ten[ii] *= val;
}

void Tensor4_2D::PY( const Tensor4_2D &input )
{
  for(int ii=0; ii<16; ++ii) ten[ii] += input(ii);
}

void Tensor4_2D::AXPY( const double &val, const Tensor4_2D &input )
{
  for(int ii=0; ii<16; ++ii) ten[ii] += val * input(ii);
}

void Tensor4_2D::add_OutProduct( const double &val, const Tensor2_2D &mleft,
            const Tensor2_2D &mright )
{
  for(int ii=0; ii<2; ++ii)
  {
    for(int jj=0; jj<2; ++jj)
    {
      for(int kk=0; kk<2; ++kk)
      {
        for(int ll=0; ll<2; ++ll)
          ten[8*ii+4*jj+2*kk+ll] += val * mleft(2*ii+jj) * mright(2*kk+ll);
      }
    }
  }
}

void Tensor4_2D::add_SymmProduct( const double &val, const Tensor2_2D &mleft,
    const Tensor2_2D &mright )
{
  for(int ii=0; ii<2; ++ii)
  {
    for(int jj=0; jj<2; ++jj)
    {
      for(int kk=0; kk<2; ++kk)
      {
        for(int ll=0; ll<2; ++ll)
        {
          ten[8*ii+4*jj+2*kk+ll] += val * 0.5 * ( mleft(2*ii+kk) * mright(2*jj+ll)
              + mleft(2*ii+ll) * mright(2*jj+kk) );
        }
      }
    }
  }
}

void Tensor4_2D::MatMult_1( const Tensor2_2D &source )
{
  double temp[16];
  for(int m=0; m<8; ++m)
  {
    temp[m] = source(0) * ten[m] + source(1) * ten[8+m];
    temp[8+m] = source(2) * ten[m] + source(3) * ten[8+m];
  }
  for(int m=0; m<16; ++m) ten[m] = temp[m];
}

void Tensor4_2D::MatMult_2( const Tensor2_2D &source )
{
  double temp[16];
  // i = 0
  for(int m=0; m<4; ++m)
  {
    temp[m] = source(0) * ten[m] + source(1) * ten[m+4];
    temp[m+4] = source(2) * ten[m] + source(3) * ten[m+4];
  }
  // i = 1
  for(int m=8; m<12; ++m)
  {
    temp[m] = source(0) * ten[m] + source(1) * ten[m+4];
    temp[m+4] = source(2) * ten[m] + source(3) * ten[m+4];
  }
  for(int m=0; m<16; ++m) ten[m] = temp[m];
}

void Tensor4_2D::MatMult_3( const Tensor2_2D &source )
{
  double temp[16];
  for(int m=0; m<=12; m=m+4)
  {
    temp[m] = source(0) * ten[m] + source(1) * ten[m+2];
    temp[m+1] = source(0) * ten[m+1] + source(1) * ten[m+3];
    temp[m+2] = source(2) * ten[m] + source(3) * ten[m+2];
    temp[m+3] = source(2) * ten[m+1] + source(3) * ten[m+3];
  }
  for(int m=0; m<16; ++m) ten[m] = temp[m];
}

void Tensor4_2D::MatMult_4( const Tensor2_2D &source )
{
  double temp[16];
  for(int m=0; m<=14; m=m+2)
  {
    temp[m] = source(0) * ten[m] + source(1) * ten[m+1];
    temp[m+1] = source(2) * ten[m] + source(3) * ten[m+1];
  }
  for(int m=0; m<16; ++m) ten[m] = temp[m];
}

void Tensor4_2D::LeftContraction( const Tensor2_2D &source, Tensor2_2D &out ) const
{
  out(0) = ten[0] * source(0) + ten[4] * source(1) + ten[8] * source(2) + ten[12] * source(3);
  out(1) = ten[1] * source(0) + ten[5] * source(1) + ten[9] * source(2) + ten[13] * source(3);
  out(2) = ten[2] * source(0) + ten[6] * source(1) + ten[10] * source(2) + ten[14] * source(3);
  out(3) = ten[3] * source(0) + ten[7] * source(1) + ten[11] * source(2) + ten[15] * source(3);
}

void Tensor4_2D::RightContraction( const Tensor2_2D &source, Tensor2_2D &out ) const
{
  out(0) = ten[0] * source(0) + ten[1] * source(1) + ten[2] * source(2) + ten[3] * source(3);
  out(1) = ten[4] * source(0) + ten[5] * source(1) + ten[6] * source(2) + ten[7] * source(3);
  out(2) = ten[8] * source(0) + ten[9] * source(1) + ten[10] * source(2) + ten[11] * source(3);
  out(3) = ten[12] * source(0) + ten[13] * source(1) + ten[14] * source(2) + ten[15] * source(3);
}

double Tensor4_2D::LnRContraction( const Tensor2_2D &Left,
            const Tensor2_2D &Right ) const
{
  return Left(0) * (ten[0]*Right(0) + ten[1] * Right(1) + ten[2] * Right(2) + ten[3] * Right(3) ) + Left(1) * (ten[4]*Right(0) + ten[5] * Right(1) + ten[6] * Right(2) + ten[7] * Right(3) ) + Left(2) * (ten[8]*Right(0) + ten[9] * Right(1) + ten[10] * Right(2) + ten[11] * Right(3) ) + Left(3) * (ten[12]*Right(0) + ten[13] * Right(1) + ten[14] * Right(2) + ten[15] * Right(3));
}

double Tensor4_2D::Ten4Contraction( const Tensor4_2D &input ) const
{
  double sum = 0.0;
  for(int ii=0; ii<16; ++ii) sum += input(ii) * ten[ii];

  return sum;
}

// EOF
