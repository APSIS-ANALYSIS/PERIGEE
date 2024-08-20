#include "SymmTensor4_3D.hpp"

SymmTensor4_3D::SymmTensor4_3D()
{
  gen_symm_id();
}

SymmTensor4_3D::SymmTensor4_3D( const std::array<double,21> &source )
{
  for( int ii=0; ii<21; ++ii ) ten[ii] = source[ii];
}

Tensor4_3D SymmTensor4_3D::full() const
{
  Tensor4_3D out;
  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      for(int kk=0; kk<3; ++kk)
      {
        out(ii, jj, kk, 0) = ten[ Voigt_notation(ii, jj, kk, 0)];
        out(ii, jj, kk, 1) = ten[ Voigt_notation(ii, jj, kk, 1)];
        out(ii, jj, kk, 2) = ten[ Voigt_notation(ii, jj, kk, 2)];
      }
    }
  }
  return out;
}

double& SymmTensor4_3D::operator()(const int &ii, const int &jj, const int &kk, const int &ll)
{
  return ten[ Voigt_notation(ii,jj,kk,ll) ];
}

const double& SymmTensor4_3D::operator()(const int &ii, const int &jj, const int &kk, const int &ll) const
{
  return ten[ Voigt_notation(ii,jj,kk,ll) ];
}

bool SymmTensor4_3D::is_identical(const Tensor4_3D &source, const double &tol) const
{
  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      for(int kk=0; kk<3; ++kk)
      {
        if( std::abs( source(ii,jj,kk,0) - ten[Voigt_notation(ii,jj,kk,0)] ) > tol ) return false;
        if( std::abs( source(ii,jj,kk,1) - ten[Voigt_notation(ii,jj,kk,1)] ) > tol ) return false;
        if( std::abs( source(ii,jj,kk,2) - ten[Voigt_notation(ii,jj,kk,2)] ) > tol ) return false;
      }
    }
  }
  return true;
}

bool SymmTensor4_3D::is_identical(const SymmTensor4_3D &source, const double &tol) const
{
  for(int ii=0; ii<21; ++ii)
    if( std::abs(source(ii) - ten[ii]) > tol ) return false;
  return true;
}

void SymmTensor4_3D::gen_rand(const double &left, const double &right)
{
  std::random_device rd;
  std::mt19937_64 gen( rd() );
  std::uniform_real_distribution<double> dis(left, right);
  for(int ii=0; ii<21; ++ii) ten[ii] = dis(gen);
}

void SymmTensor4_3D::gen_zero()
{
  for(int ii=0; ii<21; ++ii) ten[ii] = 0.0;
}

void SymmTensor4_3D::gen_symm_id()
{
  gen_zero();
  ten[0] = 1.0; ten[6] = 1.0; ten[11] = 1.0;
  ten[15] = 0.5; ten[18] = 0.5; ten[20] = 0.5;
}

SymmTensor4_3D& SymmTensor4_3D::operator= (const SymmTensor4_3D &source)
{
  if(this == &source) return *this;

  for(int ii=0; ii<21; ++ii) ten[ii] = source(ii);

  return *this;
}

void SymmTensor4_3D::gen_Ptilde( const SymmTensor2_3D &invC )
{
  gen_zero();
  add_SymmProduct( 1.0, invC, invC );
  add_OutProduct( -1.0/3.0, invC );
}

void SymmTensor4_3D::print() const
{
  std::cout<<"SymmTensor4_3D: \n";
  for(int kk=0; kk<3; ++kk)
  {
    for(int ll=0; ll<3; ++ll)
    {
      std::cout<<"k = "<<kk<<"\tl = "<<ll<<'\n';
      for(int ii=0; ii<3; ++ii)
      {
        for(int jj=0; jj<3; ++jj)
        {
          std::cout<<"i = "<<ii<<'\t'<<"j = "<<jj<<'\t'
            <<std::setprecision(6)<<ten[ Voigt_notation(ii, jj, kk, ll) ]<<'\t';
        }
        std::cout<<'\n';
      }
      std::cout<<'\n';
    }
  }  
}

void SymmTensor4_3D::print_in_mat() const
{
  std::cout<<"SymmTensor4_3D: \n\n";
  for ( int ii=0; ii<3; ++ii )
  {
    for( int jj=0; jj<3; ++jj )
    {
      for( int kk=0; kk<3; ++kk )
      { 
        for ( int ll=0; ll<3; ++ll )
        {
          std::cout << std::setprecision(6) << std::setw(12) << std::left << std::setfill(' ') 
            << ten[ Voigt_notation(ii, jj, kk, ll) ] << " ";
        }
        std::cout<<"\t";
      }  
      std::cout<<'\n';
    }
    std::cout<<"\n";      
  }
}

SymmTensor4_3D operator+( const SymmTensor4_3D &left, const SymmTensor4_3D &right)
{
  SymmTensor4_3D result;
  for(int ii=0; ii<21; ++ii) result.ten[ii] = left.ten[ii] + right.ten[ii];

  return result;
}

SymmTensor4_3D operator-( const SymmTensor4_3D &left, const SymmTensor4_3D &right)
{
  SymmTensor4_3D result;
  for(int ii=0; ii<21; ++ii) result.ten[ii] = left.ten[ii] - right.ten[ii];

  return result;
}

SymmTensor4_3D& SymmTensor4_3D::operator+=( const SymmTensor4_3D &source )
{
  for(int ii=0; ii<21; ++ii) ten[ii] += source(ii);
  return *this;
}

SymmTensor4_3D& SymmTensor4_3D::operator-=( const SymmTensor4_3D &source )
{
  for(int ii=0; ii<21; ++ii) ten[ii] -= source(ii);
  return *this;
}

SymmTensor4_3D& SymmTensor4_3D::operator*=( const double &val )
{
  for(int ii=0; ii<21; ++ii) ten[ii] *= val;
  return *this;
}

void SymmTensor4_3D::add_OutProduct( const double &val, const SymmTensor2_3D &mmat )
{
  ten[0] += val * mmat(0) * mmat(0);
  ten[1] += val * mmat(0) * mmat(1);
  ten[2] += val * mmat(0) * mmat(2);
  ten[3] += val * mmat(0) * mmat(3);
  ten[4] += val * mmat(0) * mmat(4);
  ten[5] += val * mmat(0) * mmat(5);

  ten[6]  += val * mmat(1) * mmat(1);
  ten[7]  += val * mmat(1) * mmat(2);
  ten[8]  += val * mmat(1) * mmat(3);
  ten[9]  += val * mmat(1) * mmat(4);
  ten[10] += val * mmat(1) * mmat(5);

  ten[11] += val * mmat(2) * mmat(2);
  ten[12] += val * mmat(2) * mmat(3);
  ten[13] += val * mmat(2) * mmat(4);
  ten[14] += val * mmat(2) * mmat(5);

  ten[15] += val * mmat(3) * mmat(3);
  ten[16] += val * mmat(3) * mmat(4);
  ten[17] += val * mmat(3) * mmat(5);

  ten[18] += val * mmat(4) * mmat(4);
  ten[19] += val * mmat(4) * mmat(5);

  ten[20] += val * mmat(5) * mmat(5);
}

void SymmTensor4_3D::add_OutProduct( const double &val, const SymmTensor2_3D &mleft,
    const SymmTensor2_3D &mright )
{
  ten[0] += val * mleft(0) * mright(0);
  ten[1] += val * mleft(0) * mright(1);
  ten[2] += val * mleft(0) * mright(2);
  ten[3] += val * mleft(0) * mright(3);
  ten[4] += val * mleft(0) * mright(4);
  ten[5] += val * mleft(0) * mright(5);

  ten[6]  += val * mleft(1) * mright(1);
  ten[7]  += val * mleft(1) * mright(2);
  ten[8]  += val * mleft(1) * mright(3);
  ten[9]  += val * mleft(1) * mright(4);
  ten[10] += val * mleft(1) * mright(5);

  ten[11] += val * mleft(2) * mright(2);
  ten[12] += val * mleft(2) * mright(3);
  ten[13] += val * mleft(2) * mright(4);
  ten[14] += val * mleft(2) * mright(5);

  ten[15] += val * mleft(3) * mright(3);
  ten[16] += val * mleft(3) * mright(4);
  ten[17] += val * mleft(3) * mright(5);

  ten[18] += val * mleft(4) * mright(4);
  ten[19] += val * mleft(4) * mright(5);

  ten[20] += val * mleft(5) * mright(5);
}

void SymmTensor4_3D::add_SymmOutProduct( const double &val, const Vector_3 &vec1, 
    const Vector_3 &vec2 )
{
  ten[0]  += val * 4.0 * vec1(0) * vec2(0) * vec1(0) * vec2(0);

  ten[1]  += val * 4.0 * vec1(0) * vec2(0) * vec1(1) * vec2(1);

  ten[2]  += val * 4.0 * vec1(0) * vec2(0) * vec1(2) * vec2(2);

  ten[3]  += val * 2.0 * ( vec1(0) * vec2(0) * vec1(1) * vec2(2)
      + vec1(0) * vec2(0) * vec1(2) * vec2(1) );

  ten[4]  += val * 2.0 * ( vec1(0) * vec2(0) * vec1(0) * vec2(2)
      + vec1(0) * vec2(0) * vec1(2) * vec2(0) );

  ten[5]  += val * 2.0 * ( vec1(0) * vec2(0) * vec1(0) * vec2(1)
      + vec1(0) * vec2(0) * vec1(1) * vec2(0) );

  ten[6]  += val * 4.0 * vec1(1) * vec2(1) * vec1(1) * vec2(1);

  ten[7]  += val * 4.0 * vec1(1) * vec2(1) * vec1(2) * vec2(2);

  ten[8]  += val * 2.0 * ( vec1(1) * vec2(1) * vec1(1) * vec2(2)
      + vec1(1) * vec2(1) * vec1(2) * vec2(1) );

  ten[9]  += val * 2.0 * ( vec1(0) * vec2(2) * vec1(1) * vec2(1)
      + vec1(2) * vec2(0) * vec1(1) * vec2(1) );

  ten[10] += val * 2.0 * ( vec1(1) * vec2(0) * vec1(1) * vec2(1)
      + vec1(0) * vec2(1) * vec1(1) * vec2(1) );

  ten[11] += val * 4.0 * vec1(2) * vec2(2) * vec1(2) * vec2(2);

  ten[12] += val * 2.0 * ( vec1(1) * vec2(2) * vec1(2) * vec2(2)
      + vec1(2) * vec2(1) * vec1(2) * vec2(2) );

  ten[13] += val * 2.0 * ( vec1(2) * vec2(0) * vec1(2) * vec2(2)
      + vec1(0) * vec2(2) * vec1(2) * vec2(2) );

  ten[14] += val * 2.0 * ( vec1(0) * vec2(1) * vec1(2) * vec2(2)
      + vec1(1) * vec2(0) * vec1(2) * vec2(2) );

  ten[15] += val * ( vec1(1) * vec2(2) * vec1(1) * vec2(2)
      + vec1(1) * vec2(2) * vec1(2) * vec2(1)
      + vec1(2) * vec2(1) * vec1(1) * vec2(2)
      + vec1(2) * vec2(1) * vec1(2) * vec2(1) );

  ten[16] += val * ( vec1(0) * vec2(2) * vec1(1) * vec2(2)
      + vec1(0) * vec2(2) * vec1(2) * vec2(1)
      + vec1(2) * vec2(0) * vec1(1) * vec2(2)
      + vec1(2) * vec2(0) * vec1(2) * vec2(1) );

  ten[17] += val * ( vec1(0) * vec2(1) * vec1(1) * vec2(2)
      + vec1(0) * vec2(1) * vec1(2) * vec2(1)
      + vec1(1) * vec2(0) * vec1(1) * vec2(2)
      + vec1(1) * vec2(0) * vec1(2) * vec2(1) );

  ten[18] += val * ( vec1(0) * vec2(2) * vec1(0) * vec2(2)
      + vec1(0) * vec2(2) * vec1(2) * vec2(0)
      + vec1(2) * vec2(0) * vec1(0) * vec2(2)
      + vec1(2) * vec2(0) * vec1(2) * vec2(0) );

  ten[19] += val * ( vec1(0) * vec2(1) * vec1(0) * vec2(2)
      + vec1(0) * vec2(1) * vec1(2) * vec2(0)
      + vec1(1) * vec2(0) * vec1(0) * vec2(2)
      + vec1(1) * vec2(0) * vec1(2) * vec2(0) );

  ten[20] += val * ( vec1(0) * vec2(1) * vec1(0) * vec2(1)
      + vec1(0) * vec2(1) * vec1(1) * vec2(0)
      + vec1(1) * vec2(0) * vec1(0) * vec2(1)
      + vec1(1) * vec2(0) * vec1(1) * vec2(0) );
}

void SymmTensor4_3D::add_SymmProduct( const double &val, const SymmTensor2_3D &mleft,
    const SymmTensor2_3D &mright )
{
  ten[0]  += val * ( mleft(0) * mright(0) );
  ten[1]  += val * ( mleft(5) * mright(5) );
  ten[2]  += val * ( mleft(4) * mright(4) );

  ten[3]  += val * 0.5 * ( mleft(5) * mright(4) + mleft(4) * mright(5) );
  ten[4]  += val * 0.5 * ( mleft(0) * mright(4) + mleft(4) * mright(0) );
  ten[5]  += val * 0.5 * ( mleft(0) * mright(5) + mleft(5) * mright(0) );

  ten[6]  += val * ( mleft(1) * mright(1) );
  ten[7]  += val * ( mleft(3) * mright(3) );

  ten[8]  += val * 0.5 * ( mleft(1) * mright(3) + mleft(3) * mright(1) );
  ten[9]  += val * 0.5 * ( mleft(5) * mright(3) + mleft(5) * mright(3) );
  ten[10] += val * 0.5 * ( mleft(5) * mright(1) + mleft(5) * mright(1) );

  ten[11] += val * ( mleft(2) * mright(2) );

  ten[12] += val * 0.5 * ( mleft(3) * mright(2) + mleft(3) * mright(2) );
  ten[13] += val * 0.5 * ( mleft(4) * mright(2) + mleft(4) * mright(2) );
  ten[14] += val * 0.5 * ( mleft(4) * mright(3) + mleft(4) * mright(3) );
  ten[15] += val * 0.5 * ( mleft(1) * mright(2) + mleft(3) * mright(3) );
  ten[16] += val * 0.5 * ( mleft(5) * mright(2) + mleft(4) * mright(3) );
  ten[17] += val * 0.5 * ( mleft(5) * mright(3) + mleft(4) * mright(1) );
  ten[18] += val * 0.5 * ( mleft(0) * mright(2) + mleft(4) * mright(4) );
  ten[19] += val * 0.5 * ( mleft(0) * mright(3) + mleft(4) * mright(5) );
  ten[20] += val * 0.5 * ( mleft(0) * mright(1) + mleft(5) * mright(5) );
}

void SymmTensor4_3D::add_SymmOutProduct( const double &val, const SymmTensor2_3D &mleft,
    const SymmTensor2_3D &mright )
{
  ten[0]  += val * ( mleft(0) * mright(0) + mleft(0) * mright(0) );
  ten[1]  += val * ( mleft(0) * mright(1) + mleft(1) * mright(0) );
  ten[2]  += val * ( mleft(0) * mright(2) + mleft(2) * mright(0) );
  ten[3]  += val * ( mleft(0) * mright(3) + mleft(3) * mright(0) );
  ten[4]  += val * ( mleft(0) * mright(4) + mleft(4) * mright(0) );
  ten[5]  += val * ( mleft(0) * mright(5) + mleft(5) * mright(0) );
  ten[6]  += val * ( mleft(1) * mright(1) + mleft(1) * mright(1) );
  ten[7]  += val * ( mleft(1) * mright(2) + mleft(2) * mright(1) );
  ten[8]  += val * ( mleft(1) * mright(3) + mleft(3) * mright(1) );
  ten[9]  += val * ( mleft(4) * mright(1) + mleft(1) * mright(4) );
  ten[10] += val * ( mleft(5) * mright(1) + mleft(1) * mright(5) );
  ten[11] += val * ( mleft(2) * mright(2) + mleft(2) * mright(2) );
  ten[12] += val * ( mleft(3) * mright(2) + mleft(2) * mright(3) );
  ten[13] += val * ( mleft(4) * mright(2) + mleft(2) * mright(4) );
  ten[14] += val * ( mleft(5) * mright(2) + mleft(2) * mright(5) );
  ten[15] += val * ( mleft(3) * mright(3) + mleft(3) * mright(3) );
  ten[16] += val * ( mleft(4) * mright(3) + mleft(3) * mright(4) );
  ten[17] += val * ( mleft(5) * mright(3) + mleft(3) * mright(5) );
  ten[18] += val * ( mleft(4) * mright(4) + mleft(4) * mright(4) );
  ten[19] += val * ( mleft(5) * mright(4) + mleft(4) * mright(5) );
  ten[20] += val * ( mleft(5) * mright(5) + mleft(5) * mright(5) );
}

void SymmTensor4_3D::TenPMult( const Tensor4_3D &P )
{
  double temp[21] = {0.0};

  constexpr int value_left[21]  = {0,  0,  0,  0,  0, 0, 36, 36, 36, 36, 36, 72, 72, 72, 72, 45, 45, 45, 18, 9, 9};
  constexpr int value_right[21] = {0, 36, 72, 45, 18, 9, 36, 72, 45, 18,  9, 72, 45, 18,  9, 45, 18, 9, 18, 18, 9};

  for(int kk=0; kk<21; ++kk)
  {
    temp[kk] += P( value_left[kk] ) * ten[0] * P( value_right[kk]   );
    temp[kk] += P( value_left[kk] ) * ten[1] * P( value_right[kk]+4 );
    temp[kk] += P( value_left[kk] ) * ten[2] * P( value_right[kk]+8 );

    temp[kk] += 2.0 * P( value_left[kk] ) * ten[3] * P( value_right[kk]+5 );
    temp[kk] += 2.0 * P( value_left[kk] ) * ten[4] * P( value_right[kk]+2 );
    temp[kk] += 2.0 * P( value_left[kk] ) * ten[5] * P( value_right[kk]+1 );

    temp[kk] += P( value_left[kk]+4 ) * ten[1] * P( value_right[kk]   );
    temp[kk] += P( value_left[kk]+4 ) * ten[6] * P( value_right[kk]+4 );
    temp[kk] += P( value_left[kk]+4 ) * ten[7] * P( value_right[kk]+8 );

    temp[kk] += 2.0 * P( value_left[kk]+4 ) * ten[8]  * P( value_right[kk]+5 );
    temp[kk] += 2.0 * P( value_left[kk]+4 ) * ten[9]  * P( value_right[kk]+2 );
    temp[kk] += 2.0 * P( value_left[kk]+4 ) * ten[10] * P( value_right[kk]+1 );

    temp[kk] += P( value_left[kk]+8 ) * ten[2]  * P( value_right[kk]   );
    temp[kk] += P( value_left[kk]+8 ) * ten[7]  * P( value_right[kk]+4 );
    temp[kk] += P( value_left[kk]+8 ) * ten[11] * P( value_right[kk]+8 );

    temp[kk] += 2.0 * P( value_left[kk]+8 ) * ten[12] * P( value_right[kk]+5 );
    temp[kk] += 2.0 * P( value_left[kk]+8 ) * ten[13] * P( value_right[kk]+2 );
    temp[kk] += 2.0 * P( value_left[kk]+8 ) * ten[14] * P( value_right[kk]+1 );

    temp[kk] += 2.0 * P( value_left[kk]+5 ) * ten[3]  * P( value_right[kk]   );
    temp[kk] += 2.0 * P( value_left[kk]+5 ) * ten[8]  * P( value_right[kk]+4 );
    temp[kk] += 2.0 * P( value_left[kk]+5 ) * ten[12] * P( value_right[kk]+8 );

    temp[kk] += 4.0 * P( value_left[kk]+5 ) * ten[15] * P( value_right[kk]+5 );
    temp[kk] += 4.0 * P( value_left[kk]+5 ) * ten[16] * P( value_right[kk]+2 );
    temp[kk] += 4.0 * P( value_left[kk]+5 ) * ten[17] * P( value_right[kk]+1 );

    temp[kk] += 2.0 * P( value_left[kk]+2 ) * ten[4]  * P( value_right[kk]   );
    temp[kk] += 2.0 * P( value_left[kk]+2 ) * ten[9]  * P( value_right[kk]+4 );
    temp[kk] += 2.0 * P( value_left[kk]+2 ) * ten[13] * P( value_right[kk]+8 );

    temp[kk] += 4.0 * P( value_left[kk]+2 ) * ten[16] * P( value_right[kk]+5 );
    temp[kk] += 4.0 * P( value_left[kk]+2 ) * ten[18] * P( value_right[kk]+2 );
    temp[kk] += 4.0 * P( value_left[kk]+2 ) * ten[19] * P( value_right[kk]+1 );

    temp[kk] += 2.0 * P( value_left[kk]+1 ) * ten[5]  * P( value_right[kk]   );
    temp[kk] += 2.0 * P( value_left[kk]+1 ) * ten[10] * P( value_right[kk]+4 );
    temp[kk] += 2.0 * P( value_left[kk]+1 ) * ten[14] * P( value_right[kk]+8 );

    temp[kk] += 4.0 *  P( value_left[kk]+1 ) * ten[17] * P( value_right[kk]+5 );
    temp[kk] += 4.0 *  P( value_left[kk]+1 ) * ten[19] * P( value_right[kk]+2 );
    temp[kk] += 4.0 *  P( value_left[kk]+1 ) * ten[20] * P( value_right[kk]+1 );
  }

  for(int ii=0; ii<21; ++ii) ten[ii] = temp[ii];
}

SymmTensor4_3D operator*( const double &val, const SymmTensor4_3D &input )
{
  return SymmTensor4_3D( std::array<double,21>{{ 
      val * input(0),  val * input(1),  val * input(2),  
      val * input(3),  val * input(4),  val * input(5),
      val * input(6),  val * input(7),  val * input(8),
      val * input(9),  val * input(10), val * input(11), 
      val * input(12), val * input(13), val * input(14),
      val * input(15), val * input(16), val * input(17), 
      val * input(18), val * input(19), val * input(20) }});
}

SymmTensor4_3D STen4::gen_zero()
{
  constexpr std::array<double,21> temp {{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};
  return SymmTensor4_3D(temp);
}

SymmTensor4_3D STen4::gen_symm_id()
{
  constexpr std::array<double,21> temp {{ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.5 }};
  return SymmTensor4_3D(temp);
}

SymmTensor4_3D STen4::gen_Ptilde( const SymmTensor2_3D &invC )
{     
  SymmTensor4_3D out = STen4::gen_zero();
  out.add_SymmProduct( 1.0, invC, invC );
  out.add_OutProduct( -1.0/3.0, invC );
  return out;
}

// EOF
