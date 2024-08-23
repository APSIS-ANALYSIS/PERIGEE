#include "Tensor4_3D.hpp"

Tensor4_3D::Tensor4_3D()
{
  for(int ii=0; ii<81; ++ii) ten[ii] = 0.0;

  for(int aa=0; aa<3; ++aa)
  {
    for(int bb=0; bb<3; ++bb) ten[ 27 * aa + 9 * bb + 3 * aa + bb ] = 1.0;
  }
}

Tensor4_3D::Tensor4_3D( const Tensor4_3D &source )
{
  for(int ii=0; ii<81; ++ii) ten[ii] = source(ii);
}

Tensor4_3D::Tensor4_3D( const std::array<double,81> &source )
{
  for(int ii=0; ii<81; ++ii) ten[ii] = source[ii];
}

bool Tensor4_3D::is_identical(const Tensor4_3D &source, const double &tol) const
{
  for(int ii=0; ii<81; ++ii)
    if( std::abs(source(ii) - ten[ii]) > tol ) return false;
  return true;
}

void Tensor4_3D::copy( const Tensor4_3D &source )
{
  for(int ii=0; ii<81; ++ii) ten[ii] = source(ii);
}

Tensor4_3D& Tensor4_3D::operator= (const Tensor4_3D &source)
{
  if(this == &source) return *this;

  for(int ii=0; ii<81; ++ii) ten[ii] = source(ii);
  
  return *this;
}

void Tensor4_3D::print() const
{
  std::cout<<"Tensor4_3D: \n";
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
            <<std::setprecision(6)<<ten[27*ii+9*jj+3*kk+ll]<<'\t';
        }
        std::cout<<'\n';
      }
      std::cout<<'\n';
    }
  }
}

void Tensor4_3D::print_in_mat() const
{
  std::cout<<"Tensor4_3D:\n";
  for ( int ii=0; ii<3; ii++ )
  {
    for( int jj=0; jj<3; jj++ )
    {
      for( int kk=0; kk<3; kk++ )
      { 
        for ( int ll=0; ll<3; ll++ )
        {
          std::cout << std::setprecision(6) << std::setw(12) << std::left << std::setfill(' ') 
          << ten[27*ii+9*jj+3*kk+ll] << " ";
        }
        std::cout<<"\t";
      }  
      std::cout<<'\n';
    }
    std::cout<<"\n";      
  }
}

Tensor4_3D operator+( const Tensor4_3D &left, const Tensor4_3D &right)
{
  Tensor4_3D result;
  for(int ii=0; ii<81; ++ii) result.ten[ii] = left.ten[ii] + right.ten[ii];

  return result;
}

Tensor4_3D operator-( const Tensor4_3D &left, const Tensor4_3D &right)
{
  Tensor4_3D result;
  for(int ii=0; ii<81; ++ii) result.ten[ii] = left.ten[ii] - right.ten[ii];

  return result;
}

Tensor4_3D& Tensor4_3D::operator+=( const Tensor4_3D &source )
{
  for(int ii=0; ii<81; ++ii) ten[ii] += source(ii);
  return *this;
}

Tensor4_3D& Tensor4_3D::operator-=( const Tensor4_3D &source )
{
  for(int ii=0; ii<81; ++ii) ten[ii] -= source(ii);
  return *this;
}

Tensor4_3D& Tensor4_3D::operator*=( const double &val )
{
  for(int ii=0; ii<81; ++ii) ten[ii] *= val;
  return *this;
}

void Tensor4_3D::gen_id()
{
  for(int ii=0; ii<81; ++ii) ten[ii] = 0.0;

  for(int aa=0; aa<3; ++aa)
  {
    for(int bb=0; bb<3; ++bb) ten[ 27 * aa + 9 * bb + 3 * aa + bb ] = 1.0;
  }
}

void Tensor4_3D::gen_symm_id()
{
  for(int ii=0; ii<81; ++ii) ten[ii] = 0.0;
  for(int aa=0; aa<3; ++aa)
  {
    for(int bb=0; bb<3; ++bb) 
    {
      ten[ 27 * aa + 9 * bb + 3 * aa + bb ] += 0.5;
      ten[ 27 * aa + 9 * bb + 3 * bb + aa ] += 0.5; 
    }
  }
}

void Tensor4_3D::gen_proj_dev()
{
  const double pt33 = 1.0 / 3.0;
  for(int ii=0; ii<81; ++ii) ten[ii] = 0.0;
  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      ten[ 27*ii + 9*jj + 3*ii + jj] += 1.0;
      ten[ 27*ii + 9*ii + 3*jj + jj] -= pt33;
    }
  }
}

void Tensor4_3D::gen_P( const Tensor2_3D &C, const Tensor2_3D &invC )
{
  gen_symm_id();
  add_OutProduct( -1.0 / 3.0, invC, C );
}

void Tensor4_3D::gen_Ptilde( const Tensor2_3D &invC )
{
  gen_zero();
  add_SymmProduct(1.0, invC, invC);
  add_OutProduct( -1.0 / 3.0, invC, invC );
}

void Tensor4_3D::gen_rand(const double &left, const double &right)
{
  std::random_device rd;
  std::mt19937_64 gen( rd() );
  std::uniform_real_distribution<double> dis(left, right);
  for(int ii=0; ii<81; ++ii) ten[ii] = dis(gen);
}

void Tensor4_3D::gen_zero()
{
  for(int ii=0; ii<81; ++ii) ten[ii] = 0.0;
}

void Tensor4_3D::scale( const double &val )
{
  for(int ii=0; ii<81; ++ii) ten[ii] *= val;
}

void Tensor4_3D::AXPY( const double &val, const Tensor4_3D &input )
{
  for(int ii=0; ii<81; ++ii) ten[ii] += val * input(ii);
}

void Tensor4_3D::transpose()
{
  std::array<double,81> temp;

  for(int mm=0; mm<9; ++mm)
  {
    temp[   mm] = ten[9*mm];
    temp[ 9+mm] = ten[9*mm+1];
    temp[18+mm] = ten[9*mm+2];
    temp[27+mm] = ten[9*mm+3];
    temp[36+mm] = ten[9*mm+4];
    temp[45+mm] = ten[9*mm+5];
    temp[54+mm] = ten[9*mm+6];
    temp[63+mm] = ten[9*mm+7];
    temp[72+mm] = ten[9*mm+8];
  }

  for(int ii=0; ii<81; ++ii) ten[ii] = temp[ii];
}

void Tensor4_3D::add_OutProduct( const double &val, const Tensor2_3D &mleft,
    const Tensor2_3D &mright )
{
  for(int mm=0; mm<9; ++mm)
  {
    for(int nn=0; nn<9; ++nn)
      ten[9*mm+nn] += val * mleft(mm) * mright(nn);
  }
}

void Tensor4_3D::add_OutProduct( const double &val, const SymmTensor2_3D &mleft,
    const SymmTensor2_3D &mright )
{
  for(int mm=0; mm<9; ++mm)
  {
    const int id = mleft.Voigt_notation(mm);
    ten[9*mm+0] += val * mleft( id ) * mright(0);
    ten[9*mm+1] += val * mleft( id ) * mright(5);
    ten[9*mm+2] += val * mleft( id ) * mright(4);
    ten[9*mm+3] += val * mleft( id ) * mright(5);
    ten[9*mm+4] += val * mleft( id ) * mright(1);
    ten[9*mm+5] += val * mleft( id ) * mright(3);
    ten[9*mm+6] += val * mleft( id ) * mright(4);
    ten[9*mm+7] += val * mleft( id ) * mright(3);
    ten[9*mm+8] += val * mleft( id ) * mright(2);
  }
}

void Tensor4_3D::add_OutProduct( const double &val, const Vector_3 &vec1, 
    const Vector_3 &vec2, const Vector_3 &vec3, const Vector_3 &vec4 )
{
  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      for(int kk=0; kk<3; ++kk)
      {
        for(int ll=0; ll<3; ++ll)
          ten[27*ii+9*jj+3*kk+ll] += val * vec1( ii ) * vec2( jj ) * vec3( kk ) * vec4( ll );
      }
    }
  }
}

void Tensor4_3D::add_SymmOutProduct( const double &val, const Vector_3 &vec1, 
    const Vector_3 &vec2, const Vector_3 &vec3, const Vector_3 &vec4 )
{
  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      for(int kk=0; kk<3; ++kk)
      {
        for(int ll=0; ll<3; ++ll)
        {
          ten[27*ii+9*jj+3*kk+ll] += val * ( vec1( ii ) * vec2( jj ) * vec3( kk ) * vec4( ll )
            + vec1( ii ) * vec2( jj ) * vec3( ll ) * vec4( kk )
            + vec1( jj ) * vec2( ii ) * vec3( kk ) * vec4( ll )
            + vec1( jj ) * vec2( ii ) * vec3( ll ) * vec4( kk ) );
        }
      }
    }
  }
}

void Tensor4_3D::add_SymmProduct( const double &val, const Tensor2_3D &mleft,
    const Tensor2_3D &mright )
{
  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      for(int kk=0; kk<3; ++kk)
      {
        for(int ll=0; ll<3; ++ll)
        {
          ten[27*ii+9*jj+3*kk+ll] += val * 0.5 * ( mleft(3*ii+kk) * mright(3*jj+ll)
              + mleft(3*ii+ll) * mright(3*jj+kk) );
        }
      }
    }
  }
}

void Tensor4_3D::add_SymmOutProduct( const double &val, const Tensor2_3D &mleft,
    const Tensor2_3D &mright )
{
  for(int ii=0; ii<3 ; ++ii)
  {
    for(int jj=0; jj<3 ; ++jj)
    {
      for(int kk=0; kk<3; ++kk)
      {
        for(int ll=0; ll<3; ++ll)
        {
          ten[27*ii + 9*jj + 3*kk + ll] += val *  (mleft(3*ii + jj) * mright(3*kk + ll)
              + mright(3*ii + jj) * mleft(3*kk + ll));
        }
      }
    }
  }
}

void Tensor4_3D::MatMult_1( const Tensor2_3D &source )
{
  double temp[81] {0.0};
  
  for( int m=0; m<27; ++m )
  {
    temp[m]    = source(0) * ten[m] + source(1) * ten[27+m] + source(2) * ten[54+m];
    temp[27+m] = source(3) * ten[m] + source(4) * ten[27+m] + source(5) * ten[54+m];
    temp[54+m] = source(6) * ten[m] + source(7) * ten[27+m] + source(8) * ten[54+m];
  }

  for(int ii=0; ii<81; ++ii) ten[ii] = temp[ii];
}

void Tensor4_3D::MatMult_2( const Tensor2_3D &source )
{
  // 27 I + 3 K + L ranges from 0 - 8, 27 - 35, 54 - 62
  const int index[27] { 0, 1, 2, 3, 4, 5, 6, 7, 8,
    27, 28, 29, 30, 31, 32, 33, 34, 35,
    54, 55, 56, 57, 58, 59, 60, 61, 62 };

  double temp[81] {0.0};
  
  for(int m=0; m<27; ++m)
  {
    const int loc = index[m];
    temp[loc]    = source(0) * ten[loc] + source(1) * ten[loc+9] + source(2) * ten[loc+18];
    temp[loc+9]  = source(3) * ten[loc] + source(4) * ten[loc+9] + source(5) * ten[loc+18]; 
    temp[loc+18] = source(6) * ten[loc] + source(7) * ten[loc+9] + source(8) * ten[loc+18];
  }

  for(int ii=0; ii<81; ++ii) ten[ii] = temp[ii];
}

void Tensor4_3D::MatMult_3( const Tensor2_3D &source )
{
  // 27 I + 9 J + L: ranges from [0-2], [0-2] + 9, [0-2] + 18
  // [0-2] + 27, [0-2] + 36, [0-2] + 45
  // [0-2] + 54, [0-2] + 63, [0-2] + 72
  const int index[27] { 0, 1, 2, 9, 10, 11, 18, 19, 20,
    27, 28, 29, 36, 37, 38, 45, 46, 47,
    54, 55, 56, 63, 64, 65, 72, 73, 74 };

  double temp[81] {0.0};

  for(int m=0; m<27; ++m)
  {
    const int loc = index[m];
    temp[loc]   = source(0) * ten[loc] + source(1) * ten[loc+3] + source(2) * ten[loc+6];
    temp[loc+3] = source(3) * ten[loc] + source(4) * ten[loc+3] + source(5) * ten[loc+6];
    temp[loc+6] = source(6) * ten[loc] + source(7) * ten[loc+3] + source(8) * ten[loc+6];
  } 

  for(int ii=0; ii<81; ++ii) ten[ii] = temp[ii];
}

void Tensor4_3D::MatMult_4( const Tensor2_3D &source )
{
  double temp[81] {0.0};
  
  for(int m=0; m<27; ++m)
  {
    const int loc = 3*m;
    temp[loc]   = source(0) * ten[loc] + source(1) * ten[loc+1] + source(2) * ten[loc+2];
    temp[loc+1] = source(3) * ten[loc] + source(4) * ten[loc+1] + source(5) * ten[loc+2];
    temp[loc+2] = source(6) * ten[loc] + source(7) * ten[loc+1] + source(8) * ten[loc+2];
  }

  for(int ii=0; ii<81; ++ii) ten[ii] = temp[ii];
}

Tensor2_3D Tensor4_3D::LeftContraction( const Tensor2_3D &a ) const
{
  Tensor2_3D out;
  for(int m=0; m<9; ++m)
  {
    out(m) = a(0) * ten[m] + a(3) * ten[m+27] + a(6) * ten[m+54] 
      + a(1) * ten[m+9] + a(4) * ten[m+36] + a(7) * ten[m+63]
      + a(2) * ten[m+18] + a(5) * ten[m+45] + a(8) * ten[m+72];
  }
  return out;
}

Tensor2_3D Tensor4_3D::RightContraction( const Tensor2_3D &a ) const
{
  Tensor2_3D out;
  for(int n=0; n<9; ++n)
  {
    const int loc = 9*n;
    out(n) = ten[loc] * a(0) + ten[loc+1] * a(1) + ten[loc+2] * a(2)
      + ten[loc+3] * a(3) + ten[loc+4] * a(4) + ten[loc+5] * a(5)
      + ten[loc+6] * a(6) + ten[loc+7] * a(7) + ten[loc+8] * a(8);
  }
  return out;
}

double Tensor4_3D::LnRContraction( const Tensor2_3D &Left, const Tensor2_3D &Right ) const
{
  double sum = 0.0;
  for(int n=0; n<9; ++n)
  {
    const int loc = 9*n;
    sum += ( ten[loc] * Right(0) + ten[loc+1] * Right(1) + ten[loc+2] * Right(2)
        + ten[loc+3] * Right(3) + ten[loc+4] * Right(4) + ten[loc+5] * Right(5)
        + ten[loc+6] * Right(6) + ten[loc+7] * Right(7) + ten[loc+8] * Right(8) ) * Left(n);
  }
  return sum;
}

double Tensor4_3D::Ten4Contraction( const Tensor4_3D &input ) const
{
  double sum = 0.0;
  for(int ii=0; ii<81; ++ii) sum += input(ii) * ten[ii];
  return sum;
}

// Ten[27I + 9J + 3K + L] = L[27I + 9J + 3M + N] * R[27M + 9N + 3K + L]
// Let ii = 3I + J ; jj = 3K + L ; kk = 3M + N
// then, Ten[9ii+jj] = L[9ii + kk] * R[9kk + jj] 
void Tensor4_3D::TenMult( const Tensor4_3D &tleft, const Tensor4_3D &tright )
{
  double temp[81];
  for(int ii=0; ii<9; ++ii)
  {
    for(int jj=0; jj<9; ++jj)
    {
      const int index = 9*ii + jj;
      temp[index] = 0.0;
      for(int kk=0; kk<9; ++kk)
        temp[index] += tleft(9*ii+kk) * tright(9*kk+jj);
    }
  }

  for(int ii=0; ii<81; ++ii) ten[ii] = temp[ii];
}

void Tensor4_3D::TenRMult( const Tensor4_3D &tright )
{
  double temp[81];
  for(int ii=0; ii<9; ++ii)
  {
    for(int jj=0; jj<9; ++jj)
    {
      const int index = 9*ii + jj;
      temp[index] = 0.0;
      for(int kk=0; kk<9; ++kk)
        temp[index] += ten[9*ii+kk] * tright(9*kk+jj);
    }
  }

  for(int ii=0; ii<81; ++ii) ten[ii] = temp[ii];
}

void Tensor4_3D::TenLMult( const Tensor4_3D &tleft )
{
  double temp[81];
  for(int ii=0; ii<9; ++ii)
  {
    for(int jj=0; jj<9; ++jj)
    {
      const int index = 9*ii + jj;
      temp[index] = 0.0;
      for(int kk=0; kk<9; ++kk)
        temp[index] += tleft(9*ii+kk) * ten[9*kk+jj];
    }
  }

  for(int ii=0; ii<81; ++ii) ten[ii] = temp[ii];
}

void Tensor4_3D::TenLRMult( const Tensor4_3D &tleft, const Tensor4_3D &tright )
{
  double temp[81];
  for(int ii=0; ii<9; ++ii)
  {
    for(int jj=0; jj<9; ++jj)
    {
      const int index = 9*ii+jj;
      temp[index] = 0.0;
      for(int mm=0; mm<9; ++mm)
      {
        for(int nn=0; nn<9; ++nn)
          temp[index] += tleft(9*ii+mm) * ten[9*mm+nn] * tright(9*nn+jj);
      }
    }
  }
  
  for(int ii=0; ii<81; ++ii) ten[ii] = temp[ii];
}

void Tensor4_3D::TenPMult( const Tensor4_3D &P )
{
  double temp[81];
  for(int ii=0; ii<9; ++ii)
  {
    for(int jj=0; jj<9; ++jj)
    {
      const int index = 9*ii + jj;
      temp[index] = 0.0;
      for(int mm=0; mm<9; ++mm)
      {
        for(int nn=0; nn<9; ++nn)
          temp[index] += P(9*ii+mm) * ten[9*mm+nn] * P(9*jj+nn);
      }
    }
  }
  
  for(int ii=0; ii<81; ++ii) ten[ii] = temp[ii];
}

bool Tensor4_3D::is_major_sym( const double &tol ) const
{
  for( int ii=0; ii<3; ++ii )
    for( int jj=0; jj<3; ++jj )
      for( int kk=0; kk<3; ++kk )
        for( int ll=0; ll<3; ++ll )
          if( std::abs( ten[ 27*ii+ 9*jj + 3*kk + ll ] - ten[ 27*kk+ 9*ll + 3*ii + jj ] ) >= tol ) return false;

  return true;
}

bool Tensor4_3D::is_minor_sym( const double &tol ) const
{
  for( int ii=0; ii<3; ++ii )
    for( int jj=0; jj<3; ++jj )
      for( int kk=0; kk<3; ++kk )
        for( int ll=0; ll<3; ++ll )
          if( std::abs( ten[ 27*ii+ 9*jj + 3*kk + ll ] - ten[ 27*jj+ 9*ii + 3*kk + ll ] ) >= tol
            || std::abs( ten[ 27*ii+ 9*jj + 3*kk + ll ] - ten[ 27*ii+ 9*jj + 3*ll + kk ] ) >= tol ) return false;

  return true;
}

Tensor2_3D Tensor4_3D::solve( const Tensor2_3D &B )
{
  ASSERT( is_minor_sym() == true, "The rank-four tensor needs to satisfy the minor symmetry.\n" );  

  MATH_T::Matrix_Dense<6> AA_mat;

  AA_mat(0, 0) = ten[0];
  AA_mat(0, 5) = ten[1];
  AA_mat(0, 4) = ten[2];
  AA_mat(0, 1) = ten[4];
  AA_mat(0, 3) = ten[5];
  AA_mat(0, 2) = ten[8];

  AA_mat(5, 0) = ten[9];
  AA_mat(5, 5) = ten[10];
  AA_mat(5, 4) = ten[11];
  AA_mat(5, 1) = ten[13];
  AA_mat(5, 3) = ten[14];
  AA_mat(5, 2) = ten[17];

  AA_mat(4, 0) = ten[18];
  AA_mat(4, 5) = ten[19];
  AA_mat(4, 4) = ten[20];
  AA_mat(4, 1) = ten[22];
  AA_mat(4, 3) = ten[23];
  AA_mat(4, 2) = ten[26];

  AA_mat(1, 0) = ten[36];
  AA_mat(1, 5) = ten[37];
  AA_mat(1, 4) = ten[38];
  AA_mat(1, 1) = ten[40];
  AA_mat(1, 3) = ten[41];
  AA_mat(1, 2) = ten[44];

  AA_mat(3, 0) = ten[45];
  AA_mat(3, 5) = ten[46];
  AA_mat(3, 4) = ten[47];
  AA_mat(3, 1) = ten[49];
  AA_mat(3, 3) = ten[50];
  AA_mat(3, 2) = ten[53];

  AA_mat(2, 0) = ten[72];
  AA_mat(2, 5) = ten[73];
  AA_mat(2, 4) = ten[74];
  AA_mat(2, 1) = ten[76];
  AA_mat(2, 3) = ten[77];
  AA_mat(2, 2) = ten[80];

  AA_mat.LU_fac();

  const auto out_array = AA_mat.LU_solve(
    {{ B(0),
       B(4),
       B(8),
       B(5),
       B(2),
       B(1) }} );

  return Tensor2_3D( 
          out_array[0], 
      0.5*out_array[5], 
      0.5*out_array[4],
      0.5*out_array[5], 
          out_array[1],
      0.5*out_array[3],
      0.5*out_array[4],
      0.5*out_array[3],
          out_array[2] );
}

Tensor4_3D Tensor4_3D::solve( const Tensor4_3D &BB )
{
  Tensor4_3D out = Ten4::gen_zero();

  for( int cc=0; cc<9; ++cc )
  {
    const Tensor2_3D B( BB( cc ), BB( 9  + cc ), BB( 18 + cc ),
     BB( 27 + cc ), BB( 36 + cc ), BB( 45 + cc ),
     BB( 54 + cc ), BB( 63 + cc ), BB( 72 + cc ) );

    const Tensor2_3D temp = solve( B );

    out(      cc ) = temp( 0 );
    out( 9  + cc ) = temp( 1 );
    out( 18 + cc ) = temp( 2 );

    out( 27 + cc ) = temp( 3 );
    out( 36 + cc ) = temp( 4 );
    out( 45 + cc ) = temp( 5 );

    out( 54 + cc ) = temp( 6 );
    out( 63 + cc ) = temp( 7 );
    out( 72 + cc ) = temp( 8 );
  }

  return out;
}

Tensor2_3D operator*( const Tensor4_3D &input, const Tensor2_3D &aa )
{
  Tensor2_3D out;
  for(int n=0; n<9; ++n)
  {
    const int loc = 9*n;
    out(n) = input(loc) * aa(0) + input(loc+1) * aa(1) + input(loc+2) * aa(2)
      + input(loc+3) * aa(3) + input(loc+4) * aa(4) + input(loc+5) * aa(5)
      + input(loc+6) * aa(6) + input(loc+7) * aa(7) + input(loc+8) * aa(8);
  }
  
  return out;
}

Tensor2_3D operator*( const Tensor2_3D &aa, const Tensor4_3D &input )
{
  Tensor2_3D out;
  for(int m=0; m<9; ++m)
  {
    out(m) = aa(0) * input(m) + aa(3) * input(m+27) + aa(6) * input(m+54)
      + aa(1) * input(m+9) + aa(4) * input(m+36) + aa(7) * input(m+63)
      + aa(2) * input(m+18) + aa(5) * input(m+45) + aa(8) * input(m+72);
  }
  return out;
}

Tensor4_3D operator*( const Tensor4_3D &tleft, const Tensor4_3D &tright )
{
  Tensor4_3D out;
  for(int ii=0; ii<9; ++ii)
  {
    for(int jj=0; jj<9; ++jj)
    {
      const int index = 9*ii + jj;
      out(index) = 0.0;
      for(int kk=0; kk<9; ++kk)
        out(index) += tleft(9*ii+kk) * tright(9*kk+jj);
    }
  }

  return out;
}

Tensor4_3D operator*( const double &val, const Tensor4_3D &input )
{
  Tensor4_3D out;
  for(int ii=0; ii<81; ++ii) out(ii) = val * input(ii);

  return out;
}

Tensor4_3D Ten4::gen_symm_id()
{
  constexpr std::array<double,81> temp {{ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 }};
  return Tensor4_3D(temp);
}

Tensor4_3D Ten4::gen_P( const Tensor2_3D &C, const Tensor2_3D &invC )
{
  Tensor4_3D out = Ten4::gen_symm_id();
  
  out.add_OutProduct( -1.0 / 3.0, invC, C );

  return out;
}

Tensor4_3D Ten4::gen_P( const Tensor2_3D &C )
{
  return Ten4::gen_P( C, Ten2::inverse(C) );
}

Tensor4_3D Ten4::gen_P( const SymmTensor2_3D &C, const SymmTensor2_3D &invC )
{
  Tensor4_3D out = Ten4::gen_symm_id();
  
  out.add_OutProduct( -1.0 / 3.0, invC, C );

  return out;
}

Tensor4_3D Ten4::gen_P( const SymmTensor2_3D &C )
{
  return Ten4::gen_P( C, STen2::inverse(C) );
}


Tensor4_3D Ten4::gen_Pt( const Tensor2_3D &C )
{
  return Ten4::gen_P( Ten2::inverse(C) );
}

Tensor4_3D Ten4::gen_Pt( const SymmTensor2_3D &C )
{
  return Ten4::gen_P( STen2::inverse(C) );
}

Tensor4_3D Ten4::gen_Ptilde( const Tensor2_3D &invC )
{
  Tensor4_3D out;
  out.gen_zero();
  out.add_SymmProduct(1.0, invC, invC);
  out.add_OutProduct( -1.0 / 3.0, invC, invC );
  return out;
}

Tensor4_3D Ten4::transpose( const Tensor4_3D &input )
{
  std::array<double,81> temp;

  for(int mm=0; mm<9; ++mm)
  {
    temp[   mm] = input(9*mm);
    temp[ 9+mm] = input(9*mm+1);
    temp[18+mm] = input(9*mm+2);
    temp[27+mm] = input(9*mm+3);
    temp[36+mm] = input(9*mm+4);
    temp[45+mm] = input(9*mm+5);
    temp[54+mm] = input(9*mm+6);
    temp[63+mm] = input(9*mm+7);
    temp[72+mm] = input(9*mm+8);
  }

  return Tensor4_3D(temp);
}

// EOF
