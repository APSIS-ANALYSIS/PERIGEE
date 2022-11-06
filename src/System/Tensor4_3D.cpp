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

Tensor4_3D::~Tensor4_3D()
{}

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
  std::cout<<"Tensor4_3D: \n\n";
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

void Tensor4_3D::gen_P( const Matrix_3x3 &C, const Matrix_3x3 &invC )
{
  gen_symm_id();
  add_OutProduct( -1.0 / 3.0, invC, C );
}

void Tensor4_3D::gen_Ptilde( const Matrix_3x3 &invC )
{
  gen_zero();
  add_SymmProduct(1.0, invC, invC);
  add_OutProduct( -1.0 / 3.0, invC, invC );
}

void Tensor4_3D::gen_rand()
{
  srand(time(NULL));

  for(int ii=0; ii<81; ++ii)
  {
    double value = rand() % 100000;

    ten[ii] = value * 1.0e-4 - 5.0; // range [-5, 4.9999]
  }
}

void Tensor4_3D::gen_zero()
{
  for(int ii=0; ii<81; ++ii) ten[ii] = 0.0;
}

void Tensor4_3D::scale( const double &val )
{
  for(int ii=0; ii<81; ++ii) ten[ii] *= val;
}

void Tensor4_3D::PY( const Tensor4_3D &input )
{
  for(int ii=0; ii<81; ++ii) ten[ii] += input(ii);
}

void Tensor4_3D::AXPY( const double &val, const Tensor4_3D &input )
{
  for(int ii=0; ii<81; ++ii) ten[ii] += val * input(ii);
}

void Tensor4_3D::add_OutProduct( const double &val, const Matrix_3x3 &mleft,
    const Matrix_3x3 &mright )
{
  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      for(int kk=0; kk<3; ++kk)
      {
        for(int ll=0; ll<3; ++ll)
          ten[27*ii+9*jj+3*kk+ll] += val * mleft(3*ii+jj) * mright(3*kk+ll);
      }
    }
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

void Tensor4_3D::add_SymmProduct( const double &val, const Matrix_3x3 &mleft,
    const Matrix_3x3 &mright )
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

void Tensor4_3D::add_SymmOutProduct( const double &val, const Matrix_3x3 &mleft,
    const Matrix_3x3 &mright )
{
  for(int ii=0; ii<3 ; ++ii)
  {
    for(int jj=0; jj<3 ; ++jj)
    {
      for(int kk=0; kk<3; ++kk)
      {
        for(int ll=0; ll<3; ++ll)
        {
          ten[27*ii + 9*jj + 3*kk + ll] += val *  (mleft(3*ii + jj) * mright(3*jj + kk) 
              + mright(3*ii + jj) * mleft(3*kk + ll));
        }
      }
    }
  }
}

void Tensor4_3D::MatMult_1( const Matrix_3x3 &source )
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

void Tensor4_3D::MatMult_2( const Matrix_3x3 &source )
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

void Tensor4_3D::MatMult_3( const Matrix_3x3 &source )
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

void Tensor4_3D::MatMult_4( const Matrix_3x3 &source )
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

void Tensor4_3D::LeftContraction( const Matrix_3x3 &a, Matrix_3x3 &out ) const
{
  for(int m=0; m<9; ++m)
  {
    out(m) = a(0) * ten[m] + a(3) * ten[m+27] + a(6) * ten[m+54] 
      + a(1) * ten[m+9] + a(4) * ten[m+36] + a(7) * ten[m+63]
      + a(2) * ten[m+18] + a(5) * ten[m+45] + a(8) * ten[m+72];
  }
}

void Tensor4_3D::RightContraction( const Matrix_3x3 &a, Matrix_3x3 &out ) const
{
  for(int n=0; n<9; ++n)
  {
    const int loc = 9*n;
    out(n) = ten[loc] * a(0) + ten[loc+1] * a(1) + ten[loc+2] * a(2)
      + ten[loc+3] * a(3) + ten[loc+4] * a(4) + ten[loc+5] * a(5)
      + ten[loc+6] * a(6) + ten[loc+7] * a(7) + ten[loc+8] * a(8);
  }
}

double Tensor4_3D::LnRContraction( const Matrix_3x3 &Left, const Matrix_3x3 &Right ) const
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

// EOF
