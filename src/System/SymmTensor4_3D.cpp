#include "SymmTensor4_3D.hpp"

SymmTensor4_3D::SymmTensor4_3D()
{
  gen_symm_id();
}

SymmTensor4_3D::~SymmTensor4_3D()
{}

bool SymmTensor4_3D::is_identical(const Tensor4_3D &source, const double &tol) const
{
  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      for(int kk=0; kk<3; ++kk)
      {
        for(int ll=0; ll<3; ++ll)
        {
          if( std::abs(source(27*ii+9*jj+3*kk+ll) - ten[Voigt_notation(ii,jj,kk,ll)] ) > tol )
          {
            return false;
          }
        }
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

void SymmTensor4_3D::gen_rand()
{
  srand(time(NULL));

  for(int ii=0; ii<21; ++ii)
  {
    double value = rand() % 100000;

    ten[ii] = value * 1.0e-4 - 5.0; // range [-5, 4.9999]
  }
}

void SymmTensor4_3D::gen_zero()
{
  ten = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,0};
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

void SymmTensor4_3D::add_OutProduct( const double &val, const SymmMatrix_3x3 &mmat )
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

void SymmTensor4_3D::add_SymmOutProduct( const double &val, const Vector_3 &vec1, 
  const Vector_3 &vec2 )
{
  ten[0]  += val * 4 * vec1(0) * vec2(0) * vec1(0) * vec2(0);

  ten[1]  += val * 4 * vec1(0) * vec2(0) * vec1(1) * vec2(1);

  ten[2]  += val * 4 * vec1(0) * vec2(0) * vec1(2) * vec2(2);

  ten[3]  += val * ( 2 * vec1(0) * vec2(0) * vec1(1) * vec2(2)
                   + 2 * vec1(0) * vec2(0) * vec1(2) * vec2(1) );

  ten[4]  += val * ( 2 * vec1(0) * vec2(0) * vec1(0) * vec2(2)
                   + 2 * vec1(0) * vec2(0) * vec1(2) * vec2(0) );

  ten[5]  += val * ( 2 * vec1(0) * vec2(0) * vec1(0) * vec2(1)
                   + 2 * vec1(0) * vec2(0) * vec1(1) * vec2(0) );

  ten[6]  += val * 4 * vec1(1) * vec2(1) * vec1(1) * vec2(1);

  ten[7]  += val * 4 * vec1(1) * vec2(1) * vec1(2) * vec2(2);

  ten[8]  += val * ( 2 * vec1(1) * vec2(1) * vec1(1) * vec2(2)
                   + 2 * vec1(1) * vec2(1) * vec1(2) * vec2(1) );

  ten[9]  += val * ( 2 * vec1(0) * vec2(2) * vec1(1) * vec2(1)
                   + 2 * vec1(2) * vec2(0) * vec1(1) * vec2(1) );

  ten[10] += val * ( 2 * vec1(1) * vec2(0) * vec1(1) * vec2(1)
                   + 2 * vec1(0) * vec2(1) * vec1(1) * vec2(1) );

  ten[11] += val * 4 * vec1(2) * vec2(2) * vec1(2) * vec2(2);

  ten[12] += val * ( 2 * vec1(1) * vec2(2) * vec1(2) * vec2(2)
                   + 2 * vec1(2) * vec2(1) * vec1(2) * vec2(2) );

  ten[13] += val * ( 2 * vec1(2) * vec2(0) * vec1(2) * vec2(2)
                   + 2 * vec1(0) * vec2(2) * vec1(2) * vec2(2) );

  ten[14] += val * ( 2 * vec1(0) * vec2(1) * vec1(2) * vec2(2)
                   + 2 * vec1(1) * vec2(0) * vec1(2) * vec2(2) );

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

void SymmTensor4_3D::add_SymmProduct( const double &val, const SymmMatrix_3x3 &mleft,
  const SymmMatrix_3x3 &mright )
{
  ten[0]  += val * 0.5 * 2.0 * ( mleft(0) * mright(0) );
  ten[1]  += val * 0.5 * 2.0 * ( mleft(5) * mright(5) );
  ten[2]  += val * 0.5 * 2.0 * ( mleft(4) * mright(4) );

  ten[3]  += val * 0.5 * ( mleft(5) * mright(4) + mleft(4) * mright(5) );
  ten[4]  += val * 0.5 * ( mleft(0) * mright(4) + mleft(4) * mright(0) );
  ten[5]  += val * 0.5 * ( mleft(0) * mright(5) + mleft(5) * mright(0) );

  ten[6]  += val * 0.5 * 2.0 * ( mleft(1) * mright(1) );
  ten[7]  += val * 0.5 * 2.0 * ( mleft(3) * mright(3) );

  ten[8]  += val * 0.5 * ( mleft(1) * mright(3) + mleft(3) * mright(1) );
  ten[9]  += val * 0.5 * ( mleft(5) * mright(3) + mleft(5) * mright(3) );
  ten[10] += val * 0.5 * ( mleft(5) * mright(1) + mleft(5) * mright(1) );

  ten[11] += val * 0.5 * 2.0 * ( mleft(2) * mright(2) );

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

void SymmTensor4_3D::add_SymmOutProduct( const double &val, const SymmMatrix_3x3 &mleft,
    const SymmMatrix_3x3 &mright )
{
  ten[0]  = val * ( mleft(0) * mright(0) + mleft(0) * mright(0) );
  ten[1]  = val * ( mleft(0) * mright(1) + mleft(1) * mright(0) );
  ten[2]  = val * ( mleft(0) * mright(2) + mleft(2) * mright(0) );
  ten[3]  = val * ( mleft(0) * mright(3) + mleft(3) * mright(0) );
  ten[4]  = val * ( mleft(0) * mright(4) + mleft(4) * mright(0) );
  ten[5]  = val * ( mleft(0) * mright(5) + mleft(5) * mright(0) );
  ten[6]  = val * ( mleft(1) * mright(1) + mleft(1) * mright(1) );
  ten[7]  = val * ( mleft(1) * mright(2) + mleft(2) * mright(1) );
  ten[8]  = val * ( mleft(1) * mright(3) + mleft(3) * mright(1) );
  ten[9]  = val * ( mleft(4) * mright(1) + mleft(1) * mright(4) );
  ten[10] = val * ( mleft(5) * mright(1) + mleft(1) * mright(5) );
  ten[11] = val * ( mleft(2) * mright(2) + mleft(2) * mright(2) );
  ten[12] = val * ( mleft(3) * mright(2) + mleft(2) * mright(3) );
  ten[13] = val * ( mleft(4) * mright(2) + mleft(2) * mright(4) );
  ten[14] = val * ( mleft(5) * mright(2) + mleft(2) * mright(5) );
  ten[15] = val * ( mleft(3) * mright(3) + mleft(3) * mright(3) );
  ten[16] = val * ( mleft(4) * mright(3) + mleft(3) * mright(4) );
  ten[17] = val * ( mleft(5) * mright(3) + mleft(3) * mright(5) );
  ten[18] = val * ( mleft(4) * mright(4) + mleft(4) * mright(4) );
  ten[19] = val * ( mleft(5) * mright(4) + mleft(4) * mright(5) );
  ten[20] = val * ( mleft(5) * mright(5) + mleft(5) * mright(5) );
}

int SymmTensor4_3D::Voigt_notation( const int &ii, const int &jj, const int &kk, const int &ll ) const
{
  const int index_I = Voigt_notation(ii, jj);
  const int index_J = Voigt_notation(kk, ll);

  const int mat[36] = {0, 1,  2,  3,  4,  5,
                       1, 6,  7,  8,  9,  10,
                       2, 7,  11, 12, 13, 14,
                       3, 8,  12, 15, 16, 17,
                       4, 9,  13, 16, 18, 19,
                       5, 10, 14, 17, 19, 20 };

  return mat[ 6 * index_I + index_J ];
}

int SymmTensor4_3D::Voigt_notation( const int &ii, const int &jj ) const
{
  const int mat[9] = { 0, 5, 4, 
                       5, 1, 3, 
                       4, 3, 2 };
  return mat[ 3 * ii + jj ];
}
// EOF
