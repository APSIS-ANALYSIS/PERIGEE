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
  double add_ten[21];
  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      for(int kk=0; kk<3; ++kk)
      {
        for(int ll=0; ll<3; ++ll)
        {
          add_ten[ Voigt_notation(ii,jj,kk,ll) ] = val * 0.5 * ( mleft( Voigt_notation(ii,kk) ) 
            * mright( Voigt_notation(jj,ll) ) + mleft( Voigt_notation(ii,ll) ) * mright( Voigt_notation(jj,kk) ) );
        }
      }
    }
  }
  for(int counter=0; counter<21; ++counter)
  {
    ten[counter] += add_ten[counter];
  }
}

void SymmTensor4_3D::add_SymmOutProduct( const double &val, const SymmMatrix_3x3 &mleft,
    const SymmMatrix_3x3 &mright )
{
  double add_ten[21];
  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      for(int kk=0; kk<3; ++kk)
      {
        for(int ll=0; ll<3; ++ll)
        {
          add_ten[ Voigt_notation(ii,jj,kk,ll) ] = val * ( mleft( Voigt_notation(ii,jj) ) 
            * mright( Voigt_notation(kk,ll) ) + mleft( Voigt_notation(kk,ll) ) * mright( Voigt_notation(ii,jj) ) );
        }
      }
    }
  }
  for(int counter=0; counter<21; ++counter)
  {
    ten[counter] += add_ten[counter];
  }
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
