#include "SymmTensor4_3D.hpp"

constexpr std::array<int,9> SymmTensor4_3D::map;

constexpr std::array<int,36> SymmTensor4_3D::mapper;

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

SymmTensor4_3D& SymmTensor4_3D::operator= (const SymmTensor4_3D &source)
{
  if(this != &source) ten = source.ten;

  return *this;
}

void SymmTensor4_3D::print(std::ostream& os, const std::string& delimiter) const
{
  os<<"SymmTensor4_3D: \n";
  for(int kk=0; kk<3; ++kk)
  {
    for(int ll=0; ll<3; ++ll)
    {
      os<<"k = "<<kk<<"\tl = "<<ll<<'\n';
      for(int ii=0; ii<3; ++ii)
      {
        for(int jj=0; jj<3; ++jj)
        {
          os<<"i = "<<ii<<delimiter<<"j = "<<jj<<delimiter;
          os<<std::setprecision(6)<<ten[ Voigt_notation(ii, jj, kk, ll) ]<<delimiter;
        }
        os<<'\n';
      }
      os<<'\n';
    }
  }  
}

void SymmTensor4_3D::print_in_mat() const
{
  std::cout<<"SymmTensor4_3D:\n";
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
  for(int ii=0; ii<21; ++ii) ten[ii] += source.ten[ii];
  return *this;
}

SymmTensor4_3D& SymmTensor4_3D::operator-=( const SymmTensor4_3D &source )
{
  for(int ii=0; ii<21; ++ii) ten[ii] -= source.ten[ii];
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

void SymmTensor4_3D::TenQMult( const SymmTensor4_3D &QQ )
{
  std::array<double, 21> temp = {{0.0}}; 

  const int map0 = 6*map[0];
  const int map1 = 6*map[1];
  const int map2 = 6*map[2];
  const int map4 = 6*map[4];
  const int map5 = 6*map[5];
  const int map8 = 6*map[8];

  for(int ii=0; ii<9; ++ii)
  {
    for(int jj=0; jj<9; ++jj)
    {
      const int index = 6*map[ii]+map[jj];

      temp[0]  += QQ(mapper[map0+map[ii]]) * ten[mapper[index]] * QQ(mapper[map0+map[jj]]);
      temp[1]  += QQ(mapper[map0+map[ii]]) * ten[mapper[index]] * QQ(mapper[map4+map[jj]]);
      temp[2]  += QQ(mapper[map0+map[ii]]) * ten[mapper[index]] * QQ(mapper[map8+map[jj]]);
      temp[3]  += QQ(mapper[map0+map[ii]]) * ten[mapper[index]] * QQ(mapper[map5+map[jj]]);
      temp[4]  += QQ(mapper[map0+map[ii]]) * ten[mapper[index]] * QQ(mapper[map2+map[jj]]);
      temp[5]  += QQ(mapper[map0+map[ii]]) * ten[mapper[index]] * QQ(mapper[map1+map[jj]]);
      temp[6]  += QQ(mapper[map4+map[ii]]) * ten[mapper[index]] * QQ(mapper[map4+map[jj]]);
      temp[7]  += QQ(mapper[map4+map[ii]]) * ten[mapper[index]] * QQ(mapper[map8+map[jj]]);
      temp[8]  += QQ(mapper[map4+map[ii]]) * ten[mapper[index]] * QQ(mapper[map5+map[jj]]);
      temp[9]  += QQ(mapper[map2+map[ii]]) * ten[mapper[index]] * QQ(mapper[map4+map[jj]]);
      temp[10] += QQ(mapper[map1+map[ii]]) * ten[mapper[index]] * QQ(mapper[map4+map[jj]]);
      temp[11] += QQ(mapper[map8+map[ii]]) * ten[mapper[index]] * QQ(mapper[map8+map[jj]]);
      temp[12] += QQ(mapper[map5+map[ii]]) * ten[mapper[index]] * QQ(mapper[map8+map[jj]]);
      temp[13] += QQ(mapper[map2+map[ii]]) * ten[mapper[index]] * QQ(mapper[map8+map[jj]]);
      temp[14] += QQ(mapper[map1+map[ii]]) * ten[mapper[index]] * QQ(mapper[map8+map[jj]]);
      temp[15] += QQ(mapper[map5+map[ii]]) * ten[mapper[index]] * QQ(mapper[map5+map[jj]]);
      temp[16] += QQ(mapper[map2+map[ii]]) * ten[mapper[index]] * QQ(mapper[map5+map[jj]]);
      temp[17] += QQ(mapper[map1+map[ii]]) * ten[mapper[index]] * QQ(mapper[map5+map[jj]]);
      temp[18] += QQ(mapper[map2+map[ii]]) * ten[mapper[index]] * QQ(mapper[map2+map[jj]]);
      temp[19] += QQ(mapper[map1+map[ii]]) * ten[mapper[index]] * QQ(mapper[map2+map[jj]]);
      temp[20] += QQ(mapper[map1+map[ii]]) * ten[mapper[index]] * QQ(mapper[map1+map[jj]]);
   }
  }
  ten = temp;
}

void SymmTensor4_3D::TenPMult( const SymmTensor2_3D &C )
{
  const auto invC  = STen2::inverse(C);
  
  std::array<double, 21> temp = {{0.0}}; 

  constexpr double pt33 = 1.0 / 3.0; 
  constexpr double pt11 = 1.0 / 9.0; 

  temp = ten;

  for (int mm = 0; mm < 3; ++mm)
  {
    for (int nn = 0; nn < 3; ++nn)
    {
      const int index_1 = map[3*mm + nn];
      // mapper index for IJMN
      const int map0_Ind1 = 6*map[ 0 ] + index_1; 
      const int map4_Ind1 = 6*map[ 4 ] + index_1;
      const int map8_Ind1 = 6*map[ 8 ] + index_1;
      const int map5_Ind1 = 6*map[ 5 ] + index_1;
      const int map2_Ind1 = 6*map[ 2 ] + index_1;
      const int map1_Ind1 = 6*map[ 1 ] + index_1;
      // mapper index for MNKL
      const int Ind1_map0 = 6*index_1 + map[ 0 ];
      const int Ind1_map4 = 6*index_1 + map[ 4 ];
      const int Ind1_map8 = 6*index_1 + map[ 8 ];
      const int Ind1_map5 = 6*index_1 + map[ 5 ];
      const int Ind1_map2 = 6*index_1 + map[ 2 ];
      const int Ind1_map1 = 6*index_1 + map[ 1 ];

      temp[0] -= pt33 * ten[ mapper[ map0_Ind1 ] ] * C(index_1) * invC(map[ 0 ]) + pt33 * invC(map[ 0 ]) * C(index_1) * ten[ mapper[ Ind1_map0 ] ];

      temp[1] -= pt33 * ten[ mapper[ map0_Ind1 ] ] * C(index_1) * invC(map[ 4 ]) + pt33 * invC(map[ 0 ]) * C(index_1) * ten[ mapper[ Ind1_map4 ] ];

      temp[2] -= pt33 * ten[ mapper[ map0_Ind1 ] ] * C(index_1) * invC(map[ 8 ]) + pt33 * invC(map[ 0 ]) * C(index_1) * ten[ mapper[ Ind1_map8 ] ];

      temp[3] -= pt33 * ten[ mapper[ map0_Ind1 ] ] * C(index_1) * invC(map[ 5 ]) + pt33 * invC(map[ 0 ]) * C(index_1) * ten[ mapper[ Ind1_map5 ] ];

      temp[4] -= pt33 * ten[ mapper[ map0_Ind1 ] ] * C(index_1) * invC(map[ 2 ]) + pt33 * invC(map[ 0 ]) * C(index_1) * ten[ mapper[ Ind1_map2 ] ];

      temp[5] -= pt33 * ten[ mapper[ map0_Ind1 ] ] * C(index_1) * invC(map[ 1 ]) + pt33 * invC(map[ 0 ]) * C(index_1) * ten[ mapper[ Ind1_map1 ] ];

      temp[6] -= pt33 * ten[ mapper[ map4_Ind1 ] ] * C(index_1) * invC(map[ 4 ]) + pt33 * invC(map[ 4 ]) * C(index_1) * ten[ mapper[ Ind1_map4 ] ];

      temp[7] -= pt33 * ten[ mapper[ map4_Ind1 ] ] * C(index_1) * invC(map[ 8 ]) + pt33 * invC(map[ 4 ]) * C(index_1) * ten[ mapper[ Ind1_map8 ] ];

      temp[8] -= pt33 * ten[ mapper[ map4_Ind1 ] ] * C(index_1) * invC(map[ 5 ]) + pt33 * invC(map[ 4 ]) * C(index_1) * ten[ mapper[ Ind1_map5 ] ];

      temp[9] -= pt33 * ten[ mapper[ map4_Ind1 ] ] * C(index_1) * invC(map[ 2 ]) + pt33 * invC(map[ 4 ]) * C(index_1) * ten[ mapper[ Ind1_map2 ] ];

      temp[10] -= pt33 * ten[ mapper[ map4_Ind1 ] ] * C(index_1) * invC(map[ 1 ]) + pt33 * invC(map[ 4 ]) * C(index_1) * ten[ mapper[ Ind1_map1 ] ];

      temp[11] -= pt33 * ten[ mapper[ map8_Ind1 ] ] * C(index_1) * invC(map[ 8 ]) + pt33 * invC(map[ 8 ]) * C(index_1) * ten[ mapper[ Ind1_map8 ] ];

      temp[12] -= pt33 * ten[ mapper[ map8_Ind1 ] ] * C(index_1) * invC(map[ 5 ]) + pt33 * invC(map[ 8 ]) * C(index_1) * ten[ mapper[ Ind1_map5 ] ];

      temp[13] -= pt33 * ten[ mapper[ map8_Ind1 ] ] * C(index_1) * invC(map[ 2 ]) + pt33 * invC(map[ 8 ]) * C(index_1) * ten[ mapper[ Ind1_map2 ] ];

      temp[14] -= pt33 * ten[ mapper[ map8_Ind1 ] ] * C(index_1) * invC(map[ 1 ]) + pt33 * invC(map[ 8 ]) * C(index_1) * ten[ mapper[ Ind1_map1 ] ];

      temp[15] -= pt33 * ten[ mapper[ map5_Ind1 ] ] * C(index_1) * invC(map[ 5 ]) + pt33 * invC(map[ 5 ]) * C(index_1) * ten[ mapper[ Ind1_map5 ] ];

      temp[16] -= pt33 * ten[ mapper[ map5_Ind1 ] ] * C(index_1) * invC(map[ 2 ]) + pt33 * invC(map[ 5 ]) * C(index_1) * ten[ mapper[ Ind1_map2 ] ];

      temp[17] -= pt33 * ten[ mapper[ map5_Ind1 ] ] * C(index_1) * invC(map[ 1 ]) + pt33 * invC(map[ 5 ]) * C(index_1) * ten[ mapper[ Ind1_map1 ] ];

      temp[18] -= pt33 * ten[ mapper[ map2_Ind1 ] ] * C(index_1) * invC(map[ 2 ]) + pt33 * invC(map[ 2 ]) * C(index_1) * ten[ mapper[ Ind1_map2 ] ];

      temp[19] -= pt33 * ten[ mapper[ map2_Ind1 ] ] * C(index_1) * invC(map[ 1 ]) + pt33 * invC(map[ 2 ]) * C(index_1) * ten[ mapper[ Ind1_map1 ] ];

      temp[20] -= pt33 * ten[ mapper[ map1_Ind1 ] ] * C(index_1) * invC(map[ 1 ]) + pt33 * invC(map[ 1 ]) * C(index_1) * ten[ mapper[ Ind1_map1 ] ];
      
      for (int pp = 0; pp < 3; ++pp)
      {
        for (int qq = 0; qq < 3; ++qq)
        {
          const int index_2 = map[3*pp + qq];
          // mapper index for MNPQ
          const int Ind1_Ind2 = 6*index_1 + index_2; 

          temp[0] += pt11 * invC(map[ 0 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 0 ]);

          temp[1] += pt11 * invC(map[ 0 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 4 ]);

          temp[2] += pt11 * invC(map[ 0 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 8 ]);

          temp[3] += pt11 * invC(map[ 0 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 5 ]);

          temp[4] += pt11 * invC(map[ 0 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 2 ]);

          temp[5] += pt11 * invC(map[ 0 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 1 ]);

          temp[6] += pt11 * invC(map[ 4 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 4 ]);

          temp[7] += pt11 * invC(map[ 4 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 8 ]);

          temp[8] += pt11 * invC(map[ 4 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 5 ]);

          temp[9] += pt11 * invC(map[ 4 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 2 ]);

          temp[10] += pt11 * invC(map[ 4 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 1 ]);

          temp[11] += pt11 * invC(map[ 8 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 8 ]);

          temp[12] += pt11 * invC(map[ 8 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 5 ]);

          temp[13] += pt11 * invC(map[ 8 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 2 ]);

          temp[14] += pt11 * invC(map[ 8 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 1 ]);

          temp[15] += pt11 * invC(map[ 5 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 5 ]);

          temp[16] += pt11 * invC(map[ 5 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 2 ]);

          temp[17] += pt11 * invC(map[ 5 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 1 ]);

          temp[18] += pt11 * invC(map[ 2 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 2 ]);

          temp[19] += pt11 * invC(map[ 2 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 1 ]);

          temp[20] += pt11 * invC(map[ 1 ]) * C(index_1) * ten[ mapper[ Ind1_Ind2 ] ] * C(index_2) * invC(map[ 1 ]);
        }
      }
    }
  }
  ten = temp;
}

void SymmTensor4_3D::pull_back_stiffness( const Tensor2_3D &invF )
{
  auto ten4 = full();

  ten4.MatMult_1( invF );
  ten4.MatMult_2( invF );
  ten4.MatMult_3( invF );
  ten4.MatMult_4( invF );

  ten = STen4::gen_symm_part(ten4).to_std_array();
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

SymmTensor2_3D operator*( const SymmTensor2_3D &mleft, const SymmTensor4_3D &mright )
{
  SymmTensor2_3D out = STen2::gen_zero();
  out(0) += mleft(0) * mright(0);
  out(0) += mleft(5) * mright(5);
  out(0) += mleft(4) * mright(4);
  out(0) += mleft(5) * mright(5);
  out(0) += mleft(1) * mright(1);
  out(0) += mleft(3) * mright(3);
  out(0) += mleft(4) * mright(4);
  out(0) += mleft(3) * mright(3);
  out(0) += mleft(2) * mright(2);

  out(1) += mleft(0) * mright(1);
  out(1) += mleft(5) * mright(10);
  out(1) += mleft(4) * mright(9);
  out(1) += mleft(5) * mright(10);
  out(1) += mleft(1) * mright(6);
  out(1) += mleft(3) * mright(8);
  out(1) += mleft(4) * mright(9);
  out(1) += mleft(3) * mright(8);
  out(1) += mleft(2) * mright(7);

  out(2) += mleft(0) * mright(2);
  out(2) += mleft(5) * mright(14);
  out(2) += mleft(4) * mright(13);
  out(2) += mleft(5) * mright(14);
  out(2) += mleft(1) * mright(7);
  out(2) += mleft(3) * mright(12);
  out(2) += mleft(4) * mright(13);
  out(2) += mleft(3) * mright(12);
  out(2) += mleft(2) * mright(11);

  out(3) += mleft(0) * mright(3);
  out(3) += mleft(5) * mright(17);
  out(3) += mleft(4) * mright(16);
  out(3) += mleft(5) * mright(17);
  out(3) += mleft(1) * mright(8);
  out(3) += mleft(3) * mright(15);
  out(3) += mleft(4) * mright(16);
  out(3) += mleft(3) * mright(15);
  out(3) += mleft(2) * mright(12);

  out(4) += mleft(0) * mright(4);
  out(4) += mleft(5) * mright(19);
  out(4) += mleft(4) * mright(18);
  out(4) += mleft(5) * mright(19);
  out(4) += mleft(1) * mright(9);
  out(4) += mleft(3) * mright(16);
  out(4) += mleft(4) * mright(18);
  out(4) += mleft(3) * mright(16);
  out(4) += mleft(2) * mright(13);

  out(5) += mleft(0) * mright(5);
  out(5) += mleft(5) * mright(20);
  out(5) += mleft(4) * mright(19);
  out(5) += mleft(5) * mright(20);
  out(5) += mleft(1) * mright(10);
  out(5) += mleft(3) * mright(17);
  out(5) += mleft(4) * mright(19);
  out(5) += mleft(3) * mright(17);
  out(5) += mleft(2) * mright(14);

  return out;
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

SymmTensor4_3D STen4::gen_symm_part( const Tensor4_3D &input )
{
  SymmTensor4_3D out = STen4::gen_zero();
  out(0)  += 0.125 * ( input(0)  + input(0)  + input(0)  + input(0)  + input(0)  + input(0)  + input(0)  + input(0) );
  out(1)  += 0.125 * ( input(4)  + input(4)  + input(4)  + input(4)  + input(36) + input(36) + input(36) + input(36) );
  out(2)  += 0.125 * ( input(72) + input(72) + input(72) + input(72) + input(8)  + input(8)  + input(8)  + input(8) );
  out(3)  += 0.125 * ( input(5)  + input(5)  + input(7)  + input(7)  + input(45) + input(63) + input(45) + input(63) );
  out(4)  += 0.125 * ( input(2)  + input(2)  + input(6)  + input(6)  + input(18) + input(54) + input(18) + input(54) );
  out(5)  += 0.125 * ( input(1)  + input(1)  + input(3)  + input(3)  + input(9)  + input(27) + input(9)  + input(27) );
  out(6)  += 0.125 * ( input(40) + input(40) + input(40) + input(40) + input(40) + input(40) + input(40) + input(40) );
  out(7)  += 0.125 * ( input(44) + input(44) + input(44) + input(44) + input(76) + input(76) + input(76) + input(76) );
  out(8)  += 0.125 * ( input(41) + input(41) + input(43) + input(43) + input(49) + input(67) + input(49) + input(67) );
  out(9)  += 0.125 * ( input(42) + input(42) + input(38) + input(38) + input(58) + input(22) + input(58) + input(22) );
  out(10) += 0.125 * ( input(39) + input(39) + input(37) + input(37) + input(31) + input(13) + input(31) + input(13) );
  out(11) += 0.125 * ( input(80) + input(80) + input(80) + input(80) + input(80) + input(80) + input(80) + input(80) );
  out(12) += 0.125 * ( input(79) + input(79) + input(77) + input(77) + input(71) + input(53) + input(71) + input(53) );
  out(13) += 0.125 * ( input(62) + input(26) + input(62) + input(26) + input(74) + input(78) + input(74) + input(74) );
  out(14) += 0.125 * ( input(35) + input(17) + input(35) + input(17) + input(73) + input(75) + input(73) + input(73) );
  out(15) += 0.125 * ( input(70) + input(52) + input(68) + input(50) + input(68) + input(52) + input(68) + input(50) );
  out(16) += 0.125 * ( input(69) + input(51) + input(65) + input(47) + input(59) + input(25) + input(59) + input(23) );
  out(17) += 0.125 * ( input(46) + input(64) + input(48) + input(66) + input(16) + input(32) + input(16) + input(34) );
  out(18) += 0.125 * ( input(20) + input(56) + input(24) + input(60) + input(24) + input(56) + input(24) + input(60) );
  out(19) += 0.125 * ( input(21) + input(57) + input(19) + input(55) + input(33) + input(11) + input(33) + input(15) );
  out(20) += 0.125 * ( input(10) + input(28) + input(12) + input(30) + input(12) + input(28) + input(12) + input(30) );

  return out;
}

SymmTensor4_3D STen4::gen_rand(const double &left, const double &right)
{
  std::random_device rd;
  std::mt19937_64 gen( rd() );
  std::uniform_real_distribution<double> dis(left, right);
  const std::array<double,21> temp {{ dis(gen), dis(gen), dis(gen), dis(gen), dis(gen), dis(gen), dis(gen),
    dis(gen), dis(gen), dis(gen), dis(gen), dis(gen), dis(gen), dis(gen),
    dis(gen), dis(gen), dis(gen), dis(gen), dis(gen), dis(gen), dis(gen) }};
  return SymmTensor4_3D(temp);
}

SymmTensor4_3D STen4::gen_Ptilde( const SymmTensor2_3D &invC )
{     
  SymmTensor4_3D out = STen4::gen_zero();
  out.add_SymmProduct( 1.0, invC, invC );
  out.add_OutProduct( -1.0/3.0, invC );
  return out;
}

SymmTensor4_3D STen4::gen_dyad( const SymmTensor4_3D &input )
{
  SymmTensor4_3D out = STen4::gen_zero();
  out(0) += input(0) * input(0);
  out(0) += input(5) * input(5);
  out(0) += input(4) * input(4);
  out(0) += input(5) * input(5);
  out(0) += input(1) * input(1);
  out(0) += input(3) * input(3);
  out(0) += input(4) * input(4);
  out(0) += input(3) * input(3);
  out(0) += input(2) * input(2);

  out(1) += input(0) * input(1);
  out(1) += input(5) * input(10);
  out(1) += input(4) * input(9);
  out(1) += input(5) * input(10);
  out(1) += input(1) * input(6);
  out(1) += input(3) * input(8);
  out(1) += input(4) * input(9);
  out(1) += input(3) * input(8);
  out(1) += input(2) * input(7);

  out(2) += input(0) * input(2);
  out(2) += input(5) * input(14);
  out(2) += input(4) * input(13);
  out(2) += input(5) * input(14);
  out(2) += input(1) * input(7);
  out(2) += input(3) * input(12);
  out(2) += input(4) * input(13);
  out(2) += input(3) * input(12);
  out(2) += input(2) * input(11);

  out(3) += input(0) * input(3);
  out(3) += input(5) * input(17);
  out(3) += input(4) * input(16);
  out(3) += input(5) * input(17);
  out(3) += input(1) * input(8);
  out(3) += input(3) * input(15);
  out(3) += input(4) * input(16);
  out(3) += input(3) * input(15);
  out(3) += input(2) * input(12);

  out(4) += input(0) * input(4);
  out(4) += input(5) * input(19);
  out(4) += input(4) * input(18);
  out(4) += input(5) * input(19);
  out(4) += input(1) * input(9);
  out(4) += input(3) * input(16);
  out(4) += input(4) * input(18);
  out(4) += input(3) * input(16);
  out(4) += input(2) * input(13);

  out(5) += input(0) * input(5);
  out(5) += input(5) * input(20);
  out(5) += input(4) * input(19);
  out(5) += input(5) * input(20);
  out(5) += input(1) * input(10);
  out(5) += input(3) * input(17);
  out(5) += input(4) * input(19);
  out(5) += input(3) * input(17);
  out(5) += input(2) * input(14);

  out(6) += input(1)  * input(1);
  out(6) += input(10) * input(10);
  out(6) += input(9)  * input(9);
  out(6) += input(10) * input(10);
  out(6) += input(6)  * input(6);
  out(6) += input(8)  * input(8);
  out(6) += input(9)  * input(9);
  out(6) += input(8)  * input(8);
  out(6) += input(7)  * input(7);

  out(7) += input(2)  * input(1);
  out(7) += input(14) * input(10);
  out(7) += input(13) * input(9);
  out(7) += input(14) * input(10);
  out(7) += input(7)  * input(6);
  out(7) += input(12) * input(8);
  out(7) += input(13) * input(9);
  out(7) += input(12) * input(8);
  out(7) += input(11) * input(7);

  out(8) += input(3)  * input(1);
  out(8) += input(17) * input(10);
  out(8) += input(16) * input(9);
  out(8) += input(17) * input(10);
  out(8) += input(8)  * input(6);
  out(8) += input(15) * input(8);
  out(8) += input(16) * input(9);
  out(8) += input(15) * input(8);
  out(8) += input(12) * input(7);

  out(9) += input(4)  * input(1);
  out(9) += input(19) * input(10);
  out(9) += input(18) * input(9);
  out(9) += input(19) * input(10);
  out(9) += input(9)  * input(6);
  out(9) += input(16) * input(8);
  out(9) += input(18) * input(9);
  out(9) += input(16) * input(8);
  out(9) += input(13) * input(7);

  out(10) += input(1)  * input(5);
  out(10) += input(10) * input(20);
  out(10) += input(9)  * input(19);
  out(10) += input(10) * input(20);
  out(10) += input(6)  * input(10);
  out(10) += input(8)  * input(17);
  out(10) += input(9)  * input(19);
  out(10) += input(8)  * input(17);
  out(10) += input(7)  * input(14);

  out(11) += input(2)  * input(2);
  out(11) += input(14) * input(14);
  out(11) += input(13) * input(13);
  out(11) += input(14) * input(14);
  out(11) += input(7)  * input(7);
  out(11) += input(12) * input(12);
  out(11) += input(13) * input(13);
  out(11) += input(12) * input(12);
  out(11) += input(11) * input(11);

  out(12) += input(2)  * input(3);
  out(12) += input(14) * input(17);
  out(12) += input(13) * input(16);
  out(12) += input(14) * input(17);
  out(12) += input(7)  * input(8);
  out(12) += input(12) * input(15);
  out(12) += input(13) * input(16);
  out(12) += input(12) * input(15);
  out(12) += input(11) * input(12);

  out(13) += input(2)  * input(4);
  out(13) += input(14) * input(19);
  out(13) += input(13) * input(18);
  out(13) += input(14) * input(19);
  out(13) += input(7)  * input(9);
  out(13) += input(12) * input(16);
  out(13) += input(13) * input(18);
  out(13) += input(12) * input(16);
  out(13) += input(11) * input(13);

  out(14) += input(2)  * input(5);
  out(14) += input(14) * input(20);
  out(14) += input(13) * input(19);
  out(14) += input(14) * input(20);
  out(14) += input(7)  * input(10);
  out(14) += input(12) * input(17);
  out(14) += input(13) * input(19);
  out(14) += input(12) * input(17);
  out(14) += input(11) * input(14);

  out(15) += input(3)  * input(3);
  out(15) += input(17) * input(17);
  out(15) += input(16) * input(16);
  out(15) += input(17) * input(17);
  out(15) += input(8)  * input(8);
  out(15) += input(15) * input(15);
  out(15) += input(16) * input(16);
  out(15) += input(15) * input(15);
  out(15) += input(12) * input(12);

  out(16) += input(3)  * input(4);
  out(16) += input(17) * input(19);
  out(16) += input(16) * input(18);
  out(16) += input(17) * input(19);
  out(16) += input(8)  * input(9);
  out(16) += input(15) * input(16);
  out(16) += input(16) * input(18);
  out(16) += input(15) * input(16);
  out(16) += input(12) * input(13);

  out(17) += input(3)  * input(5);
  out(17) += input(17) * input(20);
  out(17) += input(16) * input(19);
  out(17) += input(17) * input(20);
  out(17) += input(8)  * input(10);
  out(17) += input(15) * input(17);
  out(17) += input(16) * input(19);
  out(17) += input(15) * input(17);
  out(17) += input(12) * input(14);

  out(18) += input(4)  * input(4);
  out(18) += input(19) * input(19);
  out(18) += input(18) * input(18);
  out(18) += input(19) * input(19);
  out(18) += input(9)  * input(9);
  out(18) += input(16) * input(16);
  out(18) += input(18) * input(18);
  out(18) += input(16) * input(16);
  out(18) += input(13) * input(13);

  out(19) += input(4)  * input(5);
  out(19) += input(19) * input(20);
  out(19) += input(18) * input(19);
  out(19) += input(19) * input(20);
  out(19) += input(9)  * input(10);
  out(19) += input(16) * input(17);
  out(19) += input(18) * input(19);
  out(19) += input(16) * input(17);
  out(19) += input(13) * input(14);

  out(20) += input(5)  * input(5);
  out(20) += input(20) * input(20);
  out(20) += input(19) * input(19);
  out(20) += input(20) * input(20);
  out(20) += input(10) * input(10);
  out(20) += input(17) * input(17);
  out(20) += input(19) * input(19);
  out(20) += input(17) * input(17);
  out(20) += input(14) * input(14);

  return out;
}

// EOF
