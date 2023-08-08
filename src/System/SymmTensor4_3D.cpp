#include "SymmTensor4_3D.hpp"

SymmTensor4_3D::SymmTensor4_3D()
{
  gen_zero();

  ten[0] = 1.0; ten[6] = 1.0; ten[11] = 1.0;
  ten[15] = 0.5; ten[18] = 0.5; ten[20] = 0.5;
}

SymmTensor4_3D::~SymmTensor4_3D()
{}

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
  for ( int ii=0; ii<3; ii++ )
  {
    for( int jj=0; jj<3; jj++ )
    {
      for( int kk=0; kk<3; kk++ )
      { 
        for ( int ll=0; ll<3; ll++ )
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
  double add_ten[21];
  for(int ii=0; ii<3; ++ii)
  {
    for(int jj=0; jj<3; ++jj)
    {
      for(int kk=0; kk<3; ++kk)
      {
        for(int ll=0; ll<3; ++ll)
          add_ten[ Voigt_notation(ii,jj,kk,ll) ] = val * mmat( Voigt_notation(ii,jj) ) * mmat( Voigt_notation(kk,ll) );
      }
    }
  }
  for(int counter=0; counter<21; counter++)
  {
    ten[counter] += add_ten[counter];
  }
}

void SymmTensor4_3D::add_SymmOutProduct( const double &val, const Vector_3 &vec1, 
  const Vector_3 &vec2, const Vector_3 &vec3, const Vector_3 &vec4 )
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
          add_ten[ Voigt_notation(ii,jj,kk,ll) ] = val * ( vec1( ii ) * vec2( jj ) * vec3( kk ) * vec4( ll )
            + vec1( ii ) * vec2( jj ) * vec3( ll ) * vec4( kk )
            + vec1( jj ) * vec2( ii ) * vec3( kk ) * vec4( ll )
            + vec1( jj ) * vec2( ii ) * vec3( ll ) * vec4( kk ) );
        }
      }
    }
  }
  for(int counter=0; counter<21; counter++)
  {
    ten[counter] += add_ten[counter];
  }
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
  for(int counter=0; counter<21; counter++)
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
  for(int counter=0; counter<21; counter++)
  {
    ten[counter] += add_ten[counter];
  }
}

int SymmTensor4_3D::Voigt_notation( const int &ii, const int &jj, const int &kk, const int &ll ) const
{
  int index_I = 3;
  int index_J = 3;
  switch (ii)
  {
  case 0:
    switch (jj)
    {
    case 0:
     index_I = 0;
     break;

   case 1:
     index_I = 5;
     break;

   case 2:
     index_I = 4;
     break;
   }
   break;

 case 1:
  switch (jj)
  {
  case 0:
   index_I = 5;
   break;

 case 1:
   index_I = 1;
   break;

 case 2:
   index_I = 3;
   break;
 }
 break;

case 2:
  switch (jj)
  {
  case 0:
   index_I = 4;
   break;

 case 1:
   index_I = 3;
   break;

 case 2:
   index_I = 2;
   break;
 }
 break;}

 switch (kk)
 {
 case 0:
  switch (ll)
  {
  case 0:
   index_J = 0;
   break;

 case 1:
   index_J = 5;
   break;

 case 2:
   index_J = 4;
   break;
 }
 break;

case 1:
  switch (ll)
  {
  case 0:
   index_J = 5;
   break;

 case 1:
   index_J = 1;
   break;

 case 2:
   index_J = 3;
   break;
 }
 break;

case 2:
  switch (ll)
  {
  case 0:
   index_J = 4;
   break;

 case 1:
   index_J = 3;
   break;

 case 2:
   index_J = 2;
   break;
 }
 break;}

 int sum = 0;
 if (index_I <= index_J){
  for (int counter=index_I; counter>0; counter--)
  {
    sum += 6 - counter;
  }
  return sum + index_J;}
  else
  {
    for (int counter=index_J; counter>0; counter--)
    {
      sum += 6 - counter;
    }

    return sum + index_I;
  }
}

int SymmTensor4_3D::Voigt_notation( const int &ii, const int &jj ) const
{
  int index = 3;
  switch (ii)
  {
  case 0:
    switch (jj)
    {
    case 0:
     index = 0;
     break;

   case 1:
     index = 5;
     break;

   case 2:
     index = 4;
     break;
   }
   break;

 case 1:
  switch (jj)
  {
  case 0:
   index = 5;
   break;

 case 1:
   index = 1;
   break;

 case 2:
   index = 3;
   break;
 }
 break;

case 2:
  switch (jj)
  {
  case 0:
   index = 4;
   break;

 case 1:
   index = 3;
   break;

 case 2:
   index = 2;
   break;
 }
 break;}

 return index;
}
// EOF
