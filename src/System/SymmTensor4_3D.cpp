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
// EOF
