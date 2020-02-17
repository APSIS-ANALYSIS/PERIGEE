#include "ALocal_Solid_Elem_Data.hpp"

ALocal_Solid_Elem_Data::ALocal_Solid_Elem_Data( 
    const ALocal_Elem_wTag * const &in_aelem,
    const int &in_nqp )
: nqp( in_nqp ), 
  nlocalelem(in_aelem->get_nlocalele())
{
  nelem_solid = 0;

  prestress = new double * [nlocalelem];

  for(int ee=0; ee<nlocalelem; ++ee)
  {
    if( in_aelem->get_elem_tag(ee) == 0 )
    {
      prestress[ee] = nullptr;
    }
    else if( in_aelem->get_elem_tag(ee) == 1 )
    {
      nelem_solid += 1; // count the solid element
      prestress[ee] = new double [7*nqp];
      for(int ii=0; ii<7*nqp; ++ii) 
        prestress[ee][ii] = 0.0;
    }
  }
}

ALocal_Solid_Elem_Data::~ALocal_Solid_Elem_Data()
{
  for(int ee=0; ee<nlocalelem; ++ee)
  {
    if(prestress[ee] != nullptr)
      delete [] prestress[ee];
  }
  delete [] prestress;
}

void ALocal_Solid_Elem_Data::print_info()
{
  std::cout<<"ALocal_Solid_Elem_Data: \n";
  std::cout<<"solid element number is "<<nelem_solid<<std::endl;
}


void ALocal_Solid_Elem_Data::get_prestress( const int &ee,
    const int &ii, const double &S00, const double &S01,
    const double &S02, const double &S11,
    const double &S12, const double &S22,
    const double &pressure ) const
{
  if( prestress[ee] == nullptr )
  {
    S00 = 0.0;
    S01 = 0.0;
    S02 = 0.0;
    S11 = 0.0;
    S12 = 0.0;
    S22 = 0.0;
    pressure = 0.0;
  }
  else
  {
    S00 = prestress[ee][7*ii+0];
    S01 = prestress[ee][7*ii+1];
    S02 = prestress[ee][7*ii+2];
    S11 = prestress[ee][7*ii+3];
    S12 = prestress[ee][7*ii+4];
    S22 = prestress[ee][7*ii+5];
    pressure = prestress[ee][7*ii+6];
  }
}

// EOF
