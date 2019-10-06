#include "Prestress_solid.hpp"

Prestress_solid::Prestress_solid( 
    const ALocal_Elem * const &locelem, 
    const int &in_nqp_tet )
: nlocelem(locelem->get_nlocalele()), nqp(in_nqp_tet)
{
  sval.resize( nlocelem );

  for(int ee=0; ee<nlocelem; ++ee)
  {
    if( locelem->get_elem_tag(ee) == 0 )
      sval[ee].clear();
    else if( locelem->get_elem_tag(ee) == 1 )
      sval[ee].resize(nqp * 7);
    else
      SYS_T::print_fatal("Error: element tag should be 0 (fluid) or 1 (solid).\n");
  }
}


Prestress_solid::~Prestress_solid()
{
  for(int ee=0; ee<nlocelem; ++ee) VEC_T::clean(sval[ee]);

  VEC_T::clean(sval);
}


void Prestress_solid::get_sval(const int &ee, 
    std::vector<double> &esval ) const
{
  esval = sval[ee];
}


void Prestress_solid::update_sval(const int &ee,
    const std::vector<double> &in_esval )
{
  for(int ii=0; ii<7*nqp; ++ii) sval[ee][ii] += in_esval[ii];
}


void Prestress_solid::print_info() const
{
  for(int ee=0; ee<nlocelem; ++ee)
  {
    std::cout<<"ee: "<<ee<<'\t';
    VEC_T::print(sval[ee]);
    std::cout<<std::endl;
  }
}


// EOF
