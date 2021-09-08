#include "Prestress_solid.hpp"

Prestress_solid::Prestress_solid( 
    const ALocal_Elem * const &locelem, 
    const int &in_nqp_tet, const std::string &in_ps_fName )
: nlocalele(locelem->get_nlocalele()), nqp(in_nqp_tet), 
  ps_fileBaseName( in_ps_fName )
{
  sval.resize( nlocalele );

  for(int ee=0; ee<nlocalele; ++ee)
  {
    if( locelem->get_elem_tag(ee) == 0 ) sval[ee].clear();
    else if( locelem->get_elem_tag(ee) == 1 ) sval[ee].resize(nqp * 6);
    else
      SYS_T::print_fatal("Error: element tag should be 0 (fluid) or 1 (solid).\n");
  }
}

Prestress_solid::~Prestress_solid()
{
  for(int ee=0; ee<nlocalele; ++ee) VEC_T::clean(sval[ee]);

  VEC_T::clean(sval);
}

std::vector<double> Prestress_solid::get_prestress( const int &ee ) const
{
  return sval[ee];
}

void Prestress_solid::set_prestress( const int &ee, const int &ii,
    const double * const &in_esval )
{
  for(int ii=0; ii<6*nqp; ++ii) sval[ee][ii] += in_esval[ii];
}


void Prestress_solid::print_info() const
{
  for(int ee=0; ee<nlocalele; ++ee)
  {
    std::cout<<"ee: "<<ee<<'\t';
    VEC_T::print(sval[ee]);
    std::cout<<std::endl;
  }
}

void Prestress_solid::write_prestress_hdf5() const
{}


void Prestress_solid::read_prestress_hdf5() const
{}

// EOF
