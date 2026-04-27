#include "ALocal_NBC_Solid.hpp"
#include "HDF5_Reader.hpp"

ALocal_NBC_Solid::ALocal_NBC_Solid( const std::string &fileBaseName,
    const int &cpu_rank, const std::string &gname )
: ALocal_NBC( fileBaseName, cpu_rank, gname )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  auto h5r = SYS_T::make_unique<HDF5_Reader>(fName);

  read_disp_flag( h5r.get(), gname );
}

ALocal_NBC_Solid::ALocal_NBC_Solid( const HDF5_Reader * const &h5r,
    const std::string &gname )
: ALocal_NBC( h5r, gname )
{
  read_disp_flag( h5r, gname );
}

void ALocal_NBC_Solid::read_disp_flag( const HDF5_Reader * const &h5r,
    const std::string &gname )
{
  const int total_ld = VEC_T::sum( Num_LD );

  if( total_ld <= 0 )
  {
    LDN_is_disp_driven.clear();
    return;
  }

  const std::string data_path = gname + "/LDN_is_disp_driven";
  if( h5r->check_data( data_path.c_str() ) )
  {
    const auto raw_flag = h5r->read_intVector( gname.c_str(), "LDN_is_disp_driven" );
    SYS_T::print_fatal_if( int(raw_flag.size()) != total_ld,
        "Error: ALocal_NBC_Solid, LDN_is_disp_driven length does not match Num_LD.\n" );
    LDN_is_disp_driven.assign( raw_flag.begin(), raw_flag.end() );
  }
  else
  {
    LDN_is_disp_driven.assign( total_ld, false );
  }
}

// EOF
