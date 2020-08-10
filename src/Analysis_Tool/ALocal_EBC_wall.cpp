#include "ALocal_EBC_wall.hpp"

ALocal_EBC_wall::ALocal_EBC_wall( const std::string &fileBaseName,
    const int &cpu_rank, const std::string &gname )
: ALocal_EBC( fileBaseName, cpu_rank, gname )
{
  // Read in the fluid density, thickness, and young's modulus
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );
  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF5_Reader * h5r = new HDF5_Reader( file_id );
  
  fluid_density = h5r -> read_doubleScalar(gname.c_str(), "fluid_density");

  thickness.clear();
  youngsmod.clear();

  // wall has only one surface per the assumption in wall ebc  
  const int ebc_id = 0;

  if( num_local_cell[ebc_id] > 0 )
  {
    std::string subgroup_name(gname);
    subgroup_name.append("/ebcid_");
    subgroup_name.append( SYS_T::to_string(ebc_id) );

    h5r -> read_doubleVector( subgroup_name.c_str(), "thickness", thickness );
    h5r -> read_doubleVector( subgroup_name.c_str(), "youngsmod", youngsmod );
  }

  delete h5r; H5Fclose( file_id );
}


ALocal_EBC_wall::~ALocal_EBC_wall()
{
  VEC_T::clean(thickness);
  VEC_T::clean(youngsmod);
}

void ALocal_EBC_wall::print_info() const
{
  SYS_T::commPrint("---- ALocal_EBC_wall: \n");
  SYS_T::commPrint("     num_ebc = %d \n", num_ebc);
}

// EOF
