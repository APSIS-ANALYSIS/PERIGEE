#include "ALocal_EBC_wall.hpp"

ALocal_EBC_wall::ALocal_EBC_wall( const std::string &fileBaseName,
    const int &in_cpu_rank, const IQuadPts * const &quad, 
    const std::string &gname, const bool &prestress_flag)
: ALocal_EBC( fileBaseName, in_cpu_rank, gname ), cpu_rank (in_cpu_rank),
  face_nqp( quad -> get_num_quadPts() ), solve_prestress( prestress_flag )
{
  SYS_T::print_fatal_if(gname != "ebc_wall", 
      "Error: ALocal_EBC_wall data should be read from group ebc_wall.\n" );

  // Read in the fluid density, thickness, and young's modulus
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );
  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF5_Reader * h5r = new HDF5_Reader( file_id );
  
  fluid_density = h5r -> read_doubleScalar(gname.c_str(), "fluid_density");

  thickness.clear();
  youngsmod.clear();
  qua_prestress.clear();

  // wall has only one surface per the assumption in wall ebc  
  if( num_local_cell[0] > 0 )
  {
    std::string subgroup_name(gname);
    subgroup_name.append("/ebcid_0");

    h5r -> read_doubleVector( subgroup_name.c_str(), "thickness", thickness );
    h5r -> read_doubleVector( subgroup_name.c_str(), "youngsmod", youngsmod );

    // If the prestress is solved, read the prestress data
    // otherwise, just allocate a container for prestress data at each quadrature point
    if( false )
    {
      h5r -> read_doubleVector( subgroup_name.c_str(), "prestress", qua_prestress );

      SYS_T::print_fatal_if( static_cast<int>( qua_prestress.size() ) != 6 * face_nqp * num_local_cell[0],
        "ALocal_EBC_wall: size of qua_prestress is inconsistent with face_nqp. \n");
    }
    else
    {
      qua_prestress.resize( 6 * face_nqp * num_local_cell[0] );
      for(int ii=0; ii< 6 * face_nqp * num_local_cell[0]; ++ii) qua_prestress[ii] = 0.0;
    }
  }

  delete h5r; H5Fclose( file_id );
}

ALocal_EBC_wall::~ALocal_EBC_wall()
{
  VEC_T::clean(thickness);
  VEC_T::clean(youngsmod);
  VEC_T::clean(qua_prestress);
}

void ALocal_EBC_wall::get_thickness( const int &eindex,
    double * const &e_thickness ) const
{
  // Only one surface per the assumption in wall ebc
  const int nLocBas = cell_nLocBas[0];
  for( int ii = 0; ii < nLocBas; ++ii )
  {
    const int pos = local_tri_ien[0][nLocBas*eindex + ii]; 
    e_thickness[ii] = thickness[pos];
  }
}

void ALocal_EBC_wall::get_youngsmod( const int &eindex,
    double * const &e_youngsmod ) const
{
  // Only one surface per the assumption in wall ebc
  const int nLocBas = cell_nLocBas[0];
  for( int ii = 0; ii < nLocBas; ++ii )
  {
    const int pos = local_tri_ien[0][nLocBas*eindex + ii]; 
    e_youngsmod[ii] = youngsmod[pos];
  }
}

void ALocal_EBC_wall::get_prestress( const int &eindex,
    double * const &e_quaprestress ) const
{
  const int pos = 6 * eindex * face_nqp; 

  for(int ii = 0; ii < 6 * face_nqp; ++ii)
    e_quaprestress[ii] = qua_prestress[pos + ii]; 
}

void ALocal_EBC_wall::set_prestress( const int &eindex,
    const double * const &e_quaprestress )
{
  if( solve_prestress )
  {
    const int pos = 6 * eindex * face_nqp; 

    for(int ii = 0; ii < 6 * face_nqp; ++ii)
      qua_prestress[pos + ii] = e_quaprestress[ii]; 
  }
  else
    SYS_T::commPrint("Warning in ALocal_EBC_wall: set_prestress should only be called when solving for prestress. \n"); 
}

void ALocal_EBC_wall::write_prestress_hdf5( const char * FileName ) const
{
  const std::string input_fName(FileName);
  const std::string fName = SYS_T::gen_partfile_name( input_fName, cpu_rank );

  // re-open the file
  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  // open the folder at fName/ebc_wall again to append additional data 
  const int ebc_id = 0;              // num_ebc = 1 for wall elem bc and id = 0
  if( num_local_cell[ebc_id] > 0 )
  {
    hid_t group_id = H5Gopen( file_id, "ebc_wall/ebcid_0", H5P_DEFAULT );

    HDF5_Writer * h5w = new HDF5_Writer( file_id );

    h5w->write_doubleVector( group_id, "prestress", qua_prestress );

    H5Gclose( group_id );
    delete h5w; H5Gclose( group_id ); H5Fclose( file_id );
  }
}

void ALocal_EBC_wall::print_info() const
{
  SYS_T::commPrint("---- ALocal_EBC_wall: \n");
  SYS_T::commPrint("     num_ebc = %d \n", num_ebc);
}

// EOF
