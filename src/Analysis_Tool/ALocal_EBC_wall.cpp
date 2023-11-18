#include "ALocal_EBC_wall.hpp"

ALocal_EBC_wall::ALocal_EBC_wall( const std::string &fileBaseName,
    const int &in_cpu_rank, const int &in_face_nqp, 
    const std::string &gname, const std::string &input_ps_filebasename )
: ALocal_EBC( fileBaseName, in_cpu_rank, gname ), cpu_rank (in_cpu_rank),
  face_nqp( in_face_nqp ), ps_fileBaseName( input_ps_filebasename )
{
  SYS_T::print_fatal_if(gname != "ebc_wall", 
      "Error: ALocal_EBC_wall data should be read from group ebc_wall.\n" );

  // Read in the thickness, and young's modulus
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );
  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF5_Reader * h5r = new HDF5_Reader( file_id );
  
  num_local_node_on_sur = h5r -> read_intScalar( gname.c_str(), "num_local_node_on_sur" );

  thickness.clear();
  youngsmod.clear();
  springconst.clear();
  dampingconst.clear();
  qua_prestress.clear();

  const int ebc_id = 0;
  // wall has only one surface per the assumption in wall ebc  
  if( num_local_cell[ebc_id] > 0 )
  {
    std::string subgroup_name(gname);
    subgroup_name.append("/ebcid_0");

    thickness = h5r -> read_doubleVector( subgroup_name.c_str(), "thickness" );
    youngsmod = h5r -> read_doubleVector( subgroup_name.c_str(), "youngsmod" );

    springconst  = h5r -> read_doubleVector( subgroup_name.c_str(), "springconst"  );
    dampingconst = h5r -> read_doubleVector( subgroup_name.c_str(), "dampingconst" );

    // If the prestress data exist on file, read the prestress data
    // otherwise, just allocate a container for prestress data at each quadrature 
    // point with zero values
    const std::string ps_fName = SYS_T::gen_partfile_name( ps_fileBaseName, cpu_rank );
    if( SYS_T::file_exist(ps_fName) )
    {
      hid_t ps_file_id = H5Fopen(ps_fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

      HDF5_Reader * ps_h5r = new HDF5_Reader( ps_file_id );

      qua_prestress = ps_h5r -> read_doubleVector( "/", "prestress" );

      delete ps_h5r; H5Fclose(ps_file_id);

      SYS_T::print_fatal_if( static_cast<int>( qua_prestress.size() ) != 6 * face_nqp * num_local_cell[ebc_id],
        "ALocal_EBC_wall: size of qua_prestress is inconsistent with face_nqp. \n");
    
      SYS_T::commPrint("===> ALocal_EBC_wall : qua_prestress loaded from %s.\n", ps_fName.c_str());
    }
    else
    {
      qua_prestress.resize( 6 * face_nqp * num_local_cell[ebc_id] );
      
      for( auto &val : qua_prestress ) val = 0.0;
      
      SYS_T::commPrint("===> ALocal_EBC_wall : qua_prestress initialized to be zero.\n");
    }
  }

  local_node_on_sur_pos.clear();
  if( num_local_node_on_sur > 0 )
  {
    std::string subgroup_name(gname);
    subgroup_name.append("/ebcid_0");

    local_node_on_sur_pos = h5r -> read_intVector( subgroup_name.c_str(), "local_node_on_sur_pos" );
  }

  delete h5r; H5Fclose( file_id );
}

void ALocal_EBC_wall::get_thickness( const int &eindex,
    double * const &e_thickness ) const
{
  // Only one surface per the assumption in wall ebc
  const int nLocBas = cell_nLocBas[0];
  for( int ii = 0; ii < nLocBas; ++ii )
  {
    const int pos = local_cell_ien[0][nLocBas*eindex + ii]; 
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
    const int pos = local_cell_ien[0][nLocBas*eindex + ii]; 
    e_youngsmod[ii] = youngsmod[pos];
  }
}

void ALocal_EBC_wall::get_springconst( const int &eindex,
    double * const &e_springconst ) const
{
  // Only one surface per the assumption in wall ebc
  const int nLocBas = cell_nLocBas[0];
  for( int ii = 0; ii < nLocBas; ++ii )
  {
    const int pos = local_cell_ien[0][nLocBas*eindex + ii]; 
    e_springconst[ii] = springconst[pos];
  }
}

void ALocal_EBC_wall::get_dampingconst( const int &eindex,
    double * const &e_dampingconst ) const
{
  // Only one surface per the assumption in wall ebc
  const int nLocBas = cell_nLocBas[0];
  for( int ii = 0; ii < nLocBas; ++ii )
  {
    const int pos = local_cell_ien[0][nLocBas*eindex + ii]; 
    e_dampingconst[ii] = dampingconst[pos];
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
  const int pos = 6 * eindex * face_nqp; 

  for(int ii = 0; ii < 6 * face_nqp; ++ii)
    qua_prestress[pos + ii] = e_quaprestress[ii]; 
}

void ALocal_EBC_wall::write_prestress_hdf5() const
{
  const std::string fName = SYS_T::gen_partfile_name( ps_fileBaseName, cpu_rank );

  // open the folder at fName/ebc_wall again to append additional data 
  // num_ebc = 1 for wall elem bc and id = 0 
  const int ebc_id = 0;
  if( num_local_cell[ebc_id] > 0 )
  {
    // re-open the file
    hid_t file_id = H5Fcreate(fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    HDF5_Writer * h5w = new HDF5_Writer( file_id );
      
    h5w->write_doubleVector( file_id, "prestress", qua_prestress );

    delete h5w; H5Fclose( file_id );
  }
}

void ALocal_EBC_wall::print_info() const
{
  SYS_T::commPrint("---- ALocal_EBC_wall: \n");
  SYS_T::commPrint("     num_ebc = %d \n", num_ebc);
}

// EOF
