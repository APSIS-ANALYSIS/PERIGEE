#include "ALocal_EBC_wall.hpp"

ALocal_EBC_wall::ALocal_EBC_wall( const std::string &fileBaseName,
    const int &cpu_rank, const IQuadPts * const &quad, 
    const std::string &gname, const bool &prestress_flag)
: ALocal_EBC( fileBaseName, cpu_rank, gname ), face_nqp( quad -> get_num_quadPts() ),
  solve_prestress( prestress_flag )
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

   if( !solve_prestress )
   {
     h5r -> read_doubleVector( subgroup_name.c_str(), "prestress", qua_prestress );

     SYS_T::print_fatal_if( static_cast<int>( qua_prestress.size() ) != 6 * face_nqp * num_local_cell[0],
       "ALocal_EBC_wall: size of qua_prestress is inconsistent with face_nqp. \n");
   }
   else
     qua_prestress.resize( 6 * face_nqp * num_local_cell[0] );
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
    double * const &e_quaprestress )
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


void ALocal_EBC_wall::print_info() const
{
  SYS_T::commPrint("---- ALocal_EBC_wall: \n");
  SYS_T::commPrint("     num_ebc = %d \n", num_ebc);
}

// EOF
