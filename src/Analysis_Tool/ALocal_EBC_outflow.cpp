#include "ALocal_EBC_outflow.hpp"

ALocal_EBC_outflow::ALocal_EBC_outflow( const std::string &fileBaseName,
    const int &cpu_rank, const std::string &gname )
: ALocal_EBC( fileBaseName, cpu_rank, gname )
{
  // Read in the integral of basis function, LID, and outward normal vector
  // for the faces with id ebc
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );
  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF5_Reader * h5r = new HDF5_Reader( file_id );
  
  intNA.resize(  num_ebc ); LID.resize( num_ebc );
  outvec.resize( num_ebc ); num_face_nodes.resize( num_ebc );

  for(int ii=0; ii<num_ebc; ++ii)
  {
    if( num_local_cell[ii] > 0 )
    {
      std::string subgroup_name(gname);
      subgroup_name.append("/ebcid_");
      subgroup_name.append( std::to_string(ii) );

      intNA[ii]  = h5r -> read_doubleVector( subgroup_name.c_str(), "intNA" );
      LID[ii]    = h5r -> read_intVector(    subgroup_name.c_str(), "LID_all_face_nodes" );
      outvec[ii] = Vector_3( h5r -> read_Vector_3( subgroup_name.c_str(), "out_normal" ) );
    
      num_face_nodes[ii] = static_cast<int>(intNA[ii].size());
    }
    else
    {
      intNA[ii].clear();
      LID[ii].clear();
      outvec[ii] = Vector_3( 0.0, 0.0, 0.0 );
      num_face_nodes[ii] = 0;
    }
  }

  delete h5r; H5Fclose( file_id );
}

void ALocal_EBC_outflow::print_info() const
{
  SYS_T::commPrint("---- ALocal_EBC_outflow: \n");
  SYS_T::commPrint("     num_ebc = %d \n", num_ebc);
}

// EOF
