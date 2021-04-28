#include "ALocal_Ring_NodalBC.hpp"

ALocal_Ring_NodalBC::ALocal_Ring_NodalBC( 
    const std::string &fileBaseName, const int &cpu_rank )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/ring");
  
  Num_LD = h5r -> read_intScalar( gname.c_str(), "Num_LD" );

  num_caps = h5r -> read_intScalar( gname.c_str(), "num_caps" );

  h5r -> read_intVector( gname.c_str(), "cap_dominant_comp", dominant_n_comp );

  std::vector<double> outnormal_vec;
  h5r -> read_doubleVector( gname.c_str(), "cap_out_normal", outnormal_vec );

  outnormal.resize(num_caps);

  for(int ii=0; ii<num_caps; ++ii)
  {
    outnormal[ii](0) = outnormal_vec[3*ii+0];
    outnormal[ii](1) = outnormal_vec[3*ii+1];
    outnormal[ii](2) = outnormal_vec[3*ii+2];
  }

  // each cap's centroid xyz coordinates
  std::vector<double> centroid; // length is 3 x num_caps

  h5r -> read_doubleVector( gname.c_str(), "cap_centroid", centroid );

  // If this sub-domain contains local ring nodes,
  // load the LDN array, the corresponding cap ids, and the ring nodes'
  // coordinates
  std::vector<double> local_pt_xyz; // length is 3 x Num_LD
  if( Num_LD > 0 )
  {
    h5r->read_intVector( gname.c_str(), "LDN", LDN );
    h5r->read_intVector( gname.c_str(), "local_cap_id", local_cap_id );
    h5r->read_doubleVector( gname.c_str(), "local_pt_xyz", local_pt_xyz );
  }

  delete h5r; H5Fclose( file_id );
  
  // Compute the tangential vector for each local ring node
  tangential.resize(Num_LD);
  dominant_t_comp.resize(Num_LD);

  for( int node=0; node<Num_LD; ++node)
  {
    // Generate radial vector using nodal & centroidal coordinates
    Vector_3 radial_vec = Vector_3(
        local_pt_xyz[3*node]     - centroid[ 3*local_cap_id[node]    ],
        local_pt_xyz[3*node + 1] - centroid[ 3*local_cap_id[node] + 1],
        local_pt_xyz[3*node + 2] - centroid[ 3*local_cap_id[node] + 2] );

    tangential[node] = cross_product( outnormal[ local_cap_id[node] ], radial_vec ); 
    
    tangential[node].normalize();

    dominant_t_comp[node] = tangential[node].get_dominant_comp();

    SYS_T::print_fatal_if(dominant_t_comp[node] == dominant_n_comp[ local_cap_id[node] ], "Error: ALocal_Ring_NodalBC the tangential and normal vector have the same dominant component.\n");
  }

  VEC_T::clean( local_pt_xyz );
  VEC_T::clean( centroid     );
}

ALocal_Ring_NodalBC::~ALocal_Ring_NodalBC()
{
  VEC_T::clean( LDN             );
  VEC_T::clean( local_cap_id    );
  VEC_T::clean( dominant_n_comp );
  VEC_T::clean( dominant_t_comp );
  VEC_T::clean( outnormal       );
  VEC_T::clean( tangential      );
}

// EOF
