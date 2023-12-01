#include "ALocal_InflowBC.hpp"

ALocal_InflowBC::ALocal_InflowBC( 
    const std::string &fileBaseName, const int &cpu_rank )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/inflow");

  num_nbc = h5r -> read_intScalar( gname.c_str(), "num_nbc" );
    
  // Allocate the size of the member data
  Num_LD.resize(num_nbc); 
  LDN.resize(num_nbc);
  outnormal.resize(num_nbc); 
  act_area.resize(num_nbc); 
  ful_area.resize(num_nbc); 
  num_out_bc_pts.resize(num_nbc);
  outline_pts.resize(num_nbc);
  outline_pts_part_tag.resize(num_nbc);
  outline_pts_loc.resize(num_nbc);
  centroid.resize(num_nbc); 
  num_local_node.resize(num_nbc); 
  num_local_cell.resize(num_nbc); 
  cell_nLocBas.resize(num_nbc);
  local_pt_xyz.resize(num_nbc); 
  local_cell_ien.resize(num_nbc); 
  local_node_pos.resize(num_nbc);

  std::string groupbase(gname);
  groupbase.append("/nbcid_");

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    std::string subgroup_name(groupbase);
    subgroup_name.append( std::to_string(nbc_id) );

    Num_LD[nbc_id] = h5r -> read_intScalar( subgroup_name.c_str(), "Num_LD" );

    // If this sub-domain of this CPU contains local inflow bc points, load the LDN array.
    // Also load the centroid and outline_pts as they three will be used for
    // generating the flow profile on the inlet at the nodes.
    if( Num_LD[nbc_id] > 0 )
    {
      LDN[nbc_id]            = h5r -> read_intVector(    subgroup_name.c_str(), "LDN" );
    }
    else
    {
      LDN[nbc_id].clear();
    }
    
    SYS_T::print_fatal_if( Num_LD[nbc_id] != static_cast<int>( LDN[nbc_id].size() ), "Error: the LDN vector size does not match with the value of Num_LD.\n" );

    // Basic geometrical quantities of the nbc_id-th inlet surface
    act_area[nbc_id]       = h5r -> read_doubleScalar( subgroup_name.c_str(), "Inflow_active_area" );
    ful_area[nbc_id]       = h5r -> read_doubleScalar( subgroup_name.c_str(), "Inflow_full_area" );
    outline_pts[nbc_id]    = h5r -> read_doubleVector( subgroup_name.c_str(), "outline_pts" );
    num_out_bc_pts[nbc_id] = h5r -> read_intScalar(    subgroup_name.c_str(), "num_out_bc_pts" );
    outline_pts_part_tag[nbc_id] = h5r -> read_intVector( subgroup_name.c_str(), "outline_pts_part_tag" );
    outline_pts_loc[nbc_id]= h5r -> read_intVector(    subgroup_name.c_str(), "outline_pts_loc" );
    centroid[nbc_id]       = h5r -> read_Vector_3(     subgroup_name.c_str(), "centroid" );
    outnormal[nbc_id]      = h5r -> read_Vector_3(     subgroup_name.c_str(), "Outward_normal_vector" );

    // Cell-related data
    cell_nLocBas[nbc_id]   = h5r->read_intScalar( subgroup_name.c_str(), "cell_nLocBas" );
    num_local_cell[nbc_id] = h5r->read_intScalar( subgroup_name.c_str(), "num_local_cell" );
    num_local_node[nbc_id] = h5r->read_intScalar( subgroup_name.c_str(), "num_local_node" );

    // If this partitioned sub-domain contains inlet surface element,
    // load its geometrical info 
    if(num_local_cell[nbc_id] > 0)
    {
      const std::vector<double> temp_xyz = h5r->read_doubleVector( subgroup_name.c_str(), "local_pt_xyz" );

      ASSERT( VEC_T::get_size(temp_xyz) == num_local_node[nbc_id]*3, "Error: ALocal_InflowBC local_pt_xyz format is wrong.\n");

      local_pt_xyz[nbc_id] = std::vector<Vector_3> (num_local_node[nbc_id], Vector_3{ 0, 0, 0 });
      
      for(int ii {0}; ii < num_local_node[nbc_id]; ++ii)
        local_pt_xyz[nbc_id][ii] = Vector_3{ temp_xyz[3 * ii], temp_xyz[3 * ii + 1], temp_xyz[3 * ii + 2] };
      
      local_cell_ien[nbc_id] = h5r->read_intVector( subgroup_name.c_str(), "local_cell_ien" );
      local_node_pos[nbc_id] = h5r->read_intVector( subgroup_name.c_str(), "local_node_pos" );
    }
    else
    {
      local_pt_xyz[nbc_id].clear();
      local_cell_ien[nbc_id].clear();
      local_node_pos[nbc_id].clear();
    }
  } // end nbc_id-loop

  delete h5r; H5Fclose( file_id );
}

double ALocal_InflowBC::get_radius( const int &nbc_id,
    const Vector_3 &pt ) const
{
  // num_out_bc_pts is set to be zero for parition that does not contain
  // inflow boundary points (i.e. Num_LD = 0 ).
  SYS_T::print_fatal_if( num_out_bc_pts[nbc_id] == 0, "Error: ALocal_InflowBC::get_radius, this function can only be called in sub-domains which contains the inflow boundary node.\n");

  const double rc = Vec3::dist( pt, centroid[nbc_id] );

  // Now loop over the boundary points to find rb.
  double rb = Vec3::dist( pt, Vector_3( outline_pts[nbc_id][0], outline_pts[nbc_id][1], outline_pts[nbc_id][2]) );

  for(int ii=1; ii<num_out_bc_pts[nbc_id]; ++ii)
  {
    double newdist = Vec3::dist( pt, Vector_3( outline_pts[nbc_id][3*ii], outline_pts[nbc_id][3*ii+1], outline_pts[nbc_id][3*ii+2]) );

    if(newdist < rb) rb = newdist;
  }

  return rc / (rb + rc);
}

void ALocal_InflowBC::get_ctrlPts_xyz( const int &nbc_id,
    const int &eindex, double * const &ctrl_x, double * const &ctrl_y, 
    double * const &ctrl_z ) const
{
  for(int jj=0; jj<cell_nLocBas[nbc_id]; ++jj)
  {
    const int pos = local_cell_ien[nbc_id][ cell_nLocBas[nbc_id]*eindex+jj ];
    ctrl_x[jj] = local_pt_xyz[nbc_id][pos].x();
    ctrl_y[jj] = local_pt_xyz[nbc_id][pos].y();
    ctrl_z[jj] = local_pt_xyz[nbc_id][pos].z();
  }
}

void ALocal_InflowBC::get_SIEN( const int &nbc_id,
    const int &eindex, int * const &sien ) const
{
  for(int jj=0; jj<cell_nLocBas[nbc_id]; ++jj)
  {
    const int pos = local_cell_ien[nbc_id][ cell_nLocBas[nbc_id]*eindex+jj ];
    sien[jj] = local_node_pos[nbc_id][pos];
  }
}

std::vector<int> ALocal_InflowBC::get_SIEN( const int &nbc_id,
    const int &eindex ) const
{
  std::vector<int> out( cell_nLocBas[nbc_id], 0 );
  for(int jj=0; jj<cell_nLocBas[nbc_id]; ++jj)
  {
    const int pos = local_cell_ien[nbc_id][ cell_nLocBas[nbc_id]*eindex+jj ];
    out[jj] = local_node_pos[nbc_id][pos];
  }
  return out;
}

// EOF
