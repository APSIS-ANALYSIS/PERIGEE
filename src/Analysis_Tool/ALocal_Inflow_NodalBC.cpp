#include "ALocal_Inflow_NodalBC.hpp"

ALocal_Inflow_NodalBC::ALocal_Inflow_NodalBC( 
    const std::string &fileBaseName, const int &cpu_rank )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/inflow");

  num_nbc = h5r -> read_intScalar( gname.c_str(), "num_nbc" );

  act_area.resize(num_nbc); ful_area.resize(num_nbc); num_out_bc_pts.resize(num_nbc);
  num_local_node.resize(num_nbc); num_local_cell.resize(num_nbc);
  cell_nLocBas.resize(num_nbc); outnormal.resize(num_nbc);
  Num_LD.resize(num_nbc); LDN.resize(num_nbc);
  centroid.resize(num_nbc); outline_pts.resize(num_nbc);
  local_pt_xyz.resize(num_nbc); local_tri_ien.resize(num_nbc); local_node_pos.resize(num_nbc);

  std::string groupbase(gname);
  groupbase.append("/nbcid_");

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    std::string subgroup_name(groupbase);
    subgroup_name.append( SYS_T::to_string(nbc_id) );

    act_area[nbc_id]       = h5r->read_doubleScalar( subgroup_name.c_str(), "Inflow_active_area" );
    ful_area[nbc_id]       = h5r->read_doubleScalar( subgroup_name.c_str(), "Inflow_full_area" );
    num_local_node[nbc_id] = h5r->read_intScalar(    subgroup_name.c_str(), "num_local_node" );
    num_local_cell[nbc_id] = h5r->read_intScalar(    subgroup_name.c_str(), "num_local_cell" );
    cell_nLocBas[nbc_id]   = h5r->read_intScalar(    subgroup_name.c_str(), "cell_nLocBas" );
    outnormal[nbc_id]      = h5r->read_Vector_3(     subgroup_name.c_str(), "Outward_normal_vector" );
    Num_LD[nbc_id]         = h5r -> read_intScalar(  subgroup_name.c_str(), "Num_LD" );

    // If this sub-domain contains local inflow bc points, load the LDN array.
    if( Num_LD[nbc_id] > 0 )
    {
      LDN[nbc_id]            = h5r->read_intVector( subgroup_name.c_str(), "LDN" );
      num_out_bc_pts[nbc_id] = h5r->read_intScalar(    subgroup_name.c_str(), "num_out_bc_pts" );
      centroid[nbc_id]       = h5r->read_Vector_3(     subgroup_name.c_str(), "centroid" );
      outline_pts[nbc_id]    = h5r->read_doubleVector( subgroup_name.c_str(), "outline_pts" );
    }
    else
    {
      num_out_bc_pts[nbc_id] = 0;

      LDN[nbc_id].clear(); num_out_bc_pts[nbc_id].clear();
      centroid[nbc_id].clear(); outline_pts[nbc_id].clear();
    }
 
    // If this partitioned sub-domain contains inlet surface element,
    // load its geometrical info 
    if(num_local_cell[nbc_id] > 0)
    {

      local_pt_xyz[nbc_id]   = h5r->read_doubleVector( subgroup_name.c_str(), "local_pt_xyz" );
      local_tri_ien[nbc_id]  = h5r->read_intVector(    subgroup_name.c_str(), "local_tri_ien" );
      local_node_pos[nbc_id] = h5r->read_intVector(    subgroup_name.c_str(), "local_node_pos" );
    }
    else
    {
      local_pt_xyz[nbc_id].clear();   local_tri_ien[nbc_id].clear();
      local_node_pos[nbc_id].clear();
    }
  } // end nbc_id-loop
  
  delete h5r; H5Fclose( file_id );
}

ALocal_Inflow_NodalBC::~ALocal_Inflow_NodalBC()
{
  VEC_T::clean(LDN);
  VEC_T::clean(outline_pts);
  VEC_T::clean(local_pt_xyz);
  VEC_T::clean(local_tri_ien);
  VEC_T::clean(local_node_pos);
}

double ALocal_Inflow_NodalBC::get_radius( const int &nbc_id,
    const Vector_3 &pt ) const
{
  // num_out_bc_pts is set to be zero for parition that does not contain
  // inflow boundary points (i.e. Num_LD = 0 ).
  SYS_T::print_fatal_if( num_out_bc_pts[nbc_id] == 0, "Error: ALocal_Inflow_NodalBC::get_radius, this function can only be called in sub-domains which contains the inflow boundary node.\n");

  const double x = pt.x();
  const double y = pt.y();
  const double z = pt.z();

  const double rc = MATH_T::norm2( x-centroid[nbc_id](0), y-centroid[nbc_id](1),
      z-centroid[nbc_id](2) );

  // Now loop over the boundary points to find rb.
  double rb = MATH_T::norm2(x-outline_pts[nbc_id][0], y-outline_pts[nbc_id][1], 
      z-outline_pts[nbc_id][2]);

  for(int ii=1; ii<num_out_bc_pts[nbc_id]; ++ii)
  {
    double newdist = MATH_T::norm2(x-outline_pts[nbc_id][3*ii],
      y-outline_pts[nbc_id][3*ii+1], z-outline_pts[nbc_id][3*ii+2]);

    if(newdist < rb) rb = newdist;
  }

  return rc / (rb + rc);
}

void ALocal_Inflow_NodalBC::get_ctrlPts_xyz( const int &nbc_id,
    const int &eindex, double * const &ctrl_x, double * const &ctrl_y, 
    double * const &ctrl_z ) const
{
  for(int jj=0; jj<cell_nLocBas[nbc_id]; ++jj)
  {
    const int pos = local_tri_ien[nbc_id][ cell_nLocBas[nbc_id]*eindex+jj ];
    ctrl_x[jj] = local_pt_xyz[nbc_id][3*pos+0];
    ctrl_y[jj] = local_pt_xyz[nbc_id][3*pos+1];
    ctrl_z[jj] = local_pt_xyz[nbc_id][3*pos+2];
  }
}

void ALocal_Inflow_NodalBC::get_SIEN( const int &nbc_id,
    const int &eindex, int * const &sien ) const
{
  for(int jj=0; jj<cell_nLocBas[nbc_id]; ++jj)
  {
    const int pos = local_tri_ien[nbc_id][ cell_nLocBas[nbc_id]*eindex+jj ];
    sien[jj] = local_node_pos[nbc_id][pos];
  }
}

// EOF
