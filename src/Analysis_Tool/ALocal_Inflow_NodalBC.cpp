#include "ALocal_Inflow_NodalBC.hpp"

ALocal_Inflow_NodalBC::ALocal_Inflow_NodalBC( 
    const std::string &fileBaseName, const int &cpu_rank )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/inflow");

  outnormal = h5r -> read_Vector_3( gname.c_str(), "Outward_normal_vector" );
  
  act_area = h5r -> read_doubleScalar( gname.c_str(), "Inflow_active_area");
  ful_area = h5r -> read_doubleScalar( gname.c_str(), "Inflow_full_area");

  Num_LD = h5r -> read_intScalar(gname.c_str(), "Num_LD");

  // If this sub-domain contains local inflow bc points,
  // load the LDN array and the additional geometry info.
  if( Num_LD > 0 )
  {
    LDN = h5r->read_intVector( gname.c_str(), "LDN" );
    num_out_bc_pts = h5r->read_intScalar( gname.c_str(), "num_out_bc_pts" );
    centroid = h5r->read_Vector_3( gname.c_str(), "centroid" );
    outline_pts = h5r->read_doubleVector( gname.c_str(), "outline_pts" );
  }
  else
  {
    num_out_bc_pts = 0;
  }

  num_local_cell = h5r->read_intScalar( gname.c_str(), "num_local_cell" );
  num_local_node = h5r->read_intScalar( gname.c_str(), "num_local_node" );
  cell_nLocBas   = h5r->read_intScalar( gname.c_str(), "cell_nLocBas" );
 
  // If this partitioned sub-domain contains inlet surface element,
  // load its geometrical info 
  if(num_local_cell > 0)
  {
    local_pt_xyz = h5r->read_doubleVector( gname.c_str(), "local_pt_xyz" );
    local_tri_ien = h5r->read_intVector( gname.c_str(), "local_tri_ien" );
    local_node_pos = h5r->read_intVector( gname.c_str(), "local_node_pos" );
  }
  else
  {
    local_pt_xyz.clear();
    local_tri_ien.clear();
    local_node_pos.clear();
  }
  
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

double ALocal_Inflow_NodalBC::get_radius( const double &x, const double &y,
    const double &z ) const
{
  // num_out_bc_pts is set to be zero for parition that does not contain
  // inflow boundary points (i.e. Num_LD = 0 ).
  SYS_T::print_fatal_if( num_out_bc_pts == 0, "Error: ALocal_Inflow_NodalBC::get_radius, this function can only be called in sub-domains which contains the inflow boundary node.\n");

  const double rc = MATH_T::norm2( x-centroid(0), y-centroid(1), z-centroid(2) );

  // Now loop over the boundary points to find rb.
  double rb = MATH_T::norm2(x-outline_pts[0], y-outline_pts[1], z-outline_pts[2]);

  for(int ii=1; ii<num_out_bc_pts; ++ii)
  {
    double newdist = MATH_T::norm2(x-outline_pts[3*ii], y-outline_pts[3*ii+1], z-outline_pts[3*ii+2]);

    if(newdist < rb) rb = newdist;
  }

  return rc / (rb + rc);
}

double ALocal_Inflow_NodalBC::get_radius( const Vector_3 &pt ) const
{
  // num_out_bc_pts is set to be zero for parition that does not contain
  // inflow boundary points (i.e. Num_LD = 0 ).
  SYS_T::print_fatal_if( num_out_bc_pts == 0, "Error: ALocal_Inflow_NodalBC::get_radius, this function can only be called in sub-domains which contains the inflow boundary node.\n");

  const double x = pt.x();
  const double y = pt.y();
  const double z = pt.z();

  const double rc = MATH_T::norm2( x-centroid(0), y-centroid(1), z-centroid(2) );

  // Now loop over the boundary points to find rb.
  double rb = MATH_T::norm2(x-outline_pts[0], y-outline_pts[1], z-outline_pts[2]);

  for(int ii=1; ii<num_out_bc_pts; ++ii)
  {
    double newdist = MATH_T::norm2(x-outline_pts[3*ii], y-outline_pts[3*ii+1], z-outline_pts[3*ii+2]);

    if(newdist < rb) rb = newdist;
  }

  return rc / (rb + rc);
}

void ALocal_Inflow_NodalBC::get_ctrlPts_xyz( const int &eindex, 
    double * const &ctrl_x, double * const &ctrl_y, 
    double * const &ctrl_z ) const
{
  for(int jj=0; jj<cell_nLocBas; ++jj)
  {
    const int pos = local_tri_ien[ cell_nLocBas*eindex+jj ];
    ctrl_x[jj] = local_pt_xyz[3*pos];
    ctrl_y[jj] = local_pt_xyz[3*pos+1];
    ctrl_z[jj] = local_pt_xyz[3*pos+2];
  }
}

void ALocal_Inflow_NodalBC::get_SIEN( const int &eindex, int * const &sien ) const
{
  for(int jj=0; jj<cell_nLocBas; ++jj)
    sien[jj] = local_node_pos[ local_tri_ien[cell_nLocBas*eindex+jj] ];
}

// EOF
