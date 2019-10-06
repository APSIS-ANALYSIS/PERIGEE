#include "ALocal_Inflow_NodalBC.hpp"

ALocal_Inflow_NodalBC::ALocal_Inflow_NodalBC( 
    const std::string &fileBaseName, const int &cpu_rank )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/inflow");
  
  h5r -> read_doubleVector( gname.c_str(), "Outward_normal_vector", outvec );

  act_area = h5r -> read_doubleScalar( gname.c_str(), "Inflow_active_area");
  ful_area = h5r -> read_doubleScalar( gname.c_str(), "Inflow_full_area");

  Num_LD = h5r -> read_intScalar(gname.c_str(), "Num_LD");

  // If this sub-domain contains local inflow bc points,
  // load the LDN array and the additional geometry info.
  if( Num_LD > 0 )
  {
    h5r->read_intVector( gname.c_str(), "LDN", LDN );
    num_out_bc_pts = h5r->read_intScalar( gname.c_str(), "num_out_bc_pts");
    h5r->read_doubleVector( gname.c_str(), "centroid", centroid);
    h5r->read_doubleVector( gname.c_str(), "outline_pts", outline_pts );
  }
  else
  {
    num_out_bc_pts = 0;
  }

  delete h5r; H5Fclose( file_id );
}

ALocal_Inflow_NodalBC::~ALocal_Inflow_NodalBC()
{
  VEC_T::clean(LDN);
  VEC_T::clean(outvec);
  VEC_T::clean(centroid);
  VEC_T::clean(outline_pts);
}


double ALocal_Inflow_NodalBC::get_radius( const double &x, const double &y,
    const double &z ) const
{
  // num_out_bc_pts is set to be zero for parition that does not contain
  // inflow boundary points (i.e. Num_LD = 0 ).
  SYS_T::print_fatal_if( num_out_bc_pts == 0, "Error: ALocal_Inflow_NodalBC::get_radius, this function can only be called in sub-domains which contains the inflow boundary node.\n");

  const double rc = MATH_T::norm2(x-centroid[0], y-centroid[1], z-centroid[2]);

  // Now loop over the boundary points to find rb.
  double rb = MATH_T::norm2(x-outline_pts[0], y-outline_pts[1], z-outline_pts[2]);

  for(int ii=1; ii<num_out_bc_pts; ++ii)
  {
    double newdist = MATH_T::norm2(x-outline_pts[3*ii], y-outline_pts[3*ii+1], z-outline_pts[3*ii+2]);

    if(newdist < rb) rb = newdist;
  }

  return rc / (rb + rc);
}

// EOF
