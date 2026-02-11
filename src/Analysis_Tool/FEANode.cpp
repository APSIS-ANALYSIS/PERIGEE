#include "FEANode.hpp"

FEANode::FEANode( const std::string &fileBaseName, int cpu_rank )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  auto h5r = SYS_T::make_unique<HDF5_Reader>(file_id);

  ctrlPts_x = h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_x_loc");
  ctrlPts_y = h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_y_loc");
  ctrlPts_z = h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_z_loc");

  // Detect if the weights is in the h5 file, and read if yes
  if( H5Lexists(file_id, "/ctrlPts_loc/ctrlPts_w_loc", H5P_DEFAULT) )
    ctrlPts_w = h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_w_loc");
  else
    VEC_T::clean( ctrlPts_w );

  H5Fclose( file_id );
}

FEANode::FEANode( const HDF5_Reader * const &h5r )
{
  ctrlPts_x = h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_x_loc");
  ctrlPts_y = h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_y_loc");
  ctrlPts_z = h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_z_loc");

  // Detect if the weights is in the h5 file, and read if yes
  if( h5r -> check_data("/ctrlPts_loc/ctrlPts_w_loc") )
    ctrlPts_w = h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_w_loc");
  else
    VEC_T::clean( ctrlPts_w );
}

void FEANode::print_info() const
{
  std::cout<<"\n ctrlPts_x: \n";
  VEC_T::print(ctrlPts_x);
  std::cout<<"\n ctrlPts_y: \n";
  VEC_T::print(ctrlPts_y);
  std::cout<<"\n ctrlPts_z: \n";
  VEC_T::print(ctrlPts_z);
  std::cout<<"\n ctrlPts_w: \n";
  VEC_T::print(ctrlPts_w);
}

void FEANode::get_ctrlPts_xyz( 
    int num, const int * const &index,
    double * const &ctrl_x, double * const &ctrl_y, double * const &ctrl_z ) const
{
  for(int ii=0; ii<num; ++ii)
  {
    ctrl_x[ii] = ctrlPts_x[index[ii]];
    ctrl_y[ii] = ctrlPts_y[index[ii]];
    ctrl_z[ii] = ctrlPts_z[index[ii]];
  }
}

std::array<std::vector<double>, 3> FEANode::get_ctrlPts_xyz( 
    const std::vector<int> &index ) const
{
  // Allocate the size for the return object
  std::array<std::vector<double>, 3> out {{ std::vector<double>( index.size(), 0.0),
    std::vector<double>( index.size(), 0.0), 
    std::vector<double>( index.size(), 0.0) }};
  
  // Assign values for the return object
  for(unsigned int ii=0; ii<index.size(); ++ii)
  {
    out[0][ii] = ctrlPts_x[ index[ii] ];
    out[1][ii] = ctrlPts_y[ index[ii] ];
    out[2][ii] = ctrlPts_z[ index[ii] ];
  }
  return out;
}

void FEANode::get_ctrlPts_xyzw( 
    int num, const int * const &index,
    double * const &ctrl_x, double * const &ctrl_y, 
    double * const &ctrl_z, double * const &ctrl_w ) const
{
  for(int ii=0; ii<num; ++ii)
  {
    ctrl_x[ii] = ctrlPts_x[index[ii]];
    ctrl_y[ii] = ctrlPts_y[index[ii]];
    ctrl_z[ii] = ctrlPts_z[index[ii]];
    ctrl_w[ii] = ctrlPts_w[index[ii]];
  }
}

void FEANode::get_ctrlPts_xyw( 
    int num, const int * const &index,
    double * const &ctrl_x, double * const &ctrl_y, 
    double * const &ctrl_w ) const
{
  for(int ii=0; ii<num; ++ii)
  {
    ctrl_x[ii] = ctrlPts_x[index[ii]];
    ctrl_y[ii] = ctrlPts_y[index[ii]];
    ctrl_w[ii] = ctrlPts_w[index[ii]];
  }
}

void FEANode::get_ctrlPts_xy( 
    int num, const int * const &index,
    double * const &ctrl_x, double * const &ctrl_y ) const
{
  for(int ii=0; ii<num; ++ii)
  {
    ctrl_x[ii] = ctrlPts_x[index[ii]];
    ctrl_y[ii] = ctrlPts_y[index[ii]];
  }
}

double FEANode::get_memory_usage() const
{
  unsigned int total_length = 0;
  total_length += ctrlPts_x.size();
  total_length += ctrlPts_y.size();
  total_length += ctrlPts_z.size();
  total_length += ctrlPts_w.size();
  return double(total_length) * 8.0;
}

// EOF
