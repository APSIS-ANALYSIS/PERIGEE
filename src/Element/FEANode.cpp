#include "FEANode.hpp"

FEANode::FEANode( const std::string &fileBaseName, const int &cpu_rank )
{
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_x_loc", ctrlPts_x);
  h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_y_loc", ctrlPts_y);
  h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_z_loc", ctrlPts_z);

  // Detect if the weights is in the h5 file, and read if yes
  if( H5Lexists(file_id, "/ctrlPts_loc/ctrlPts_w_loc", H5P_DEFAULT) )
    h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_w_loc", ctrlPts_w);
  else
    VEC_T::clean( ctrlPts_w );

  delete h5r; H5Fclose( file_id );
}


FEANode::~FEANode()
{}


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
    const int &num, const int * const &index,
    double * const &ctrl_x, double * const &ctrl_y, double * const &ctrl_z ) const
{
  int pos;
  for(int ii=0; ii<num; ++ii)
  {
    pos = index[ii];
    ctrl_x[ii] = ctrlPts_x[pos];
    ctrl_y[ii] = ctrlPts_y[pos];
    ctrl_z[ii] = ctrlPts_z[pos];
  }
}


void FEANode::get_ctrlPts_xyzw( 
    const int &num, const int * const &index,
    double * const &ctrl_x, double * const &ctrl_y, 
    double * const &ctrl_z, double * const &ctrl_w ) const
{
  int pos;
  for(int ii=0; ii<num; ++ii)
  {
    pos = index[ii];
    ctrl_x[ii] = ctrlPts_x[pos];
    ctrl_y[ii] = ctrlPts_y[pos];
    ctrl_z[ii] = ctrlPts_z[pos];
    ctrl_w[ii] = ctrlPts_w[pos];
  }
}


void FEANode::get_ctrlPts_xyw( 
    const int &num, const int * const &index,
    double * const &ctrl_x, double * const &ctrl_y, 
    double * const &ctrl_w ) const
{
  int pos;
  for(int ii=0; ii<num; ++ii)
  {
    pos = index[ii];
    ctrl_x[ii] = ctrlPts_x[pos];
    ctrl_y[ii] = ctrlPts_y[pos];
    ctrl_w[ii] = ctrlPts_w[pos];
  }
}


void FEANode::get_ctrlPts_xy( 
    const int &num, const int * const &index,
    double * const &ctrl_x, double * const &ctrl_y ) const
{
  int pos;
  for(int ii=0; ii<num; ++ii)
  {
    pos = index[ii];
    ctrl_x[ii] = ctrlPts_x[pos];
    ctrl_y[ii] = ctrlPts_y[pos];
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
