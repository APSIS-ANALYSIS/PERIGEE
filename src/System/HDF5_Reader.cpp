#include "HDF5_Reader.hpp"

HDF5_Reader::HDF5_Reader( const hid_t &in_file_id )
: file_id(in_file_id)
{}

HDF5_Reader::~HDF5_Reader()
{}

bool HDF5_Reader::check_data( const char * const &name ) const
{
  return H5Lexists(file_id, name, H5P_DEFAULT);
}

void HDF5_Reader::read_intArray(const char * const &group_name,
    const char * const &data_name,
    hid_t &data_rank, hsize_t * &data_dims, int * &data  ) const
{
  // open group file and data file
  hid_t group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);
  hid_t data_id = H5Dopen(group_id, data_name, H5P_DEFAULT);

  // retrive dataspace of the dataset
  hid_t data_space = H5Dget_space( data_id );

  // get the dimension of the array in data_space
  data_rank = H5Sget_simple_extent_ndims( data_space );

  // allocating the array that store the size in each dimension
  data_dims = new hsize_t [data_rank];

  // get the length in each direction
  herr_t status = H5Sget_simple_extent_dims( data_space ,data_dims, NULL );

  check_error(status, "read_intArray");

  // define the memory allocation for reading
  hid_t mem_space = H5Screate_simple(data_rank, data_dims, NULL);

  // setup the buffer
  hsize_t dSize = 1;

  for(hid_t ii=0; ii<data_rank; ++ii)
    dSize = dSize * data_dims[ii];

  // the data will pass out in a one-dimensional array
  data = new int [dSize];

  status = H5Dread( data_id, H5T_NATIVE_INT, mem_space, data_space,
      H5P_DEFAULT, data );

  check_error(status, "read_intArray");

  // clean out
  H5Sclose( mem_space );
  H5Sclose( data_space );
  H5Dclose( data_id );
  H5Gclose( group_id );
}

void HDF5_Reader::read_doubleArray(const char * const &group_name,
    const char * const &data_name,
    hid_t &data_rank, hsize_t * &data_dims, double * &data  ) const
{
  // open group file and data file
  hid_t group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);
  hid_t data_id  = H5Dopen(group_id, data_name, H5P_DEFAULT);

  // retrive dataspace of the dataset
  hid_t data_space = H5Dget_space( data_id );

  // get the dimension of the array in data_space
  data_rank = H5Sget_simple_extent_ndims( data_space );

  // allocating the array the store the size in each dimension
  data_dims = new hsize_t [data_rank];

  // get the length in each direction
  herr_t status = H5Sget_simple_extent_dims(data_space ,data_dims, NULL);

  check_error(status, "read_doubleArray : H5Sget_simple_extent_dims");

  // define the memory allocation for reading
  hid_t mem_space = H5Screate_simple(data_rank, data_dims, NULL);

  hsize_t dSize = 1;

  for(hid_t ii=0; ii<data_rank; ++ii)
    dSize = dSize * data_dims[ii];

  data = new double [dSize];

  status = H5Dread(data_id, H5T_NATIVE_DOUBLE, mem_space, data_space, H5P_DEFAULT, data);

  check_error(status, "read_doubleArray : H5Dread");

  // clean out
  H5Sclose( mem_space );
  H5Sclose( data_space );
  H5Dclose( data_id );
  H5Gclose( group_id );
}

int HDF5_Reader::read_intScalar( const char * const &group_name,
    const char * const &data_name ) const
{
  hid_t drank;
  hsize_t * ddims;
  int * intdata;

  read_intArray( group_name, data_name, drank, ddims, intdata );

  if( drank != 1 || ddims[0] != 1 )
  {
    std::ostringstream oss;
    oss<<"Error: HDF5_Reader::read_intScalar read data at "<<group_name;
    oss<<" with name "<<data_name<<" is not a scalar! \n";
    SYS_T::print_fatal( oss.str().c_str() );
  }

  int outdata = intdata[0];
  delete [] ddims; delete [] intdata;

  return outdata;
}

double HDF5_Reader::read_doubleScalar( const char * const &group_name,
    const char * const &data_name ) const
{
  hid_t drank;
  hsize_t * ddims;
  double * ddata;

  read_doubleArray( group_name, data_name, drank, ddims, ddata );

  if( drank != 1 || ddims[0] != 1 )
  {
    std::ostringstream oss;
    oss<<"Error: HDF5_Reader::read_doubleScalar read data at "<<group_name;
    oss<<" with name "<<data_name<<" is not a scalar! \n";
    SYS_T::print_fatal( oss.str().c_str() );
  }

  double outdata = ddata[0];
  delete [] ddims; delete [] ddata;

  return outdata;
}

void HDF5_Reader::read_intVector( const char * const &group_name,
    const char * const &data_name,
    std::vector<int> &out ) const
{
  hid_t drank;
  hsize_t * ddims;
  int * intdata;

  read_intArray( group_name, data_name, drank, ddims, intdata );

  if( drank != 1 )
  {
    std::ostringstream oss;
    oss<<"Error: HDF5_Reader::read_intVector read data at "<<group_name;
    oss<<" with name "<<data_name<<" is not a 1D vector! \n";
    SYS_T::print_fatal( oss.str().c_str() );
  }

  VEC_T::fillArray( out, intdata, ddims[0] );

  delete [] ddims; delete [] intdata;
}

void HDF5_Reader::read_doubleVector( const char * const &group_name,
    const char * const &data_name,
    std::vector<double> &out ) const
{
  hid_t drank;
  hsize_t * ddims;
  double * ddata;

  read_doubleArray( group_name, data_name, drank, ddims, ddata );

  if( drank != 1 )
  {
    std::ostringstream oss;
    oss<<"Error: HDF5_Reader::read_doubleVector read data at "<<group_name;
    oss<<" with name "<<data_name<<" is not a 1D vector! \n";
    SYS_T::print_fatal( oss.str().c_str() );
  }

  VEC_T::fillArray( out, ddata, ddims[0] );

  delete [] ddims; delete [] ddata;
}

void HDF5_Reader::read_intMatrix( const char * const &group_name,
    const char * const &data_name,
    std::vector<int> &out, int &num_row, int &num_col ) const
{
  hid_t drank;
  hsize_t * ddims;
  int * intdata;

  read_intArray( group_name, data_name, drank, ddims, intdata );

  if( drank != 2 )
  {
    std::ostringstream oss;
    oss<<"Error: HDF5_Reader::read_intMatrix read data at "<<group_name;
    oss<<" with name "<<data_name<<" is not a 2D matrix! \n";
    SYS_T::print_fatal( oss.str().c_str() );
  }

  num_row = ddims[0];
  num_col = ddims[1];

  VEC_T::fillArray( out, intdata, num_row * num_col );

  delete [] ddims; delete [] intdata;
}

void HDF5_Reader::read_doubleMatrix( const char * const &group_name,
    const char * const &data_name,
    std::vector<double> &out, int &num_row, int &num_col ) const
{
  hid_t drank;
  hsize_t * ddims;
  double * ddata;

  read_doubleArray( group_name, data_name, drank, ddims, ddata );

  if( drank != 2 )
  {
    std::ostringstream oss;
    oss<<"Error: HDF5_Reader::read_doubleMatrix read data at "<<group_name;
    oss<<" with name "<<data_name<<" is not a 2D matrix! \n";
    SYS_T::print_fatal( oss.str().c_str() );
  }

  num_row = ddims[0];
  num_col = ddims[1];

  VEC_T::fillArray( out, ddata, num_row * num_col );

  delete [] ddims; delete [] ddata;
}

void HDF5_Reader::read_string( const char * const &group_name,
    const char * const &data_name,
    std::string &string_out ) const
{
  // open group file and data file
  hid_t group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);
  hid_t data_id  = H5Dopen(group_id, data_name, H5P_DEFAULT);

  hid_t filetype = H5Dget_type( data_id );
  size_t sdim = H5Tget_size( filetype );
  sdim++;

  hid_t data_space = H5Dget_space( data_id );

  hid_t data_rank = H5Sget_simple_extent_ndims( data_space );

  SYS_T::print_fatal_if(data_rank!=1,
      "Error: HDF5_Reader::read_string file format is wrong.\n");

  hsize_t data_dims[1];

  H5Sget_simple_extent_dims( data_space ,data_dims, NULL );

  char * rdata = new char [sdim * data_dims[0]];

  hid_t memtype = H5Tcopy(H5T_C_S1);
  H5Tset_size(memtype, sdim);

  H5Dread(data_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);

  string_out.assign(rdata);

  delete [] rdata; rdata = NULL;
  H5Tclose( memtype );
  H5Tclose( filetype );
  H5Sclose( data_space );
  H5Dclose( data_id );
  H5Gclose( group_id );
}

// EOF
