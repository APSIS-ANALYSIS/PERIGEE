#include "HDF5_Writer.hpp"

void HDF5_Writer::write_intScalar( const hid_t &group_id,
    const char * const &data_name, const int &value ) const
{
  write_intScalar_impl( group_id, data_name, value );
}

void HDF5_Writer::write_intScalar( const char * const &data_name,
    const int &value ) const
{
  write_intScalar_impl( file_id, data_name, value );
}

void HDF5_Writer::write_int64Scalar( const char * const &data_name,
    const int64_t &value ) const
{
  hid_t dataspace, dataset;
  hsize_t dims[1];
  herr_t status;
  dims[0] = 1;

  dataspace = H5Screate_simple(1, dims, NULL);
  dataset   = H5Dcreate( file_id, data_name, H5T_NATIVE_LLONG,
     dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

  status = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &value );

  check_error(status, "write_int64Scalar");

  H5Dclose( dataset );
  H5Sclose( dataspace );
}

void HDF5_Writer::write_uintScalar( 
    const hid_t &group_id, const char * const &data_name,
    const unsigned int &value ) const
{
  hid_t dataspace, dataset;
  hsize_t dims[1];
  herr_t status;
  dims[0] = 1;

  dataspace = H5Screate_simple(1, dims, NULL);
  dataset   = H5Dcreate( group_id, data_name, H5T_NATIVE_UINT,
     dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

  status = H5Dwrite( dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &value );

  check_error(status, "write_uintScalar");

  H5Dclose( dataset );
  H5Sclose( dataspace );
}

void HDF5_Writer::write_doubleScalar( const hid_t &group_id,
    const char * const &data_name, const double &value ) const
{
  write_doubleScalar_impl( group_id, data_name, value );
}

void HDF5_Writer::write_doubleScalar( const char * const &data_name,
    const double &value ) const
{
  write_doubleScalar_impl( file_id, data_name, value );
}

void HDF5_Writer::write_intVector( const hid_t &group_id,
    const char * const &data_name, const std::vector<int> &value ) const
{
  write_intVector_impl( group_id, data_name, value.data(), VEC_T::get_size( value ) );
}

void HDF5_Writer::write_intVector( const char * const &data_name, 
    const std::vector<int> &value ) const
{
  write_intVector_impl( file_id, data_name, value.data(), VEC_T::get_size( value ) );
}

void HDF5_Writer::write_uintVector( const hid_t &group_id,
    const char * const &data_name, const std::vector<unsigned int> &value ) const
{
  hsize_t dims[1]; dims[0] = value.size();
  if(dims[0] > 0)
  {
    hid_t dataspace, dataset;
    herr_t status;

    dataspace = H5Screate_simple(1, dims, NULL);
    dataset   = H5Dcreate( group_id, data_name, H5T_NATIVE_UINT,
        dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    status = H5Dwrite( dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
       H5P_DEFAULT, &value[0] );
   
    check_error(status, "write_uintVector");

    H5Dclose( dataset );
    H5Sclose( dataspace );
  }
}

void HDF5_Writer::write_intVector( const hid_t &group_id,
    const char * const &data_name, const int * const &value, 
    const int &length ) const
{
  write_intVector_impl( group_id, data_name, value, length );
}

void HDF5_Writer::write_intVector( const char * const &data_name, 
    const int * const  &value,
    const int &length ) const
{
  write_intVector_impl( file_id, data_name, value, length );
}

void HDF5_Writer::write_int64Vector( const hid_t &group_id,
    const char * const &data_name, 
    const int64_t * const &value, const int64_t &length ) const
{
  hsize_t dims[1]; dims[0] = length;
  if(dims[0] > 0)
  {
    hid_t dataspace, dataset;
    herr_t status;

    dataspace = H5Screate_simple(1, dims, NULL);
    dataset   = H5Dcreate( group_id, data_name, H5T_NATIVE_LLONG,
        dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    status = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, value );

    check_error(status, "write_int64Vector");
    H5Dclose( dataset );
    H5Sclose( dataspace );
  }
}

void HDF5_Writer::write_int64Vector( const char * const &data_name, 
    const int64_t * const  &value,
    const int64_t &length ) const
{
  hsize_t dims[1]; dims[0] = length;
  if(dims[0] > 0)
  {
    hid_t dataspace, dataset;
    herr_t status;

    dataspace = H5Screate_simple(1, dims, NULL);
    dataset   = H5Dcreate( file_id, data_name, H5T_NATIVE_LLONG,
        dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    status = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, value );

    check_error(status, "write_int64Vector");

    H5Dclose( dataset );
    H5Sclose( dataspace );
  }
}

void HDF5_Writer::write_doubleVector( const hid_t &group_id,
    const char * const &data_name,
    const double * const &value, const int &length ) const
{
  write_doubleVector_impl( group_id, data_name, value, length );
}

void HDF5_Writer::write_doubleVector( const char * const &data_name,
    const double * const &value, const int &length ) const
{
  write_doubleVector_impl( file_id, data_name, value, length );
}

void HDF5_Writer::write_doubleVector( const hid_t &group_id,
    const char * const &data_name, const std::vector<double> &value ) const
{
  write_doubleVector_impl( group_id, data_name, value.data(), value.size() );
}

void HDF5_Writer::write_doubleVector( const char * const &data_name, 
    const std::vector<double> &value ) const
{
  write_doubleVector_impl( file_id, data_name, value.data(), value.size() );
}

void HDF5_Writer::write_Vector_3( const hid_t &group_id, const char * const &data_name,
    const Vector_3 &value ) const
{
  // First convert Vector_3 to a double array
  const double val[3] = { value.x(), value.y(), value.z() };

  hsize_t dims[1] = { 3 };
  hid_t dataspace = H5Screate_simple(1, dims, NULL);
  hid_t dataset   = H5Dcreate( group_id, data_name, H5T_NATIVE_DOUBLE,
      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  herr_t status = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &val[0] );

  check_error(status, "write_Vector_3");

  H5Dclose( dataset );
  H5Sclose( dataspace );
}

void HDF5_Writer::write_Vector_3( const char * const &data_name, const Vector_3 &value ) const
{
  // First convert Vector_3 to a double array
  const double val[3] = { value.x(), value.y(), value.z() };

  hsize_t dims[1] = { 3 };
  hid_t dataspace = H5Screate_simple(1, dims, NULL);
  hid_t dataset   = H5Dcreate( file_id, data_name, H5T_NATIVE_DOUBLE,
      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  herr_t status = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &val[0] );

  check_error(status, "write_Vector_3");

  H5Dclose( dataset );
  H5Sclose( dataspace );
}

void HDF5_Writer::write_Tensor2_3D( const hid_t &group_id, const char * const &data_name,
    const Tensor2_3D &value ) const
{
  // First convert Tensor2_3D to a double array
  const double val[9] = { value(0), value(1), value(2), value(3), value(4), value(5), value(6), value(7), value(8) };

  hsize_t dims[1] = { 9 };
  hid_t dataspace = H5Screate_simple(1, dims, NULL);
  hid_t dataset   = H5Dcreate( group_id, data_name, H5T_NATIVE_DOUBLE,
      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  herr_t status = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &val[0] );

  check_error(status, "write_Tensor2_3D");

  H5Dclose( dataset );
  H5Sclose( dataspace );
}

void HDF5_Writer::write_Tensor2_3D( const char * const &data_name, const Tensor2_3D &value ) const
{
  // First convert Tensor2_3D to a double array
  const double val[9] = { value(0), value(1), value(2), value(3), value(4), value(5), value(6), value(7), value(8) };

  hsize_t dims[1] = { 9 };
  hid_t dataspace = H5Screate_simple(1, dims, NULL);
  hid_t dataset   = H5Dcreate( file_id, data_name, H5T_NATIVE_DOUBLE,
      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  herr_t status = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &val[0] );

  check_error(status, "write_Tensor2_3D");

  H5Dclose( dataset );
  H5Sclose( dataspace );
}

void HDF5_Writer::write_intMatrix( const hid_t &group_id,
    const char * const &data_name, const std::vector<int> &value,
    const int &row_num, const int &col_num ) const
{
  hsize_t dims[2]; dims[0] = row_num; dims[1] = col_num;
  if( int(value.size()) != row_num * col_num )
  {
    std::cerr<<"ERROR: the matrix size is incompatible with the given row/column size. \n";
    exit(1);
  }

  if(value.size() > 0)
  {
    hid_t dataspace, dataset;
    herr_t status;

    dataspace = H5Screate_simple(2, dims, NULL);
    dataset = H5Dcreate( group_id, data_name, H5T_NATIVE_INT, dataspace,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    status = H5Dwrite( dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, &value[0] );

    check_error(status, "write_intMatrix");

    H5Dclose( dataset );
    H5Sclose( dataspace );
  }
}

void HDF5_Writer::write_doubleMatrix( const hid_t &group_id,
    const char * const &data_name, const std::vector<double> &value,
    const int &row_num, const int &col_num ) const
{
  hsize_t dims[2]; dims[0] = row_num; dims[1] = col_num;
  if( int(value.size()) != row_num * col_num )
  {
    std::cerr<<"ERROR: the matrix size is incompatible with the given row/column size. \n";
    exit(1);
  }

  if(value.size() > 0)
  {
    hid_t dataspace, dataset;
    herr_t status;
    dataspace = H5Screate_simple(2, dims, NULL);
    dataset = H5Dcreate( group_id, data_name, H5T_NATIVE_DOUBLE, dataspace,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    status = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, &value[0] );

    check_error(status, "write_doubleMatrix");

    H5Dclose( dataset );
    H5Sclose( dataspace );
  }
}

void HDF5_Writer::write_string( const char * const &data_name, 
    const std::string &string_input ) const
{
  write_string_impl(file_id, data_name, string_input);
}

void HDF5_Writer::write_string( const hid_t &group_id,
    const char * const &data_name, 
    const std::string &string_input ) const
{
  write_string_impl(group_id, data_name, string_input);
}

void HDF5_Writer::write_string_impl(hid_t location_id, 
    const char * const &data_name,
    const std::string& string_input ) const
{
  hsize_t dims[1] = { 1 };

  hid_t datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, string_input.size());

  hid_t dataspace = H5Screate_simple(1, dims, NULL);

  hid_t dataset = H5Dcreate(location_id, data_name, datatype, 
      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, string_input.data() );

  H5Dclose( dataset );
  H5Sclose( dataspace );
  H5Tclose( datatype );
}

void HDF5_Writer::write_intScalar_impl( hid_t location_id, 
    const char * const &data_name, 
    const int &value ) const
{
  hsize_t dims[1] = { 1 };

  hid_t dataspace = H5Screate_simple(1, dims, NULL);
  hid_t dataset   = H5Dcreate( location_id, data_name, H5T_NATIVE_INT,
      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

  herr_t status = H5Dwrite( dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &value );

  check_error(status, "write_intScalar");

  H5Dclose( dataset );
  H5Sclose( dataspace );
}

void HDF5_Writer::write_doubleScalar_impl( hid_t location_id, 
    const char * const &data_name, 
    const double &value ) const
{
  hsize_t dims[1] = { 1 };

  hid_t dataspace = H5Screate_simple(1, dims, NULL);
  hid_t dataset   = H5Dcreate( location_id, data_name, H5T_NATIVE_DOUBLE,
      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  herr_t status = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &value );

  check_error(status, "write_doubleScalar");

  H5Dclose( dataset );
  H5Sclose( dataspace );
}

void HDF5_Writer::write_intVector_impl( hid_t location_id,
    const char * const &data_name, 
    const int * const &value, const int &length ) const
{
  hsize_t dims[1] = { length };
  if(dims[0] > 0)
  {
    hid_t dataspace = H5Screate_simple(1, dims, NULL);
    hid_t dataset   = H5Dcreate( location_id, data_name, H5T_NATIVE_INT,
        dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    herr_t status = H5Dwrite( dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, &value[0] );
    
    check_error(status, "write_intVector");
    
    H5Dclose( dataset );
    H5Sclose( dataspace );
  }
}

void HDF5_Writer::write_doubleVector_impl( hid_t location_id,
    const char * const &data_name,
    const double * const &value, const int &length ) const
{
  hsize_t dims[1] = { length };
  if(dims[0] > 0)
  {
    hid_t dataspace = H5Screate_simple(1, dims, NULL);
    hid_t dataset   = H5Dcreate( location_id, data_name, H5T_NATIVE_DOUBLE,
        dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    herr_t status = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, value );

    check_error(status, "write_doubleVector");
    H5Dclose( dataset );
    H5Sclose( dataspace );
  }    
}

// EOF
