#ifndef HDF5_TOOLS_HPP
#define HDF5_TOOLS_HPP
// ============================================================================
// HDF5_Tools.hpp
//
// HDF5_T is a namespace that provides some convenient tools to write or read a
// hdf5 file.
//
// Author: Ju Liu
// Date: Dec 26 2021
// ============================================================================
#include "HDF5_Writer.hpp"
#include "HDF5_Reader.hpp"

namespace HDF5_T
{
  inline int read_intScalar( const std::string &filename, 
      const std::string groupname, const std::string &dataname )
  {
    hid_t file_id = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    HDF5_Reader * h5r = new HDF5_Reader( file_id );

    const int output = h5r -> read_intScalar( groupname.c_str(), dataname.c_str() );

    delete h5r; H5Fclose( file_id );

    return output;
  }
}

#endif
