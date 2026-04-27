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
#include "HDF5_Reader.hpp"

namespace HDF5_T
{
  inline int read_intScalar( const std::string &filename,
      const std::string &groupname, const std::string &dataname )
  {
    std::unique_ptr<HDF5_Reader> h5r = SYS_T::make_unique<HDF5_Reader>( filename );

    const int output = h5r -> read_intScalar( groupname.c_str(), dataname.c_str() );

    return output;
  }

  inline double read_doubleScalar( const std::string &filename,
      const std::string &groupname, const std::string &dataname )
  {
    std::unique_ptr<HDF5_Reader> h5r = SYS_T::make_unique<HDF5_Reader>( filename );

    const double output = h5r -> read_doubleScalar( groupname.c_str(), dataname.c_str() );


    return output;
  }

  inline std::vector<int> read_intVector( const std::string &filename,
      const std::string &groupname, const std::string &dataname )
  {
    std::unique_ptr<HDF5_Reader> h5r = SYS_T::make_unique<HDF5_Reader>( filename );

    std::vector<int> output = h5r -> read_intVector( groupname.c_str(), dataname.c_str() );

    output.shrink_to_fit();

    return output;
  }
}

#endif
