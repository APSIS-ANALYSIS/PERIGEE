#ifndef HDF5_GROUP_HPP
#define HDF5_GROUP_HPP
// ============================================================================
// HDF5_Group.hpp
//
// RAII wrapper for HDF5 group handles. The wrapped group is automatically
// closed with H5Gclose when the object leaves scope.
//
// Author: Ju Liu
// Date: Mar. 21 2026
// ============================================================================
#include "Sys_Tools.hpp"
#include "hdf5.h"

class HDF5_Group
{
  public:
    explicit HDF5_Group( hid_t in_group_id = -1 ) noexcept
      : group_id(in_group_id) {}

    ~HDF5_Group() noexcept
    {
      if( group_id >= 0 ) H5Gclose( group_id );
    }

    HDF5_Group( const HDF5_Group & ) = delete;
    HDF5_Group &operator=( const HDF5_Group & ) = delete;

    HDF5_Group( HDF5_Group &&other ) noexcept
      : group_id(other.group_id)
    {
      other.group_id = -1;
    }

    HDF5_Group &operator=( HDF5_Group &&other ) noexcept
    {
      if( this != &other )
      {
        if( group_id >= 0 ) H5Gclose( group_id );

        group_id = other.group_id;
        other.group_id = -1;
      }

      return *this;
    }

    static HDF5_Group create( hid_t location_id, const char * group_name )
    {
      const hid_t new_group_id = H5Gcreate( location_id, group_name,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

      check_group_id( new_group_id, "create", group_name );
      return HDF5_Group( new_group_id );
    }

    static HDF5_Group create( hid_t location_id, const std::string &group_name )
    {
      return create( location_id, group_name.c_str() );
    }

    static HDF5_Group open( hid_t location_id, const char * group_name )
    {
      const hid_t new_group_id = H5Gopen( location_id, group_name, H5P_DEFAULT );

      check_group_id( new_group_id, "open", group_name );
      return HDF5_Group( new_group_id );
    }

    static HDF5_Group open( hid_t location_id, const std::string &group_name )
    {
      return open( location_id, group_name.c_str() );
    }

    hid_t id() const noexcept { return group_id; }
    operator hid_t() const noexcept { return group_id; }

    bool is_valid() const noexcept { return group_id >= 0; }

  private:
    hid_t group_id;

    static void check_group_id( hid_t in_group_id,
        const char * operation, const char * group_name )
    {
      if( in_group_id < 0 )
      {
        std::ostringstream oss;
        oss<<"Error: HDF5_Group::"<<operation
          <<" cannot access group "<<group_name<<".\n";
        SYS_T::print_fatal( oss.str().c_str() );
      }
    }
};

#endif
