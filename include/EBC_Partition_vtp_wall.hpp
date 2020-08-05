#ifndef EBC_PARTITION_VTP_WALL_HPP
#define EBC_PARTITION_VTP_WALL_HPP
// ==================================================================
// EBC_Partition_vtp_wall.hpp
//
// Element boundary condition partition for the vessel wall
// in the coupled momentum method (CMM).
// 
// If the partition owns the wall domain, we will record the  
// radius, thickness, and Young's modulus.
//
// These information will be used for the wall coupling in CMM.
// 
// Author: Ju Liu
// Date: Jul. 25 2020
// ==================================================================
#include "EBC_Partition_vtp.hpp"

class EBC_Partition_vtp_wall : public EBC_Partition_vtp
{
  public:
    // The input ElemBC should be ElemBC_3D_tet_wall
    EBC_Partition_vtp_wall( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const ElemBC * const &ebc );

    virtual ~EBC_Partition_vtp_wall();

    // write the data to hdf5 file in folder /ebc/ebc_wall 
    virtual void write_hdf5( const char * FileName ) const;

    // write the data to hdf5 file in folder /GroupName/ebc_wall 
    virtual void write_hdf5( const char * FileName, const char * GroupName ) const;

  protected:
    // fluid density used in Young's modulus calculation
    const double fluid_density;

    // Length is num_local_node[0] 
    std::vector<double> part_thickness;
    std::vector<double> part_youngsmod;
};

#endif
