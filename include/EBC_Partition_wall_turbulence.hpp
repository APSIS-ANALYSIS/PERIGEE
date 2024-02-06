#ifndef EBC_PARTITION_WALL_TURBULENCE_HPP
#define EBC_PARTITION_WALL_TURBULENCE_HPP
// ==================================================================
// EBC_Partition_wall_turbulence.hpp
//
// Element boundary condition partition for the weak imposition of
// no-slip boundary condition with wall-function applied.
// 
// Author: Xuanming Huang
// Date: Oct. 20th 2023
// ==================================================================
#include "EBC_Partition.hpp"

class EBC_Partition_wall_turbulence : public EBC_Partition
{
  public:
    // The input ElemBC should be ElemBC_3D_wall_turbulence
    EBC_Partition_wall_turbulence( const IPart * const &part,
        const Map_Node_Index * const &mnidex,
        const ElemBC * const &ebc );

    virtual ~EBC_Partition_wall_turbulence();

    // write the data to hdf5 file in folder /weak
    virtual void write_hdf5( const std::string &FileName ) const;

  protected:
    const int wall_model_type;

    // stores the local volume element id in this part
    std::vector<int> part_vol_ele_id {};

    // stores the face id of the volume element
    std::vector<int> ele_face_id {};

};

#endif
