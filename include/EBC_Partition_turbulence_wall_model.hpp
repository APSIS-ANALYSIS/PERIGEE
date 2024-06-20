#ifndef EBC_PARTITION_TURBULENCE_WALL_MODEL_HPP
#define EBC_PARTITION_TURBULENCE_WALL_MODEL_HPP
// ==================================================================
// EBC_Partition_turbulence_wall_model.hpp
//
// Element boundary condition partition for the weak imposition of
// no-slip boundary condition with wall-function applied.
// 
// Author: Xuanming Huang
// Date: Oct. 20th 2023
// ==================================================================
#include "EBC_Partition.hpp"

class EBC_Partition_turbulence_wall_model : public EBC_Partition
{
  public:
    // The input ElemBC should be ElemBC_3D_turbulence_wall_model
    EBC_Partition_turbulence_wall_model( const IPart * const &part,
        const Map_Node_Index * const &mnidex,
        const ElemBC * const &ebc );

    virtual ~EBC_Partition_turbulence_wall_model();

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
