#ifndef EBC_PARTITION_WEAK_HPP
#define EBC_PARTITION_WEAK_HPP
// ==================================================================
// EBC_Partition_weak.hpp
//
// Element boundary condition partition for the weak imposition of
// no-slip boundary condition.
// 
// Author: Xuanming Huang
// Date: Oct. 20th 2023
// ==================================================================
#include "EBC_Partition.hpp"

class EBC_Partition_weak : public EBC_Partition
{
  public:
    // The input ElemBC should be ElemBC_3D_weak
    EBC_Partition_weak( const IPart * const &part,
        const Map_Node_Index * const &mnidex,
        const ElemBC * const &ebc );

    virtual ~EBC_Partition_weak();

    // write the data to hdf5 file in folder /weak
    virtual void write_hdf5( const std::string &FileName ) const;

  protected:
    const int weak_bc_type;

    const double C_bI;

    // stores the local volume element id in this part
    std::vector< std::vector<int> > part_vol_ele_id;

    // stores the face id of the volume element
    std::vector< std::vector<int> > ele_face_id;

};

#endif