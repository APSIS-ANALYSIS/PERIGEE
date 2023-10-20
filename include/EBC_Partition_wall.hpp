#ifndef EBC_PARTITION_WALL_HPP
#define EBC_PARTITION_WALL_HPP
// ==================================================================
// EBC_Partition_wall.hpp
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
#include "EBC_Partition.hpp"

class EBC_Partition_wall : public EBC_Partition
{
  public:
    // The input ElemBC should be ElemBC_3D_wall
    EBC_Partition_wall( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const ElemBC * const &ebc );

    virtual ~EBC_Partition_wall();

    // write the data to hdf5 file in folder /ebc_wall 
    virtual void write_hdf5( const std::string &FileName ) const;

  protected:
    // Length is num_local_node[0] 
    std::vector<double> part_thickness;
    std::vector<double> part_youngsmod;
    std::vector<double> part_springconst;
    std::vector<double> part_dampingconst;

    // The number of surface nodes beloging to this subdomain
    int num_local_node_on_sur;

    // The position of all local nodes on the wall surface in the 
    // local_to_global array ( of the volumetric nodal indices).
    // Here, partition nodes are all surface nodes belonging to the CPU's subdomain,
    // which may not be associated with any cell in the local partition.
    // It is used for nodal update of surface solutions only.
    // Length is  num_local_node_on_sur
    std::vector<int> local_node_on_sur_pos;
    
    // This function is NOT allowed for wall ebc. 
    virtual void write_hdf5( const std::string &FileName, 
        const std::string &GroupName ) const;
};

#endif
