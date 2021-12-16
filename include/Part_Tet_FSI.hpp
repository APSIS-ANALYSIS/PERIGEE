#ifndef PART_TET_FSI_HPP
#define PART_TET_FSI_HPP
// ==================================================================
// Part_Tet_FSI.hpp
//
// Object: Partition 3D tet mesh into subdomains. The
//         mesh contains sub-domains that has a physical tag.
//         In the partitioned file, there will be a physical tag, an
//         integer array with the length equaling to the number of 
//         local elements.
//
// Date Created: July 27 2017
// ==================================================================
#include "Part_Tet.hpp"

class Part_Tet_FSI : public Part_Tet
{
  public:
    Part_Tet_FSI( const IMesh * const &mesh,
        const IGlobal_Part * const &gpart,
        const Map_Node_Index * const &mnindex,
        const IIEN * const &IEN,
        const std::vector<double> &ctrlPts,
        const std::vector<int> &phytag,
        const std::vector<int> &node_f,
        const std::vector<int> &node_s,
        const int &field,
        const int &in_start_idx,
        const int &in_cpu_rank, 
        const int &in_cpu_size,
        const int &in_elemType,
        const bool &in_is_geo_field );

    virtual ~Part_Tet_FSI();

    virtual void write( const char * inputFileName ) const;

  protected:
    std::vector<int> elem_phy_tag;
    
    // The location in the node_loc for the local fluid/solid nodes
    // e.g. node_loc_fluid.size() gives the number of nodes belonging
    // to a fluid domain in this partition, and node_loc_fluid[ii] gives
    // the location in node_loc for that node's index.
    // Storing the index in node_loc array makes things more flexible,
    // because sometimes, we manage arrays using local partition's numbering
    // sometimes, we directly set entries using the global numbering.
    std::vector<int> node_loc_fluid, node_loc_solid; 

    int nlocalnode_fluid, nlocalnode_solid;

    // DOF mapper
    const int start_idx;

    // Flag that determines if the field is a geometry field
    const bool is_geo_field;
};

#endif
