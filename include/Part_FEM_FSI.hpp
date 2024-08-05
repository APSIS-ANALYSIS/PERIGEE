#ifndef PART_FEM_FSI_HPP
#define PART_FEM_FSI_HPP
// ============================================================================
// Part_FEM_FSI.hpp
//
// Object: Partition 3D tet mesh into subdomains. 
//         The mesh is partitioned into sub-domains that are tagged with a 
//         physical tag. In the partitioned file, there will be a physical tag,
//         which is an integer array with the length equaling to the number of 
//         local elements.
//
// Date Created: July 27 2017
// ============================================================================
#include "Part_FEM.hpp"

class Part_FEM_FSI : public Part_FEM
{
  public:
    Part_FEM_FSI( const IMesh * const &mesh,
        const IGlobal_Part * const &gpart,
        const Map_Node_Index * const &mnindex,
        const IIEN * const &IEN,
        const std::vector<double> &ctrlPts,
        const std::vector<int> &phytag,
        const std::vector<int> &node_f,
        const std::vector<int> &node_s,
        const int &in_cpu_rank, 
        const int &in_cpu_size,
        const int &in_elemType,
        const int &in_start_idx,
        const Field_Property &in_fp );

    virtual ~Part_FEM_FSI() = default;

    virtual void write( const std::string &inputFileName ) const;

  protected:
    std::vector<int> elem_phy_tag {};
    
    // The location in the node_loc for the local fluid/solid nodes
    // e.g. node_loc_fluid.size() gives the number of nodes belonging
    // to a fluid domain in this partition, and node_loc_fluid[ii] gives
    // the location in node_loc for that node's index.
    // Storing the index in node_loc array makes things more flexible,
    // because sometimes, we manage arrays using local partition's numbering
    // sometimes, we directly set entries using the global numbering.
    std::vector<int> node_loc_fluid {};
    std::vector<int> node_loc_solid {}; 

    int nlocalnode_fluid, nlocalnode_solid;

    // DOF mapper
    const int start_idx;
};

#endif
