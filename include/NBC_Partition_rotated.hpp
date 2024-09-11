#ifndef NBC_PARTITION_ROTATED_HPP
#define NBC_PARTITION_ROTATED_HPP
// ==================================================================
// NBC_Partition_rotated.hpp
//
// Rotated Nodal Boundary condition partition implementation for 3D
// meshes. 
//
// This NBC partition code is specifically designed for the class of
// NodalBC_3D_rotated
// 
// The data recorded in the HDF5 file by this class will be loaded
// in the ALocal_RotatedBC class in the analysis code.
//
// Date crated: Sep. 11 2024
// Author: Yujie Sun

// ==================================================================
#include "IPart.hpp"
#include "INodalBC.hpp"
#include "Map_Node_Index.hpp"

class NBC_Partition_rotated
{
  public:
    NBC_Partition_rotated(const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const INodalBC * const &nbc );

    virtual ~NBC_Partition_rotated() = default;

    virtual void write_hdf5( const std::string &FileName ) const;

  protected:
    const int cpu_rank;
 
    // Local Dirichlet nodes that belong to rotated boundary surface
    // Length: NumLD
    std::vector<int> LDN;

    // Local Dirichlet nodes' coordinates that belong to rotated boundary surface
    std::vector<double> LDN_pt_xyz;

    // Number of Local Dirichlet nodes for rotated boundary surface
    int Num_LD;
 
    // number of local node / cell, element number of nodes
    int num_local_node, num_local_cell, cell_nLocBas;

    // local nodes' coordinates
    // Length: 3 x num_local_node
    std::vector<double> local_pt_xyz;

    // local cell IEN array
    // Length: cell_nLocBas x num_local_cell
    std::vector<int> local_cell_ien;

    // local node's global index
    // Length: num_local_node
    std::vector<int> local_global_node;

    // local node's position in the local_to_global_array
    // Length: num_local_node
    std::vector<int> local_node_pos;

    // local cell's global index
    // Length: num_local_cell
    std::vector<int> local_global_cell;
 };

#endif
