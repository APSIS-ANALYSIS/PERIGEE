#ifndef NBC_PARTITION_MOVING_HPP
#define NBC_PARTITION_MOVING_HPP
// ==================================================================
// NBC_Partition_moving.hpp
//
// mvoing Nodal Boundary condition partition implementation for 3D
// meshes. 
//
// This NBC partition code is specifically designed for the class of
// NodalBC_3D_moving
// 
// The data recorded in the HDF5 file by this class will be loaded
// in the ALocal_MovingBC class in the analysis code.
//
// Date crated: Jun. 13 2024
// Author: Ju Liu
// ==================================================================
#include "IPart.hpp"
#include "INodalBC.hpp"
#include "Map_Node_Index.hpp"

class NBC_Partition_moving
{
  public:
    NBC_Partition_moving(const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const INodalBC * const &nbc );

    virtual ~NBC_Partition_moving() = default;

    virtual void write_hdf5( const std::string &FileName ) const;

  protected:
    const int cpu_rank, num_nbc;
 
    // Local Dirichlet nodes that belong to different inlet surfaces
    // Length: num_nbc x NumLD[ii], for 0 <= ii < num_nbc
    std::vector< std::vector<int> > LDN;

    // Number of Local Dirichlet nodes for each inlet surface
    // Length: num_nbc
    std::vector<int> Num_LD;
 
    // number of local node / cell, element number of nodes. Length num_nbc
    std::vector<int> num_local_node, num_local_cell, cell_nLocBas;

    // local nodes' coordinates
    // num_nbc times (3 x num_local_node[ii])
    std::vector< std::vector<double> > local_pt_xyz;

    // local cell IEN array
    // num_nbc times (cell_nLocBas[ii] x num_local_cell[ii])
    std::vector< std::vector<int> > local_cell_ien;

    // local node's global index
    // num_nbc x num_local_node[ii]
    std::vector< std::vector<int> > local_global_node;

    // local node's position in the local_to_global_array
    // num_nbc x num_local_node[ii]
    std::vector< std::vector<int> > local_node_pos;

    // local cell's global index
    // num_nbc x num_local_cell[ii]
    std::vector< std::vector<int> > local_global_cell;
 };

#endif
