#ifndef NBC_PARTITION_RING_HPP
#define NBC_PARTITION_RING_HPP
// ==================================================================
// NBC_Partition_ring.hpp
//
// Ring Nodal Boundary condition partition implementation for 3D
// meshes. 
//
// This NBC partition code is specifically designed for the class of
// NodalBC_3D_ring, which contains three additional pieces of information:
// each dirichlet node's cap id, the cap outward normal vectors, and
// their dominant component indices.
//
// Date crated: Apr. 7 2021
// Author: Ju Liu, Ingrid Lan
// ==================================================================
#include "IPart.hpp"
#include "Map_Node_Index.hpp"
#include "INodalBC.hpp"

class NBC_Partition_ring
{
  public:
    NBC_Partition_ring( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const INodalBC * const &nbc );

    virtual ~NBC_Partition_ring() = default;

    virtual void write_hdf5( const std::string &FileName ) const;

  private:
    const int cpu_rank;

    // type = 0 : all dof of ring nodes are set to be essential bc;
    // type = 1 : the dominant normal components of ring nodes are set to be essential bc;
    // type = 2 : the dominant normal & tangential components of ring nodes are set to be
    //            essential bc.
    const int ring_bc_type;

    // Number of caps (inlets, outlets)
    const int num_caps; 

    // The number of ring nodes that belong to this partition.
    // Length is Num_LD
    std::vector<int> LDN;

    // The number of Dirichlet nodes that belong to this partition.
    int Num_LD;

    // Store corresponding cap ID: [0, num_caps)
    // length Num_LD
    std::vector<int> local_cap_id;

    // Each cap's 3x3 rotation matrix for skew boundary conditions
    // length 9 x num_caps. Order of components: 11, 12, 13, 21, 22, 23, 31, 32, 33
    std::vector<double> Q;

    // Each cap's unit normal vector, length 3 x num_caps
    std::vector<double> outnormal;

    // ------------------------------------------------------------------------
    // This function is NOT allowed for ring nbc
    virtual void write_hdf5( const std::string &FileName,
        const std::string &GroupName ) const;
};

#endif
