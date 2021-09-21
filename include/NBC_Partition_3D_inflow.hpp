#ifndef NBC_PARTITION_3D_INFLOW_HPP
#define NBC_PARTITION_3D_INFLOW_HPP
// ==================================================================
// NBC_Partition_3D_inflow.hpp
//
// Inflow Nodal Boundary condition partition implementation for 3D
// meshes. 
//
// This NBC partition code is specifically designed for the class of
// NodalBC_3D_inflow, which contains two additional information:
// inflow surface area and the inflow surface outward normal vector.
// 
// The data recorded in the HDF5 file by this class will be loaded
// in the ALocal_Inflow_NodalBC class in the analysis code.
//
// Date crated: Aug. 9 2017
// Author: Ju Liu
// ==================================================================
#include "NBC_Partition_3D.hpp"

class NBC_Partition_3D_inflow : public NBC_Partition_3D
{
  public:
    NBC_Partition_3D_inflow(const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const INodalBC * const &nbc );

    virtual ~NBC_Partition_3D_inflow();

    virtual void write_hdf5( const char * FileName ) const;

  private:
    double actarea, facearea;

    // unit outward normal vector on this surface
    std::vector<double> outvec;

    // number of boundary points of this surface
    int num_out_bc_pts;

    // spatial coordinates of the centroid
    Vector_3 centroid;

    // coordindates of the boundary points
    std::vector<double> outline_pts;

    // number of local node / cell, element number of nodes
    int num_local_node, num_local_cell, cell_nLocBas;

    // local nodes' coordinates, length 3 x num_local_node
    std::vector<double> local_pt_xyz;

    // local cell IEN array, size cell_nLocBas x num_local_cell
    std::vector<int> local_tri_ien;

    // local node's global index, length num_local_node
    std::vector<int> local_global_node;

    // local node's postion in the local_to_global_array
    std::vector<int> local_node_pos;

    // local cell's global index, length num_local_cell
    std::vector<int> local_global_cell;
};

#endif
