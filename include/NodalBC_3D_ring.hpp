#ifndef NODALBC_3D_RING_HPP
#define NODALBC_3D_RING_HPP
// ============================================================================
// NodalBC_3D_ring.hpp
// 
// This is an instantiation of INodalBC for 3D ring boundary conditions
// enabling in-plane motion for the inlet and outlet ring/boundary
// nodes in CMM.
//
// For ring boundary conditions, 1. there is no periodic type
// boundary condition; 2. For a given inlet/outlet, the ring nodes
// should only be constrained in the dominant component of the
// unit normal.
//
// Author: Ju Liu, Ingrid Lan
// Date: Apr. 7 2021
// ============================================================================
#include "INodalBC.hpp"
#include "Tet_Tools.hpp"
#include "Vector_3.hpp"

class NodalBC_3D_ring : public INodalBC
{
  public:
    // ------------------------------------------------------------------------
    // Generate an empty ring type boundary condition class
    // ------------------------------------------------------------------------
    NodalBC_3D_ring(const int &nFunc);
    
    // ------------------------------------------------------------------------
    // Generate the ring bc given by the inflow and outflow files, 
    // their unit normals, and the wall file.
    // ------------------------------------------------------------------------
    NodalBC_3D_ring( const std::string &inflow_file,
        const Vector_3 &inlet_outnormal,
        const std::string &wallfile,
        const std::vector<std::string> &outflow_files,
        const std::vector< Vector_3 > &outlet_outnormal,
        const int &nFunc,
        const int &in_ring_bc_type,
        const int &elemtype = 501 );

    virtual ~NodalBC_3D_ring() {};

    virtual int get_num_caps() const { return num_caps; }

    virtual int get_ring_bc_type() const { return ring_bc_type; }

    virtual std::vector<int> get_cap_id() const { return cap_id; }

    virtual std::vector<int> get_dominant_n_comp() const { return dominant_n_comp; }

    virtual std::vector<int> get_dominant_t_comp() const { return dominant_t_comp; }

    virtual std::vector<double> get_outnormal() const { return outnormal; } 

    virtual std::vector<double> get_rotation_matrix() const { return Q; } 

  private:
    NodalBC_3D_ring() : ring_bc_type(0) {};

    // Ring node bc type
    // type = 0 : all dof of ring nodes are set to be essential bc;
    // type = 1 : the dominant normal components of ring nodes are set to be essential bc;
    // type = 2 : the dominant normal & tangential components of ring nodes are set to be
    //            essential bc.
    // type = 3 : all dof of inlet ring nodes are set to be essential bc. Only the dominant
    //            normal components of outlet ring nodes are set to be essential bc.
    // type = 4 : all dof of inlet ring nodes are set to be essential bc. Only the dominant
    //            normal & tangential components of outlet ring nodes are set to be essential bc.
    // type = 5 : all dof of a single ring node on each cap are set as essential bc. For all
    //            remaining ring nodes, only the dominant normal components are set as essential bc.
    const int ring_bc_type;

    // Number of caps (inlets, outlets)
    int num_caps;

    // Store corresponding cap ID: [0, num_caps)
    // Inlet cap surface has ID 0, followed by the outlet caps
    // length num_dir_nodes
    std::vector<int> cap_id;

    // Each cap's 3x3 global-to-local transformation matrix for skew boundary conditions
    // 9 components per cap: 11, 12, 13, 21, 22, 23, 31, 32, 33
    std::vector<double> Q;

    // Dominant component index of each cap's unit normal vector: 0, 1, or 2
    // length num_caps
    std::vector<int> dominant_n_comp;

    // Dominant component index of each ring node's unit tangential vector: 0, 1, or 2
    // length num_dir_nodes
    std::vector<int> dominant_t_comp;

    // Each cap's unit normal vector 
    // length 3 x num_caps
    std::vector<double> outnormal;

    // Each ring node's unit tangential vector 
    // length 3 x num_dir_nodes
    std::vector<double> tangential;

    // Compute centroid coordinates given a cap's nodal coordinates
    void compute_cap_centroid( const std::vector<double> &pts, Vector_3 &centroid ) const;
};

#endif
