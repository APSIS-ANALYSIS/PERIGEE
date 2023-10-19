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
// should either be fully clamped (type 0), or their rotated-x dof
// (corresponding to the normal component) should be assigned as
// essential bc for in-plane motion (type 1).
//
// Author: Ju Liu, Ingrid Lan
// Date: Apr. 7 2021
// ============================================================================
#include "INodalBC.hpp"
#include "VTK_Tools.hpp"

class NodalBC_3D_ring : public INodalBC
{
  public:
    // ------------------------------------------------------------------------
    // Generate the ring bc given by the inflow and outflow files, 
    // their unit normals, and the wall file.
    // ------------------------------------------------------------------------
    NodalBC_3D_ring( const std::vector<std::string> &inflow_files,
        const std::vector< Vector_3 > &inlet_outnormal,
        const std::string &wallfile,
        const std::vector<std::string> &outflow_files,
        const std::vector< Vector_3 > &outlet_outnormal,
        const int &nFunc,
        const int &in_ring_bc_type,
        const int &elemtype );

    virtual ~NodalBC_3D_ring() = default;

    virtual unsigned int get_dir_nodes(const unsigned int &ii) const
    {return dir_nodes[ii];}

    virtual unsigned int get_per_slave_nodes(const unsigned int &ii) const
    {
      SYS_T::print_fatal("Error: periodic nodes are not defined in NodalBC_3D_ring.\n");
      return 0;
    }

    virtual unsigned int get_per_master_nodes(const unsigned int &ii) const
    {
      SYS_T::print_fatal("Error: periodic nodes are not defined in NodalBC_3D_ring.\n");
      return 0;
    }

    virtual unsigned int get_num_dir_nodes() const {return num_dir_nodes;}

    virtual unsigned int get_num_per_nodes() const {return 0;}

    virtual int get_num_caps() const { return num_caps; }

    virtual int get_ring_bc_type() const { return ring_bc_type; }

    virtual std::vector<int> get_cap_id() const { return cap_id; }

    virtual std::vector<double> get_outnormal() const { return outnormal; } 

    virtual std::vector<double> get_rotation_matrix() const { return Q; } 

  private:
    std::vector<unsigned int> dir_nodes;
    unsigned int num_dir_nodes;

    NodalBC_3D_ring() = delete; 

    // Ring node bc type
    // type = 0 : all dof of ring nodes are assigned as essential bc (clamped)
    // type = 1 : rotated-x dof of ring nodes are assigned as essential bc (in-plane motion) 
    const int ring_bc_type;

    // Number of caps (inlets, outlets)
    int num_caps;

    // Store corresponding cap ID: [0, num_caps)
    // Inlet cap surface has ID 0, followed by the outlet caps
    // length num_dir_nodes
    std::vector<int> cap_id;

    // Each cap's 3x3 global-to-local transformation matrix for skew boundary conditions
    // 9 components per cap: 11, 12, 13, 21, 22, 23, 31, 32, 33
    // Transforms from Cartesian x-y-z frame to normal-radial-tangential frame
    std::vector<double> Q;

    // Each cap's unit normal vector 
    // length 3 x num_caps
    std::vector<double> outnormal;

    // Compute centroid coordinates given a cap's nodal coordinates
    Vector_3 compute_cap_centroid( const std::vector<double> &pts ) const;
};

#endif
