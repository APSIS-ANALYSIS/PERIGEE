#ifndef NODALBC_3D_ROTATED_HPP
#define NODALBC_3D_ROTATED_HPP
// ============================================================================
// NodalBC_3D_rotated.hpp
// 
// This is an instantiation of INodalBC for 3D rotated type boundary conditions.
//
// For rotated boundary conditions, 1. there is no periodic type boundary 
// condition; 
//
// Author: Yujie Sun
// Date: Sep. 10 2024
// ============================================================================

#include "INodalBC.hpp"
#include "Tet_Tools.hpp"
#include "Hex_Tools.hpp"

class NodalBC_3D_rotated : public INodalBC
{
  public:
    // ------------------------------------------------------------------------
    // Generate the rotated bc bc given by a list of inffiles.
    // ------------------------------------------------------------------------
    NodalBC_3D_rotated( const std::string &inffile,
        const int &nFunc,
        const int &elemtype );

    virtual ~NodalBC_3D_rotated() = default;

    virtual unsigned int get_num_per_nodes() const {return 0;}

    // get the dirichlet-type nodal index on different nbc_id surfaces
    virtual unsigned int get_dir_nodes_on_rotated_surface( const unsigned int &ii ) const 
    { return dir_nodes_on_rotated_surface[ii]; }

    virtual unsigned int get_num_dir_nodes_on_rotated_surface() const
    { return num_dir_nodes_on_rotated_surface; }    

    // Access to num_node
    virtual int get_num_node() const {return num_node;}

    // Access to num_cell
    virtual int get_num_cell() const {return num_cell;}

    // Access to (surface) nLocBas
    virtual int get_nLocBas() const {return nLocBas;}

    // Access to (surface) ien
    virtual int get_ien(const int &cell, const int &lnode) const
    {return sur_ien[ nLocBas * cell + lnode ];}

    // Access to point coordinates, 
    // node = 0, ..., num_node-1.
    // dir  = 0, 1, 2.
    virtual double get_pt_xyz(const int &node, const int &dir) const
    {return pt_xyz[3*node+dir];}

    // Access to volumetric nodal index
    virtual int get_global_node(const int &node_idx) const
    {return global_node[node_idx];}

    // Access to volumetric cell index
    virtual int get_global_cell(const int &cell_idx) const
    {return global_cell[cell_idx];}

  private:
    // Disallow default constructor
    NodalBC_3D_rotated() = delete;

    // The dirichlet nodes on rotated surface
    // length = num_dir_nodes_on_rotated_surface
    std::vector<unsigned int> dir_nodes_on_rotated_surface;
    unsigned int num_dir_nodes_on_rotated_surface; 

    // number of moving surfaces and element type
    const int num_nbc, elem_type;

    // number of nodes and cells on rotated surface.
    // Note: num_node equal num_dir_nodes in this class.
    // nLocBas is 3, 6, 4 or 9.
    int num_node, num_cell, nLocBas;

    // IEN for each surface.
    // size = nLocBas x num_cell
    std::vector<int> sur_ien;

    // coordinates of all nodes on each surface
    // size = 3 x num_node
    std::vector<double> pt_xyz;

    // Surface nodes' volumetric mesh ID.
    // size = num_node
    std::vector<int> global_node;

    // Surface cell's volumetric mesh ID
    // size = num_cell
    std::vector<int> global_cell;
};

#endif