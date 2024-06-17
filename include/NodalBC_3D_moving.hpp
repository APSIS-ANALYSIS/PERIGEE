#ifndef NODALBC_3D_MOVING_HPP
#define NODALBC_3D_MOVING_HPP
// ============================================================================
// NodalBC_3D_moving.hpp
// 
// This is an instantiation of INodalBC for 3D moving type boundary conditions.
//
// For moving boundary conditions, 1. there is no periodic type boundary 
// condition; 
//
// Author: Yujie Sun
// Date: Jun. 13 2024
// ============================================================================

#include "INodalBC.hpp"
#include "Tet_Tools.hpp"
#include "Hex_Tools.hpp"

class NodalBC_3D_moving : public INodalBC
{
  public:
    // ------------------------------------------------------------------------
    // Generate the moving bc bc given by a list of inffiles.
    // ------------------------------------------------------------------------
    NodalBC_3D_moving( const std::vector<std::string> &inffileList,
        const int &nFunc,
        const int &elemtype );

    virtual ~NodalBC_3D_moving() = default;

    virtual unsigned int get_dir_nodes(const unsigned int &ii) const
    {return dir_nodes[ii];}

    virtual unsigned int get_num_dir_nodes() const {return num_dir_nodes;}   

  	virtual unsigned int get_per_slave_nodes(const unsigned int &ii) const
    {
      SYS_T::print_fatal("Error: NodalBC_3D_moving::get_per_slave_nodes: periodic nodes are not defined.\n");
      return 0;
    }

    virtual unsigned int get_per_master_nodes(const unsigned int &ii) const
    {
      SYS_T::print_fatal("Error: NodalBC_3D_moving::get_per_master_nodes: periodic nodes are not defined.\n");
      return 0;
    }

    virtual unsigned int get_num_per_nodes() const {return 0;}

    // get the dirichlet-type nodal index on different nbc_id surfaces
    virtual unsigned int get_dir_nodes_on_surface( const int &nbc_id, 
        const unsigned int &ii ) const { return dir_nodes_on_surface[nbc_id][ii]; }

    virtual unsigned int get_num_dir_nodes_on_surface( const int &nbc_id ) const
    { return num_dir_nodes_on_surface[nbc_id]; }    

    // Access to num_nbc
    virtual int get_num_nbc() const {return num_nbc;}

    // Access to num_node
    virtual int get_num_node(const int &nbc_id) const {return num_node[nbc_id];}

    // Access to num_cell
    virtual int get_num_cell(const int &nbc_id) const {return num_cell[nbc_id];}

    // Access to (surface) nLocBas
    virtual int get_nLocBas(const int &nbc_id) const {return nLocBas[nbc_id];}

    // Access to (surface) ien
    virtual int get_ien(const int &nbc_id, const int &cell, const int &lnode) const
    {return sur_ien[nbc_id][ nLocBas[nbc_id] * cell + lnode ];}

    // Access to point coordinates, 
    // node = 0, ..., num_node-1.
    // dir  = 0, 1, 2.
    virtual double get_pt_xyz(const int &nbc_id, const int &node, const int &dir) const
    {return pt_xyz[nbc_id][3*node+dir];}

    // Access to volumetric nodal index
    virtual int get_global_node(const int &nbc_id, const int &node_idx) const
    {return global_node[nbc_id][node_idx];}

    // Access to volumetric cell index
    virtual int get_global_cell(const int &nbc_id, const int &cell_idx) const
    {return global_cell[nbc_id][cell_idx];}

  private:
    // Disallow default constructor
    NodalBC_3D_moving() = delete;

    // The dirichlet nodes on each moving surface
    // length num_nbc x num_dir_nodes_on_surface[ii]
    std::vector< std::vector<unsigned int> > dir_nodes_on_surface;
    std::vector<unsigned int> num_dir_nodes_on_surface;

    // All dirichlet nodes stored in a row and the total number of dirichlet
    // nodes
    std::vector<unsigned int> dir_nodes;
    unsigned int num_dir_nodes;   

    // number of moving surfaces and element type
    const int num_nbc, elem_type;

    // number of nodes and cells on each surface. Length num_nbc.
    // Note: num_node[ii] equal num_dir_nodes[ii] in this class.
    // nLocBas[ii] is 3, 6, 4 or 9.
     std::vector<int> num_node, num_cell, nLocBas;

    // IEN for each surface.
    // num_nbc times ( nLocBas[ii] x num_cell[ii] ) in size.
    std::vector< std::vector<int> > sur_ien;

    // coordinates of all nodes on each surface
    // num_nbc times ( 3 x num_node[ii] ) in size.
    std::vector< std::vector<double> > pt_xyz;

    // Surface nodes' volumetric mesh ID.
    // num_nbc times num_node[ii] in size.
    std::vector< std::vector<int> > global_node;

    // Surface cell's volumetric mesh ID
    // num_nbc times num_cell[ii] in size.
    std::vector< std::vector<int> > global_cell;
};

#endif
