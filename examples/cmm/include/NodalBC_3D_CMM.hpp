#ifndef NODALBC_3D_CMM_HPP
#define NODALBC_3D_CMM_HPP
// ============================================================================
// NodalBC_3D_CMM.hpp
//
// This is an instantiation of INodalbc for 3D problems by reading the
// dirichlet nodes from the vtp file.
//
// This class is designed to handle the mesh for unstructural tetrahedral
// mesh generated by automatic mesher like tetgen.
// 
// The data contained in this class include:
// dir_nodes : the nodal indices for the Dirichlet nodes;
// num_dir_nodes : the number of the Dirichlet nodes, i.e., the length
//                 of the dir_nodes array;
// per_slave_nodes : the nodal indices for the slave nodes;
// per_master_nodes : the nodal indices for the master nodes;
// num_per_nodes : the number of periodic-type nodes, i.e., the length 
//                 of the per_slave_nodes / per_master_nodes.
//
// ID : the vector for the ID array, which is generated based on the 
//      nFunc, the number of total basis functions, and the dir_nodes
//      and per_slave_nodes/per_master_nodes. Once the dir_nodes, per_
//      xxx_nodes, and nFunc are given, the ID array will be generated
//      by the function create_ID.
//
// Date: Jan. 6 2017
// Author: Ju Liu
// ============================================================================
#include "INodalBC.hpp"
#include "Tet_Tools.hpp"

class NodalBC_3D_CMM : public INodalBC
{
  public:
    // ------------------------------------------------------------------------ 
    // Default constructor: clear the dir_nodes, per_slave_nodes,
    // per_master_nodes; set num_dir_nodes, num_per_nodes to be zero;
    // set ID based on the above "no-nodal bc" setting.
    // ------------------------------------------------------------------------ 
    NodalBC_3D_CMM( const int &nFunc, const bool &is_all_node = false );
 
    NodalBC_3D_CMM( const INodalBC * const &nbc_inflow,
        const INodalBC * const &nbc_ring,
        const INodalBC * const &nbc_wall,
        const int &comp, const int &nFunc,
        const int &cmm_bc_type );  

    virtual ~NodalBC_3D_CMM() {};

    virtual unsigned int get_dir_nodes(const unsigned int &ii) const
    {return dir_nodes[ii];}

    virtual unsigned int get_per_slave_nodes(const unsigned int &ii) const
    {
      SYS_T::print_fatal("Error: periodic nodes are not defined in NodalBC_3D_CMM.\n");
      return 0;
    }

    virtual unsigned int get_per_master_nodes(const unsigned int &ii) const
    {
      SYS_T::print_fatal("Error: periodic nodes are not defined in NodalBC_3D_CMM.\n");
      return 0;
    }

    virtual unsigned int get_num_dir_nodes() const {return num_dir_nodes;}

    virtual unsigned int get_num_per_nodes() const {return 0;}

  private:
    std::vector<unsigned int> dir_nodes;
    unsigned int num_dir_nodes;

    NodalBC_3D_CMM() {};
};

#endif
