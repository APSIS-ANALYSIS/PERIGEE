#ifndef APART_NODE_HPP
#define APART_NODE_HPP
// ============================================================================
// APart_Node.hpp
// Class storing information of partitioned node indices, including:
// 1. The re-ordered global indices of {local, ghost, local+ghost} nodes;
// 2. number of local, ghost, local+ghost nodes, etc.
// 3. dof : dofNum in the preprocessor, i.e., the degrees of freedom attached to
//    the nodes.
// 
// In the local_to_global array, there are 
//             nlocghonode = nlocalnode + nghostnode
// entries. The first nlocalnode entries store the node_loc array;
// the following nghostnode entries store the node_ghost array.
//
// Author: Ju Liu
// Date: Nov. 26th 2013
// ============================================================================
#include "HDF5_Reader.hpp"

class APart_Node
{
  public:
    // ------------------------------------------------------------------------
    // Constructor: Read the HDF5 file from mesh partitioning.
    //              fbasename : the base name of the h5 file
    //              rank      : the cpu rank
    // ------------------------------------------------------------------------
    APart_Node( const std::string &fbasename, int rank );

    APart_Node( const HDF5_Reader * const &h5r );

    // ------------------------------------------------------------------------
    // Destructor    
    // ------------------------------------------------------------------------
    virtual ~APart_Node() = default;

    // ------------------------------------------------------------------------
    // This returns dofNum in the preprocessor, the total number of dof.
    // In segregated-type algorithms, this is different from the matrix
    // problem's degrees-of-freedom, which is typically obtained by
    // ALocal_NBC->get_dof_LID();
    // ------------------------------------------------------------------------
    virtual int get_dof() const {return dof;}

    virtual int get_nlocalnode() const {return nlocalnode;}
    
    virtual int get_nghostnode() const {return nghostnode;}
    
    virtual int get_nbadnode() const {return nbadnode;}
    
    virtual int get_nlocghonode() const {return nlocghonode;}
    
    virtual int get_ntotalnode() const {return ntotalnode;}

    // ------------------------------------------------------------------------
    // local_to_global is a mapping that maps from the local node index
    // to the global/volumetric mesh's nodal index
    // Input: local nodal index with ranges 
    //        0 <= ii < nlocghonode == nlocalnode + nghostnode
    // ------------------------------------------------------------------------
    virtual int get_local_to_global(int ii) const 
    {return local_to_global[ii];}

    // ------------------------------------------------------------------------
    // node_ghost maps from [0, nghostnode) to their global/volume mesh index
    // 0 <= ii < nghostnode
    // ------------------------------------------------------------------------
    virtual int get_node_ghost(int ii) const {return node_ghost[ii];}

    // ------------------------------------------------------------------------
    // node_loc maps from [0, nlocalnode) to their global/volume mesh index
    // 0 <= index < nlocalnode
    // ------------------------------------------------------------------------
    virtual int get_node_loc(int ii) const {return node_loc[ii];}

    // ------------------------------------------------------------------------
    // Determine if a global mesh node with index belongs to this subdomain.
    // Input: ii is a global/volumetric mesh nodal index. 
    // ------------------------------------------------------------------------
    virtual bool is_node_local(int ii) const
    {
      return ( find(node_loc.begin(), node_loc.end(), ii) != node_loc.end() );
    }

    // ------------------------------------------------------------------------
    // Return this subdomain rank
    // ------------------------------------------------------------------------
    virtual int get_rank() const {return cpu_rank;}

    virtual void print_info() const;

    // ------------------------------------------------------------------------
    // Virtual functions for the _FSI derived class
    // ------------------------------------------------------------------------
    virtual int get_nlocalnode_solid() const
    {
      SYS_T::print_fatal("Error: APart_Node::get_nlocalnode_solid is not implemented.\n");
      return -1;
    }

    virtual int get_node_loc_solid(int index) const
    {
      SYS_T::print_fatal("Error: APart_Node::get_node_loc_solid is not implemented.\n");
      return -1;
    }

    virtual int get_nlocalnode_fluid() const
    {
      SYS_T::print_fatal("Error: APart_Node::get_nlocalnode_fluid is not implemented.\n");
      return -1;
    }

    virtual int get_node_loc_fluid(int index) const
    {
      SYS_T::print_fatal("Error: APart_Node::get_node_loc_fluid is not implemented.\n");
      return -1;
    }
    
    // ------------------------------------------------------------------------
    // Virtual functions for the _Rotated derived class
    // ------------------------------------------------------------------------
    virtual int get_nlocalnode_rotated() const
    {
      SYS_T::print_fatal("Error: APart_Node::get_nlocalnode_rotated is not implemented.\n");
      return -1;
    }

    virtual int get_node_loc_rotated(int index) const
    {
      SYS_T::print_fatal("Error: APart_Node::get_node_loc_rotated is not implemented.\n");
      return -1;
    }

    virtual int get_nlocalnode_fixed() const
    {
      SYS_T::print_fatal("Error: APart_Node::get_nlocalnode_fixed is not implemented.\n");
      return -1;
    }

    virtual int get_node_loc_fixed(int index) const
    {
      SYS_T::print_fatal("Error: APart_Node::get_node_loc_fixed is not implemented.\n");
      return -1;
    }

  protected:
    // ------------------------------------------------------------------------
    // rank of the CPU that identifies the subdomain
    // ------------------------------------------------------------------------
    const int cpu_rank;

    // ------------------------------------------------------------------------
    // nlocghonode = nlocalnode + nghostnode
    // ntotalnode = nlocghonnode + nbadnode
    // dof here reads the dofNum from the Global_Mesh_Info group.
    // ------------------------------------------------------------------------
    int nlocalnode, nghostnode, nbadnode, nlocghonode, ntotalnode, dof;
    
    // ------------------------------------------------------------------------
    // local_to_global = [ node_loc ] appended by [ node_ghost ].
    // The three vectors have lengths nlocghonode, nlocalnode, nghostnode, resp.
    // ------------------------------------------------------------------------
    std::vector<int> local_to_global {}, node_ghost {}, node_loc {};

    APart_Node() = delete;
};

#endif
