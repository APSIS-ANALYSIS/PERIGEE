#ifndef APART_NODE_HPP
#define APART_NODE_HPP
// ============================================================================
// APart_Node.hpp
// Interface for partitioned node indices, including:
// 1. re-ordered global indices of nodes;
// 2. number of local nodes, ghost nodes, etc.
// 3. dof : dofNum in the preprocessor. The total degrees of freedom
//    in the physical problem.
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
    APart_Node( const std::string &fbasename, const int &rank );

    // ------------------------------------------------------------------------
    // Destructor    
    // ------------------------------------------------------------------------
    virtual ~APart_Node();

    // ------------------------------------------------------------------------
    // This returns dofNum in the preprocessor, the total number of dof.
    // In segregated-type algorithms, this is different from dofMat, the matrix
    // problem's degrees-of-freedom, which is typically obtained by
    // ALocal_NodalBC->get_dofMat();
    // ------------------------------------------------------------------------
    virtual int get_dof() const {return dof;}

    virtual int get_nlocalnode() const {return nlocalnode;}
    
    virtual int get_nghostnode() const {return nghostnode;}
    
    virtual int get_nbadnode() const {return nbadnode;}
    
    virtual int get_nlocghonode() const {return nlocghonode;}
    
    virtual int get_ntotalnode() const {return ntotalnode;}

    // ------------------------------------------------------------------------
    // local_to_global is a mapping that maps from the local node index
    // to the whole mesh's nodal index
    // 0 <= index < nlocghonode == nlocalnode + nghostnode
    // ------------------------------------------------------------------------
    virtual int get_local_to_global(const int &index) const 
    {return local_to_global[index];}

    // ------------------------------------------------------------------------
    // node_ghost maps from [0, nghostnode) to their global mesh index
    // 0 <= index < nghostnode
    // ------------------------------------------------------------------------
    virtual int get_node_ghost(const int &index) const 
    {return node_ghost[index];}

    // ------------------------------------------------------------------------
    // node_loc maps from [0, nlocalnode) to their global mesh index
    // 0 <= index < nlocalnode
    // ------------------------------------------------------------------------
    virtual int get_node_loc(const int &index) const 
    {return node_loc[index];}

    // ------------------------------------------------------------------------
    // Determine if a global mesh node with index belongs to this subdomain 
    // ------------------------------------------------------------------------
    virtual bool is_node_local(const int &index) const
    {
      std::vector<int>::const_iterator it = find(node_loc.begin(), node_loc.end(), index);
      return ( it != node_loc.end() );
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
      SYS_T::print_fatal("Error: APart_Node::get_nlocalnode_solid is not implemented. \n");
      return -1;
    }

    virtual int get_node_loc_solid(const int &index) const
    {
      SYS_T::print_fatal("Error: APart_Node::get_node_loc_solid is not implemented. \n");
      return -1;
    }

    virtual int get_nlocalnode_fluid() const
    {
      SYS_T::print_fatal("Error: APart_Node::get_nlocalnode_fluid is not implemented. \n");
      return -1;
    }

    virtual int get_node_loc_fluid(const int &index) const
    {
      SYS_T::print_fatal("Error: APart_Node::get_node_loc_fluid is not implemented. \n");
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
    // local_to_global = node_loc appended by node_ghost
    // and their lengths are nlocghonode, nlocalnode, nghostnode, respectively
    // ------------------------------------------------------------------------
    std::vector<int> local_to_global, node_ghost, node_loc;
};

#endif
