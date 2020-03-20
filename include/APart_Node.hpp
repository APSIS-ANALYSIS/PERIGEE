#ifndef APART_NODE_HPP
#define APART_NODE_HPP
// ==================================================================
// APart_Node.hpp
// Interface for partitioned node indices, including:
// 1. re-ordered global indeices of nodes;
// 2. number of local nodes, ghost nodes, etc.
// 3. dof : the dof from the dofNum in the preprocessor. The total
//    degrees of freedom in the physical problem.
// 
// In local_to_global array, there are 
//             nlocghonode = nlocalnode + nghostnode
// entries. The first nlocalnode entries store the node_loc array;
// the following nghostnode entries store the node_ghost array.
//
// Author: Ju Liu
// Date: Nov. 26th 2013
// ==================================================================
#include "HDF5_Reader.hpp"

class APart_Node
{
  public:
    // --------------------------------------------------------------
    // Constructor: Read HDF5 file
    //              fbasename : the base name of the h5 file
    //              rank : the cpu rank
    // --------------------------------------------------------------
    APart_Node( const std::string &fbasename, const int &rank );
    
    virtual ~APart_Node();

    // This returns dofNum in the preprocessor, the total number of dof.
    // In segregated-type algorithms, this one is different from dofMat,
    // or e.g. ALocal_NodalBC->get_dofMat();
    virtual int get_dof() const {return dof;}

    virtual int get_nlocalnode() const {return nlocalnode;}
    
    virtual int get_nghostnode() const {return nghostnode;}
    
    virtual int get_nbadnode() const {return nbadnode;}
    
    virtual int get_nlocghonode() const {return nlocghonode;}
    
    virtual int get_ntotalnode() const {return ntotalnode;}

    virtual int get_local_to_global(const int &index) const 
    {return local_to_global[index];}

    virtual int get_node_ghost(const int &index) const 
    {return node_ghost[index];}

    virtual int get_node_loc(const int &index) const 
    {return node_loc[index];}

    virtual bool is_node_local(const int &index) const
    {
      std::vector<int>::const_iterator it = find(node_loc.begin(), node_loc.end(), index);
      return ( it != node_loc.end() );
    }

    virtual int get_rank() const {return cpu_rank;}

    virtual void print_info() const;

    // Virtual functions for the _FSI derived class
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
    int cpu_rank, dof;
    int nlocalnode, nghostnode, nbadnode, nlocghonode, ntotalnode;
    std::vector<int> local_to_global, node_ghost, node_loc;
};

#endif
