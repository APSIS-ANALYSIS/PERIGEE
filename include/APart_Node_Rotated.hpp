#ifndef APART_NODE_ROTATED_HPP
#define APART_NODE_ROTATED_HPP
// ==================================================================
// APart_Node_Rotated.hpp
//
// Partitioned node indices for ALE_ROTATED problem. In addition to the data
// in APart_Node class, this class has a local node list for the fixed/rotated
// domain.
// 
// node_loc_rotated is all the indices that also belongs to node_loc. 
// The element in rotated domain is tagged by a number, and the number 
// by default should be 1. (In ALE_ROTATED problems, 0 domain is
// fixed domain and 1 domain is rotated domain.)
//
// Usage:
// local_to_global[ node_loc_rotated[ii] ] gives the nodal index for the
// rotated domain.
// local_to_global[ node_loc_fixed[ii] ] gives the nodal index for the
// fixed domain.
// 
// There are empty intersection for the local rotated node list and 
// the local fixed node list.
//
// Author: Yujie Sun
// Date: Aug. 14 2024
// ==================================================================
#include "APart_Node.hpp"

class APart_Node_Rotated : public APart_Node
{
  public:
    APart_Node_Rotated(const std::string &fileBaseName, const int &rank );

    virtual ~APart_Node_Rotated() = default;

    virtual void print_info() const;

    // Get the solid local node number and indices
    virtual int get_nlocalnode_rotated() const 
    {return nlocalnode_rotated;}

    virtual int get_node_loc_rotated(const int &index) const
    {return node_loc_rotated[index];}

    // Get the fluid local node number and indices
    virtual int get_nlocalnode_fixed() const 
    {return nlocalnode_fixed;}

    virtual int get_node_loc_fixed(const int &index) const
    {return node_loc_fixed[index];}

  private:
    int nlocalnode_rotated, nlocalnode_fixed;
    
    std::vector<int> node_loc_rotated, node_loc_fixed;
};

#endif
