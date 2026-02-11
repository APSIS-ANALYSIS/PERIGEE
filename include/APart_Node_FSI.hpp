#ifndef APART_NODE_FSI_HPP
#define APART_NODE_FSI_HPP
// ==================================================================
// APart_Node_FSI.hpp
//
// Partitioned node indices for FSI problem. In addition to the data
// in APart_Node class, this class has a local node list for the solid
// domain.
// 
// node_loc_solid is all the indices in LIEN[ee in solid] that also
// belongs to node_loc. The element in solid is tagged by a number, and
// the number by default should be 1. (In FSI problems, 0 domain is
// fluid domain and 1 domain is solid domain.)
//
// Usage:
// local_to_global[ node_loc_solid[ii] ] gives the nodal index for the
// solid subdomain.
// local_to_global[ node_loc_fluid[ii] ] gives the nodal index for the
// fluid subdomain.
// 
// There are non-empty intersection for the local solid node list and 
// the local fluid node list.
//
// Author: Ju Liu
// Date: Aug. 10 2017
// ==================================================================
#include "APart_Node.hpp"

class APart_Node_FSI : public APart_Node
{
  public:
    APart_Node_FSI(const std::string &fileBaseName, int rank );

    APart_Node_FSI(const HDF5_Reader * const &h5r);

    virtual ~APart_Node_FSI() = default;

    virtual void print_info() const;

    // Get the solid local node number and indices
    virtual int get_nlocalnode_solid() const 
    {return nlocalnode_solid;}

    virtual int get_node_loc_solid(int index) const
    {return node_loc_solid[index];}

    // Get the fluid local node number and indices
    virtual int get_nlocalnode_fluid() const 
    {return nlocalnode_fluid;}

    virtual int get_node_loc_fluid(int index) const
    {return node_loc_fluid[index];}

  private:
    int nlocalnode_solid, nlocalnode_fluid;
    
    std::vector<int> node_loc_solid, node_loc_fluid;
};

#endif
