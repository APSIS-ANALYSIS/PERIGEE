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
// local_to_global[ node_loc_solid ] gives the node index for the
// solid subdomain.
// 
// The node_loc_fluid is defined similarily for the fluid elements.
// There may be overlaps for the local solid node list and the local
// fluid node list.
//
// Author: Ju Liu
// Date: Aug. 10 2017
// ==================================================================
#include "APart_Node.hpp"
#include "ALocal_Elem_wTag.hpp"
#include "ALocal_IEN.hpp"

class APart_Node_FSI : public APart_Node
{
  public:
    APart_Node_FSI(const std::string &fileBaseName, const int &rank,
        const ALocal_Elem * const &lelem,
        const ALocal_IEN * const &lIEN );

    virtual ~APart_Node_FSI();

    virtual void print_info() const;

    // Get the solid local node number and indices
    virtual int get_nlocalnode_solid() const {return nlocalnode_solid;}

    virtual int get_node_loc_solid(const int &index) const
    {return node_loc_solid[index];}

    // Get the fluid local node number and indices
    virtual int get_nlocalnode_fluid() const {return nlocalnode_fluid;}

    virtual int get_node_loc_fluid(const int &index) const
    {return node_loc_fluid[index];}

  private:
    int nlocalnode_solid;
    std::vector<int> node_loc_solid; 

    int nlocalnode_fluid;
    std::vector<int> node_loc_fluid;
};

#endif
