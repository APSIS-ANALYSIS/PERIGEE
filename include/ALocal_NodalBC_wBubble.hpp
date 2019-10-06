#ifndef ALOCAL_NODALBC_WBUBBLE_HPP
#define ALOCAL_NODALBC_WBUBBLE_HPP
// ==================================================================
// ALocal_NodalBC_wBubble.hpp
// 
// This re-generate the local ID array based on the bubble enrichment.
//
// The LID, LDN, LPSN, LPMN, LocalMaster and LocalMasterSlave
// will be updated by the geo2phy mapping. These vectors store the
// geometry nodal indices. With bubble enrichment, their physical
// indices will be re-organized.
//
// Since the bubble functions always have a trace zero over the 
// boundaries, we only need to append the LID array with the bubble
// nodal indices.
//
// ntotnode = pNode -> get_ntotalnode();
//
// Author: Ju Liu
// Date: Nov. 27 2017
// ==================================================================
#include "ALocal_NodalBC.hpp"
#include "APart_Node.hpp"

class ALocal_NodalBC_wBubble : public ALocal_NodalBC
{
  public:
    ALocal_NodalBC_wBubble( const std::string &fileBaseName,
        const int &cpu_rank, const APart_Node * const &pNode );

    virtual ~ALocal_NodalBC_wBubble();

    virtual void print_info() const;

    virtual int get_LID(const int &dof_index, const int &node) const
    {return LID[dof_index * ntotnode + node];}

  private:
    const int ntotnode;
};

#endif
