#ifndef ALOCAL_IEN_WBUBBLE_HPP
#define ALOCAL_IEN_WBUBBLE_HPP
// ==================================================================
// ALocal_IEN_wBubble.hpp
//
// Local IEN for analysis with bubble enrichment.
// 
// This class reads in the original (before bubble enrichment) local
// IEN from HDF5 files.
//
// Then do the following update: 
// Copy LIEN to a temp vector
// 
// Read in the value of nlocghonode from Local_Node group in the part
// file.
//
// Redefine nLocBas += nbubble_per_cell
// 
// for original nodes, 0 <= ii < nLocBas_ori
// LIEN[ee*nLocBas + ii] = temp[ee*nLocBas_ori + ii];
// 
// for bubble nodes, nLocBas_ori <= ii < nLocBas
// LIEN[ee*nLocBas + ii] = nlocghonode + ee * nbubble_per_cell + ii
// 
// Author: Ju Liu 
// Date: Nov. 27 2017
// ==================================================================
#include "ALocal_IEN.hpp"

class ALocal_IEN_wBubble : public ALocal_IEN
{
  public:
    ALocal_IEN_wBubble( const std::string &fileBaseName, 
        const int &cpu_rank, const int &nbubble_per_cell );


    virtual ~ALocal_IEN_wBubble();


    virtual void print_info() const;
};

#endif
