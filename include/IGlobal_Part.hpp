#ifndef IGLOBAL_PART_HPP
#define IGLOBAL_PART_HPP
// ==================================================================
// IGlobal_Part.hpp
// Object:
// Get the global partition information. 
//
// example:
// Use METIS to get the global epart and npart. Then these are used 
// to create the partition for each processor. This helps us avoid 
// rerunning the same METIS call.
//
// Author: Ju Liu
// Date: Oct. 2 2013
// ==================================================================
#include "metis.h"

class IGlobal_Part
{
  public:
    IGlobal_Part() = default;
    
    virtual ~IGlobal_Part() = default;

    virtual idx_t get_epart( const int &ee ) const = 0;
    
    virtual idx_t get_npart( const int &nn, const int &field = 0 ) const = 0;

    virtual bool get_isMETIS() const = 0;
    
    virtual bool get_isDual() const = 0;
    
    virtual int get_dual_edge_ncommon() const = 0;

    virtual bool is_serial() const = 0;
};

#endif
