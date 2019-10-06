#ifndef I_GLOBAL_PART_HPP
#define I_GLOBAL_PART_HPP
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
// Date:
// Oct. 2 2013
// ==================================================================
#include "IsoPETSc3D_System.hpp"
#include "metis.h"

class IGlobal_Part
{
  public:
    IGlobal_Part(){};
    virtual ~IGlobal_Part(){};

    virtual idx_t get_epart( s_int e ) const = 0;
    virtual idx_t get_npart( s_int n ) const = 0;

    virtual bool get_isMETIS() const = 0;
    virtual bool get_isDual() const = 0;
    virtual int get_dual_edge_ncommon() const = 0;
};

#endif
