#ifndef ALOCAL_WEAKBC_HPP
#define ALOCAL_WEAKBC_HPP
// ============================================================================
// ALocal_WeakBC.hpp
//
// FEM-analysis-use Local subdomain's Weak enforced no-slip Boundary Conditions.
// It is like a combination of ringBC and EBC.
//
// Author: Xuanming Huang
// Date Created: Oct. 18th 2023
// ============================================================================
#include "HDF5_Reader.hpp"

class ALocal_WeakBC
{
  public:
    ALocal_WeakBC( const std::string &fileBaseName, const int &cpu_rank );

    virtual ~ALocal_WeakBC(){};

    // Get the type of weak enforced Dirichlet BC
    virtual int get_weakbc_type() const { return weakbc_type; }



  protected:
    // type = 0: the whole wall is set to be strongly no-slip BC
    // type = 1: all dof of wall nodes are set to be weakly essential BC;
    // type = 2: the normal direction of wall nodes (the rotated-x component)
    //           are set to be strongly no-slip, and the tanget direction 
    //           of wall nodes are set to be weakly essential BC.       
    int weakbc_type;


    // ------------------------------------------------------------------------
    // Disallow default constructor
    ALocal_WeakBC() = delete;

};

#endif