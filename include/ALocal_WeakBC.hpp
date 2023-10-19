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

    // Get the global index of a Dirichlet node in the local partition
    // para node ranges [ 0, num_LD )
    virtual int get_LDN( const int &node ) const { return LDN[node]; }

    // Get the number of Dirichlet nodes in the local partition
    virtual int get_num_LD() const { return num_LD; }

    // Get the number of surface cells in the local partition
    virtual int get_num_local_cell() const { return num_local_cell; }


  protected:
    // type = 0: the whole wall is set to be strongly no-slip BC
    // type = 1: all dof of wall nodes are set to be weakly essential BC;
    // type = 2: the normal direction of wall nodes (the rotated-x component)
    //           are set to be strongly no-slip, and the tanget direction 
    //           of wall nodes are set to be weakly essential BC.       
    int weakbc_type;

    // Number of local Dirichlet nodes of the weak BC type
    int num_LD;

    // Number of local cell
    int num_local_cell;

    // Number of nodes in a cell
    int cell_nLocBas;

    // Local cell's IEN array
    std::vector<int> local_cell_ien {};

    // Local cell coordinates, the length is 3 x num_LD
    std::vector<double> LD_xyz {};

    // Store the global indices of the local Dirichlet nodes of the weak BC type,
    // the length is num_LD
    std::vector<int> LDN {};

    // If weakbc_type = 2, store the rotation matrices at nodes,
    // the length is num_LD
    std::vector<Tensor2_3D> Q {};

    // ------------------------------------------------------------------------
    // Disallow default constructor
    ALocal_WeakBC() = delete;

};

#endif