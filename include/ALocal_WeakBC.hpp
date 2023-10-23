#ifndef ALOCAL_WEAKBC_HPP
#define ALOCAL_WEAKBC_HPP
// ============================================================================
// ALocal_WeakBC.hpp
//
// FEM-analysis-use Local subdomain's Weak enforced no-slip Boundary Conditions.
//
// Author: Xuanming Huang
// Date Created: Oct. 18th 2023
// ============================================================================
#include "HDF5_Reader.hpp"

class ALocal_WeakBC
{
  public:
    ALocal_WeakBC( const std::string &fileBaseName, const int &cpu_rank );

    virtual ~ALocal_WeakBC();

    // Get the type of weak enforced Dirichlet BC
    virtual int get_weakbc_type() const { return weakbc_type; }

    virtual int get_num_weak_boundary() const { return num_weak_boundary; }

    virtual int get_C_bI() const { return C_bI; }

    virtual std::vector<int> get_part_vol_ele_id(const int &ii) { return part_vol_ele_id[ii]; }

    virtual std::vector<int> get_ele_face_id(const int &ii) { return ele_face_id[ii]; }

  protected:
    // type = 0: the whole wall is set to be strongly no-slip BC
    // type = 1: all dof of wall nodes are set to be weakly essential BC;
    // type = 2: the normal direction of wall nodes (the rotated-x component)
    //           are set to be strongly no-slip, and the tanget direction 
    //           of wall nodes are set to be weakly essential BC.       
    int weakbc_type;

    // Number of boundaries that assigned weak BC
    int num_weak_boundary;

    // Coefficient for weak BC
    int C_bI;

    // num_sur_ele[ii] gives the ii-th weakbc's surface element number.
    std::vector<int> num_sur_ele {};

    // Store the local volume element id this this part
    std::vector< std::vector<int> > part_vol_ele_id;

    // Store the face id of volume element
    std::vector< std::vector<int> > ele_face_id;

    // num_sur_node[ii] gives the ii-th weakbc's surface node number.
    std::vector<int> num_sur_node {};

    // Store the global id of surface node
    std::vector< std::vector<int> > sur_node_id;

    // Store the rotation matrices at nodes
    std::vector< std::vector<Tensor2_3D> > Q;

    // ------------------------------------------------------------------------
    // Disallow default constructor
    ALocal_WeakBC() = delete;

};

#endif