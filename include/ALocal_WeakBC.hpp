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
    virtual int get_wall_model_type() const { return wall_model_type; }

    virtual int get_num_ele() const { return num_sur_ele; }

    // Return the local volume element index
    // 0 <= ee < get_num_ele()
    virtual int get_part_vol_ele_id( const int &ee ) const { return part_vol_ele_id[ee]; }

    // Return the face id of the surface with respect to its volumetric element
    // 0 <= ee < get_num_ele()
    virtual int get_ele_face_id( const int &ee ) const { return ele_face_id[ee]; }

    virtual void print_info() const;

  protected:
    // type = 0: the whole wall is set to be strongly no-slip BC
    // type = 1: all dof of wall nodes are set to be weakly essential BC;
    // type = 2: the normal direction of wall nodes (the rotated-x component)
    //           are set to be strongly no-slip, and the tanget direction 
    //           of wall nodes are set to be weakly essential BC.       
    int wall_model_type;

    // stores the number of surface elements.
    int num_sur_ele;

    // Store the local volume element id in this part
    std::vector<int> part_vol_ele_id;

    // Store the face id of volume element
    std::vector<int> ele_face_id;

    // ------------------------------------------------------------------------
    // Disallow default constructor
    ALocal_WeakBC() = delete;
};

#endif
