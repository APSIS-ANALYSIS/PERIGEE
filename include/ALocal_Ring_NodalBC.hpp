#ifndef ALOCAL_RING_NODALBC_HPP
#define ALOCAL_RING_NODALBC_HPP
// ============================================================================
// ALocal_Ring_NodalBC.hpp
//
// Analysis-use, ring nodal indices.
//
// Author: Ingrid S. Lan
// Date Created: April 15 2021
// ============================================================================
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "HDF5_Reader.hpp"

class ALocal_Ring_NodalBC
{
  public:
    ALocal_Ring_NodalBC( const std::string &fileBaseName, const int &cpu_rank );

    virtual ~ALocal_Ring_NodalBC();

    // Get the type of ring nodal BC
    virtual int get_ringbc_type() const { return ringbc_type; }

    // Get global index of a Dirichlet node in the local partition
    // para node ranges [ 0, Num_LD ).
    virtual int get_LDN( const int &node ) const { return LDN[node]; }

    // Get the number of Dirichlet nodes in the local partition
    virtual int get_Num_LD() const { return Num_LD; }

    // Get the number of caps
    virtual int get_num_caps() const { return num_caps; }

    // Get the Dirichlet node's corresponding cap ID
    // para node ranges [ 0, Num_LD ).
    // output value is [ 0, num_caps ), where ID 0 corresponds to the inlet.
    virtual int get_cap_id( const int &node ) const { return local_cap_id[node]; }

    // get the Dirichlet node's outward normal vector components.
    // parameter node ranges [ 0, Num_LD )
    // comp=0 : x-component; comp=1 : y-component; comp=2 : z-component
    virtual double get_outvec( const int &node, const int &comp ) const
    { return outnormal[ local_cap_id[node] ]( comp ); }

    // get the Dirichlet node's rotation matrix for skew boundary conditions.
    // parameter node ranges [ 0, Num_LD )
    virtual Matrix_3x3 get_rotation_matrix( const int &node ) const
    { return Q[ local_cap_id[node] ]; }

    // determine whether a given index belongs to the LDN vector
    virtual bool is_inLDN( const int &ii) const
    { return VEC_T::is_invec(LDN, ii); }

  private:
    // type = 0 : all dof of ring nodes are set to be essential bc;
    // type = 1 : the dominant normal components of ring nodes are set to be essential bc;
    // type = 2 : the dominant normal & tangential components of ring nodes are set to be
    //            essential bc.
    int ringbc_type;

    // Number of caps (equals the number of inlets + the number of outlets)
    int num_caps;

    // Number of local Dirichlet nodes of the ring BC type
    int Num_LD;

    // If Num_LD > 0, store the global indices of the local Dirichlet nodes 
    // of the ring BC type
    std::vector<int> LDN;

    // If Num_LD > 0, store corresponding cap ID, which ranges [0, num_caps)
    // vector length is Num_LD
    std::vector<int> local_cap_id;

    // Each cap's unit normal vector
    // vector length is num_caps
    std::vector<Vector_3> outnormal;

    // Each cap's rotation matrix for skew boundary conditions
    // vector length is num_caps
    std::vector<Matrix_3x3> Q;
};

#endif
