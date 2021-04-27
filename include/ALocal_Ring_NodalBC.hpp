#ifndef ALOCAL_RING_NODALBC_HPP
#define ALOCAL_RING_NODALBC_HPP
// ==================================================================
// ALocal_Ring_NodalBC.hpp
//
// Analysis-use, ring nodal indices.
//
// Author: Ingrid S. Lan
// Date Created: April 15 2021
// ==================================================================
#include "Vector_3.hpp"
#include "HDF5_Reader.hpp"

class ALocal_Ring_NodalBC
{
  public:
    ALocal_Ring_NodalBC( const std::string &fileBaseName, const int &cpu_rank );

    virtual ~ALocal_Ring_NodalBC();

    // Get global index of a Dirichlet node in the local partition
    // para node ranges [ 0 , Num_LD ).
    virtual int get_LDN( const int &node ) const { return LDN[node]; }

    // Get the number of Dirichlet nodes in the local partition
    virtual int get_Num_LD() const { return Num_LD; }

    // get the Dirichlet node's dominant component index, where the
    // dominant component of the unit normal vector is the largest in
    // absolute value
    // parameter node ranges [ 0, Num_LD )
    // output value is 0, 1, or 2.
    virtual int get_dominant_comp( const int &node ) const
    { return dominant_comp[ cap_id[node] ]; }

    // get the Dirichlet node's outward normal vector components.
    // parameter node ranges [ 0, Num_LD )
    // comp=0 : x-component; comp=1 : y-component; comp=2 : z-component
    virtual double get_outvec( const int &node, const int &comp ) const
    { return outnormal[ cap_id[node] ]( comp ); }

    // determine whether a given index belongs to the LDN vector
    virtual bool is_inLDN( const int &ii) const
    { return VEC_T::is_invec(LDN, ii); }

  private:
    // Number of local Dirichlet nodes of the ring BC type
    int Num_LD;

    // Global indices of the local Dirichlet nodes of the ring BC type
    std::vector<int> LDN;

    // Number of caps (equals the number of inlets + the number of outlets)
    int num_caps;

    // If Num_LD > 0, store corresponding cap ID, which ranges [0, num_caps)
    // vector length is Num_LD
    std::vector<int> cap_id;

    // Dominant component index of each cap's unit normal vector: 0, 1, or 2
    // vector length is num_caps
    std::vector<int> dominant_comp;

    // Each cap's unit normal vector,
    // vector length is num_caps
    std::vector<Vector_3> outnormal;
};

#endif
