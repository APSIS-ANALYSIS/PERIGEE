#ifndef ALOCAL_NODALBC_P2P1_3D_HPP
#define ALOCAL_NODALBC_P2P1_3D_HPP
// ==================================================================
// ALocal_NodalBC_P2P1_3D.hpp
//
// This class generate the working Nodal BC and LID info for 3D
// P2/P1 Taylor-Hood element. The first degree of freedom is assumed
// to be the pressure variable, and does not have nodal BC. The next
// three dof are for the velocity or displacement.
//
// Note: Periodic BC is not considered.
//
// Author: Ju Liu
// Date: Feb. 15 2018
// ==================================================================
#include "ALocal_NodalBC_2x2Block.hpp"
#include "APart_Node.hpp"

class ALocal_NodalBC_P2P1_3D : public ALocal_NodalBC_2x2Block
{
  public:
    ALocal_NodalBC_P2P1_3D( const std::string &fileBaseName,
        const int &cpu_rank, const APart_Node * const &pNode_pres );

    virtual ~ALocal_NodalBC_P2P1_3D();

    virtual void print_info() const;

    // LID array
    virtual int get_LID_0( const int &node ) const {return LID_p[node];}

    virtual int get_LID_1( const int &dof_idx, const int &node ) const
    {return LID_u[dof_idx * nlgnode + node];}
    
    // Local Dirichlet Node index
    virtual int get_LDN_0( const int &node ) const {return LDN[node];}

    virtual int get_LDN_1( const int &dof_idx, const int &node ) const
    {return LDN[ LD_offset_u[dof_idx] + node ]; }

    // Number of Local Dirichlet nodes
    virtual int get_Num_LD_0() const {return Num_LD_p;}

    virtual int get_Num_LD_1( const int &dof_idx ) const
    {return Num_LD_u[dof_idx];}

  private:
    int nlgnode; // number of local and ghost node in P2 configuration
    int Num_LD_p;
    int Num_LD_u[3];
    int LD_offset_u[3];

    std::vector<int> LID_p, LID_u, LDN;
};

#endif
