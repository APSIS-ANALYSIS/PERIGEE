#ifndef PDNSOLUTION_TET4_ALE_NS_3D_HPP
#define PDNSOLUTION_TET4_ALE_NS_3D_HPP
// ==================================================================
// PDNSolution_Tet4_ALE_NS_3D.hpp
//
// This is the solution class for the ALE formulation for NS 3D 
// solver based on Tet4 elements.
// 
// The solution vector has 7 dofs. The first three are the mesh motion
// dofs.
//
// Author: Ju Liu
// Date Created: Aug. 5 2017
// ==================================================================
#include "PDNSolution.hpp"
#include "FEANode.hpp"
#include "ALocal_Inflow_NodalBC.hpp"

class PDNSolution_Tet4_ALE_NS_3D : public PDNSolution
{
  public:
    PDNSolution_Tet4_ALE_NS_3D( const APart_Node * const &pNode,
        const FEANode * const &fNode_ptr,
        const ALocal_Inflow_NodalBC * const &infbc,
        const int &type, const bool &isprint = true );

    PDNSolution_Tet4_ALE_NS_3D( const APart_Node * const &pNode, 
        const int &type, const bool &isprint = true );

    virtual ~PDNSolution_Tet4_ALE_NS_3D();

    // case 0 : generate full zero solution
    void Init_zero( const APart_Node * const &pNode_ptr );

    // case 1: generate flow parabolic for an arbitrary inlet face
    //         with the unit flow rate.
    void Init_flow_parabolic(
        const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr,
        const ALocal_Inflow_NodalBC * const &infbc );

  private:
    const bool is_print;
};

#endif
