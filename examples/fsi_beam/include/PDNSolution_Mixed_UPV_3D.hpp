#ifndef PDNSOLUTION_MIXED_UPV_3D_HPP
#define PDNSOLUTION_MIXED_UPV_3D_HPP
// ==================================================================
// PDNSolution_Mixed_UPV_3D.hpp
// This is the solution class for the fully coupled mixed formulation
// with u-v kinematics. It contains 7 dofs:
//  [ disp_x, disp_y, disp_z, pres, velo_x, velo_y, velo_z ].
//
// Author: Ju Liu
// Date: Aug. 11 2017
// ==================================================================
#include "PDNSolution.hpp"
#include "FEANode.hpp" 
#include "Math_Tools.hpp"

class PDNSolution_Mixed_UPV_3D : public PDNSolution
{
  public:
    PDNSolution_Mixed_UPV_3D( 
        const class APart_Node * const &pNode,
        const class FEANode * const &fNode_ptr,
        const int &type );

    virtual ~PDNSolution_Mixed_UPV_3D();

    // case 0: generate full zero vector 
    void Init_zero(const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr );

    // case 1: generate inflow profile 
    void Init_inflow(const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr );
};

#endif
