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
#include "ALocal_InflowBC.hpp"

class PDNSolution_Mixed_UPV_3D : public PDNSolution
{
  public:
    PDNSolution_Mixed_UPV_3D( const APart_Node * const &pNode,
        const FEANode * const &fNode_ptr,
        const int &type, const bool &isprint = true );

    PDNSolution_Mixed_UPV_3D( const APart_Node * const &pNode,
        const FEANode * const &fNode_ptr,
        const ALocal_InflowBC * const &infbc,
        const int &type, const bool &isprint = true );

    PDNSolution_Mixed_UPV_3D( const APart_Node * const &pNode,
        const int &type, const bool &isprint = true );

    virtual ~PDNSolution_Mixed_UPV_3D();

    // --------------------------------------------------------------
    // case 0: generate full zero vector 
    // --------------------------------------------------------------
    void Init_zero( const APart_Node * const &pNode_ptr );

    // --------------------------------------------------------------
    // case 1: generate flow parabolic for an arbitrary inlet
    //         face with unit flow rate. This solution is used as the
    //         base solution and will be scaled by the flow rate to 
    //         adjust the inlet Dirichlet nodal values.
    // --------------------------------------------------------------
    void Init_flow_parabolic( const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr,
        const ALocal_InflowBC * const &infbc );
  
    // --------------------------------------------------------------
    // case 2: generate pressure with a prescribed value
    //         for whole continuum body, and zero velocity and displacement
    // --------------------------------------------------------------
    void Init_pressure( const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr );

  private:
    const bool is_print;
};

#endif
