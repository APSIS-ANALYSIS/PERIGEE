#ifndef PDNSOLUTION_TET4_NS_3D_HPP
#define PDNSOLUTION_TET4_NS_3D_HPP
// ==================================================================
// PDNSolution_Tet4_NS_3D.hpp
//
// This is the solution class for the NS 3D solver based on Tet4
// elements.
//
// Author: Ju Liu
// Date Created: June 19 2017
// ==================================================================
#include "PDNSolution.hpp"
#include "FEANode.hpp"
#include "ALocal_Inflow_NodalBC.hpp"

class PDNSolution_Tet4_NS_3D : public PDNSolution
{
  public:
    PDNSolution_Tet4_NS_3D(  const APart_Node * const &pNode,
        const FEANode * const &fNode_ptr,
        const ALocal_Inflow_NodalBC * const &infbc,
        const double &in_vol_rate_Q,
        const int &type, const bool &isprint = true );

    PDNSolution_Tet4_NS_3D(  const APart_Node * const &pNode,
        const FEANode * const &fNode_ptr,
        const int &type, const bool &isprint = true );

    virtual ~PDNSolution_Tet4_NS_3D();

    // case 0 : generate full zero solution
    void Init_zero(const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr );

    // case 1 : generate inflow parabolic for a tube with raidus 2cm,
    //          and center axis is z-axis.
    void Init_inflow_parabolic_2cm(const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr,
        const ALocal_Inflow_NodalBC * const &infbc,
        const double &vol_rate_Q );

    // case 2 : generate flow parabolic for an arbitrary inlet face with
    //          the prescribed volumetric flow rate.
    void Init_flow_parabolic(
        const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr,
        const ALocal_Inflow_NodalBC * const &infbc,
        const double &vol_rate_Q );
    
    // case 3 : generate plug flow for a tube with radius 2cm and 
    //          30cm long. Pressure is zero.
    void Init_flow_plug(
        const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr,
        const ALocal_Inflow_NodalBC * const &infbc,
        const double &vol_rate_Q );
    
    // case 5 : generate inflow parabolic for a tube with raidus 0.6cm,
    //          and center axis is z-axis.
    void Init_inflow_parabolic_z_0d6cm(
        const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr,
        const ALocal_Inflow_NodalBC * const &infbc,
        const double &vol_rate_Q );

  private:
    const bool is_print;
};

#endif
