#ifndef PDNSOLUTION_HEAT_HPP
#define PDNSOLUTION_HEAT_HPP
// ==================================================================
// PDNSolution_Heat.hpp
//
// This is the solution class for Heat Solver. 
// 
// The solution vector has 1 dofs.
//
// Author: Ju Liu
// Date Created: May 12 2020
// ==================================================================
#include "PDNSolution.hpp"
#include "FEANode.hpp"
#include "ALocal_InflowBC.hpp"

class PDNSolution_Heat : public PDNSolution
{
  public:
    PDNSolution_Heat( const APart_Node * const &pNode,
        const FEANode * const &fNode_ptr,
        const ALocal_InflowBC * const &infbc,
        const int &type, const bool &isprint = true );

    PDNSolution_Heat( const APart_Node * const &pNode, 
        const int &type, const bool &isprint = true );

    virtual ~PDNSolution_Heat();

    // case 0 : generate full zero solution
    void Init_zero( const APart_Node * const &pNode_ptr );

    // case 1: generate flow parabolic for an arbitrary inlet face
    //         with the unit flow rate.
    void Init_flow_parabolic( const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr,
        const ALocal_InflowBC * const &infbc );

  private:
    const bool is_print;
};

#endif
