#ifndef PDNSOLUTION_NS_HPP
#define PDNSOLUTION_NS_HPP
// ==================================================================
// PDNSolution_NS.hpp
//
// This is the solution class for NS Solver. 
// 
// The solution vector has 4 dofs, pressure, u, v, w. 
//
// Author: Ju Liu
// Date Created: Jan 7 2020
// ==================================================================
#include "PDNSolution.hpp"
#include "FEANode.hpp"
#include "ALocal_InflowBC.hpp"

class PDNSolution_NS : public PDNSolution
{
  public:
    PDNSolution_NS( const APart_Node * const &pNode,
        const FEANode * const &fNode_ptr,
        const ALocal_InflowBC * const &infbc,
        const int &type, const bool &isprint = true );

    PDNSolution_NS( const APart_Node * const &pNode, 
        const int &type, const bool &isprint = true );

    virtual ~PDNSolution_NS();

    // case 0 : generate full zero solution
    void Init_zero( const APart_Node * const &pNode_ptr );

    // case 1: generate flow parabolic for an arbitrary inlet face
    //         with the unit flow rate.
    void Init_flow_parabolic( const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr,
        const ALocal_InflowBC * const &infbc );

    // case 2: generate flow parabolic for an arbitrary inlet face
    //         with the unit flow rate.
    void Init_pipe_parabolic( const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr );

    // v = v_max * (1 - r^3 / R^3), v_max = 5/3 * v_ave
    void Init_flow_cubic( const APart_Node * const &pNode_ptr, 
        const FEANode * const &fNode_ptr,
        const ALocal_InflowBC * const &infbc );

  private:
    const bool is_print;
};

#endif
