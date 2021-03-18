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
#include <complex_bessel.h>
#include "PDNSolution.hpp"
#include "FEANode.hpp"
#include "ALocal_Inflow_NodalBC.hpp"

class PDNSolution_NS : public PDNSolution
{
  public:
    PDNSolution_NS( const APart_Node * const &pNode,
        const FEANode * const &fNode_ptr,
        const ALocal_Inflow_NodalBC * const &infbc,
        const int &type, const bool &isprint = true );

    PDNSolution_NS( const APart_Node * const &pNode, 
        const int &type, const bool &isprint = true );

    // ==== WOMERSLEY CHANGES BEGIN ====
    PDNSolution_NS( const APart_Node * const &pNode, 
        const FEANode * const &fNode_ptr,
        const double &rho,
        const double &vis_mu,
        const int &type, const bool &isprint = true );
    // ==== WOMERSLEY CHANGES END ====

    virtual ~PDNSolution_NS();

    // case 0 : generate full zero solution
    void Init_zero( const APart_Node * const &pNode_ptr );

    // case 1: generate flow parabolic for an arbitrary inlet face
    //         with the unit flow rate.
    void Init_flow_parabolic( const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr,
        const ALocal_Inflow_NodalBC * const &infbc );

    // case 2: generate flow parabolic for an arbitrary inlet face
    //         with the unit flow rate.
    void Init_pipe_parabolic( const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr );


    // ==== WOMERSLEY CHANGES BEGIN ====

    // case 3 : generate full Womersley solution
    void Init_womersley( const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr,
        const double &rho, const double &vis_mu );

    // case 4 : generate full Womersley dot solution
    void Init_womersley_dot( const APart_Node * const &pNode_ptr, 
        const FEANode * const &fNode_ptr,
        const double &rho, const double &vis_mu );

    // ==== WOMERSLEY CHANGES END ====

  private:
    const bool is_print;
};

#endif
