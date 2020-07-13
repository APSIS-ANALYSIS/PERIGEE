#ifndef PDNSOLUTION_HEATEQN_HPP
#define PDNSOLUTION_HEATEQN_HPP
// ==================================================================
// PDNSolution_heatEqn.hpp
// This is the class that we define initial solution for heat equation
//
// Date: Nov. 26 2013
// ==================================================================
#include "Sys_Tools.hpp"
#include "PDNSolution.hpp"
#include "APart_Node.hpp"
#include "IALocal_BC.hpp"

class PDNSolution_heatEqn : public PDNSolution
{
  public:
    PDNSolution_heatEqn( const class APart_Node * const &pNode,
       const class IALocal_BC * const &lbc, int type );
    virtual ~PDNSolution_heatEqn();

    // Initial solution setting
    // case 0: all zero for initial temperature
    void Init_ZeroTemp( const class IALocal_BC * const &LBC );

    // case 1: interior values 1.0, boundary values 0.0
    void Init_OneTemp( const class IALocal_BC * const &LBC );
};
#endif
