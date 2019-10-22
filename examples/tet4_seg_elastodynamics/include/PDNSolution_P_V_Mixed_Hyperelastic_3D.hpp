#ifndef PDNSOLUTION_P_V_MIXED_HYPERELASTIC_3D_HPP
#define PDNSOLUTION_P_V_MIXED_HYPERELASTIC_3D_HPP
// ==================================================================
// PDNSolution_P_V_Mixed_Hyperelastic_3D.hpp
//
// This is the solution class for the pressure-velocity field in the 
// mixed formulation for hyperelasticity. This is part of the whole 
// solution U_P_V system.
//
// Author: Ju Liu
// Date: May 17 2017
// ==================================================================
#include "PDNSolution.hpp"
#include "IPLocAssem.hpp"

class PDNSolution_P_V_Mixed_Hyperelastic_3D : public PDNSolution
{
  public:
    PDNSolution_P_V_Mixed_Hyperelastic_3D(
        const class APart_Node * const &pNode,
        const class IPLocAssem * const &locassem_ptr,
        const class FEANode * const &fNode_ptr,
        const int &type );

    virtual ~PDNSolution_P_V_Mixed_Hyperelastic_3D();

    // Test case 1: -1
    void Init_test_1(const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

    // Test case 2: -2
    // Generate random numbers in the slots
    void Init_test_2(const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

    // Case 0 : 0
    // Generate full zero vector
    void Init_zero(const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );
};

#endif
