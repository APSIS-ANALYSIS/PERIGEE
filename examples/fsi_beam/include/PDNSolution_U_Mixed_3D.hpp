#ifndef PDNSOLUTION_U_MIXED_3D_HPP
#define PDNSOLUTION_U_MIXED_3D_HPP
// ==================================================================
// PDNSolution_U_Mixed_3D.hpp
//
// This is the solution class for the displacement field in the mixed
// formulation for hyperelasticity. This class is designed for the
// displacement part of the whole solution in the U-P-V system.
//
// Author: Ju Liu
// Date Created: Aug 11 2017
// ==================================================================
#include "PDNSolution.hpp"
#include "FEANode.hpp"

class PDNSolution_U_Mixed_3D : public PDNSolution
{
  public:
    PDNSolution_U_Mixed_3D(
        const class APart_Node * const &pNode,
        const class FEANode * const &fNode_ptr,
        const int &type );

    virtual ~PDNSolution_U_Mixed_3D();

    // Test case: -1, user specified vector value
    void Init_test_1(const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr );

    // Test case: -2 random numbers in all slots
    void Init_test_2(const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr );

    // Case 0: Full zero vector
    void Init_zero(const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr );
};

#endif
