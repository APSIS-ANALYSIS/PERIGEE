#ifndef PDNSOLUTION_U_MIXED_HYPERELASTIC_3D_HPP
#define PDNSOLUTION_U_MIXED_HYPERELASTIC_3D_HPP
// ==================================================================
// PDNSolution_U_Mixed_Hyperelastic_3D.hpp
//
// This is the solution class for the displacement field in the mixed
// formulation for hyperelasticity. This is part of the whole solution
// V-P-U system.
//
// Author: Ju Liu
// Date: May 17 2017
// ==================================================================

#include "PDNSolution.hpp"
#include "IPLocAssem.hpp"

class PDNSolution_U_Mixed_Hyperelastic_3D : public PDNSolution
{
  public:
    PDNSolution_U_Mixed_Hyperelastic_3D(
        const class APart_Node * const &pNode,
        const class IPLocAssem * const &locassem_ptr,
        const class FEANode * const &fNode_ptr,
        const int &type );

    virtual ~PDNSolution_U_Mixed_Hyperelastic_3D();

    void Init_test_1(const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

    void Init_zero(const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );
};

#endif
