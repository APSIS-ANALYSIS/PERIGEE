#ifndef PDNSOLUTION_MIXED_LINEARELASTIC_3D_HPP
#define PDNSOLUTION_MIXED_LINEARELASTIC_3D_HPP
// ==================================================================
// PDNSolution_Mixed_LinearElastic_3D.hpp
// This is the class that we define the solution for 3D Elasticity
// in mixed formulation.
//
// Author: Ju Liu
// Date: Dec. 22 2017
// ==================================================================
#include "PDNSolution.hpp"
#include "IPLocAssem.hpp"

class PDNSolution_Mixed_LinearElastic_3D : public PDNSolution
{
  public:
    PDNSolution_Mixed_LinearElastic_3D( 
        const APart_Node * const &pNode,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr,
        const int &type );

    ~PDNSolution_Mixed_LinearElastic_3D();

    void Init_zero_solu( const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

};

#endif
