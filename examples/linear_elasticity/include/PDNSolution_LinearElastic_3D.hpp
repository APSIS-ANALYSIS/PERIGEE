#ifndef PDNSOLUTION_LINEARELASTIC_3D_HPP
#define PDNSOLUTION_LINEARELASTIC_3D_HPP
// ==================================================================
// PDNSolution_LinearElastic_3D.hpp
// This is the class that we define the initial condition for the 3D
// Linear Elastic problems.
//
// Note: If the solution is for static calculations, make sure the
// Dirichlet boundary condition is correctly specified in this solution
// vector.
//
// Author: Ju Liu
// Date: May 7 2017
// ==================================================================
#include "PDNSolution.hpp"
#include "IPLocAssem.hpp"

class PDNSolution_LinearElastic_3D : public PDNSolution
{
  public:
    PDNSolution_LinearElastic_3D( const class APart_Node * const &pNode,
        const class IPLocAssem * const &locassem_ptr,
        const class FEANode * const &fNode_ptr,
        const int &type );
    
    ~PDNSolution_LinearElastic_3D();

    void Init_test_1( const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

    void Init_zero_solu( const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );
};

#endif
