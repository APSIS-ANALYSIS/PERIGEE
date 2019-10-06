#ifndef PDNSOLUTION_LINEARELASTIC_2D_HPP
#define PDNSOLUTION_LINEARELASTIC_2D_HPP
// ==================================================================
// PDNSolution_LinearElastic_2D.hpp
// This is the class that we define the initial condition for the 2D
// Linear Elastic problems.
//
// Note: If the solution is for static calculations, make sure the
// Dirichlet boundary condition is correctly specified in this solution
// vector.
//
// Author: Ju Liu
// Date: Sept. 12 2016
// ==================================================================

#include "PDNSolution.hpp"
#include "IPLocAssem.hpp"

class PDNSolution_LinearElastic_2D : public PDNSolution
{
  public:
    PDNSolution_LinearElastic_2D( const class APart_Node * const &pNode,
        const class IPLocAssem * const &locassem_ptr,
        const class FEANode * const &fNode_ptr,
        const int &type );
    
    ~PDNSolution_LinearElastic_2D();

    void Init_zero_solu( const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );
};


#endif
