#ifndef PDNSOLUTION_DISP_3D_HPP
#define PDNSOLUTION_DISP_3D_HPP
// ==================================================================
// PDNSolution_Disp_3D.hpp
// This is a class that stores the Displacement solution vector for
// 3D Elasticity.
//
// Author: Ju Liu
// Date: Jan 16 2018
// ==================================================================
#include "PDNSolution.hpp"
#include "FEANode.hpp"

class PDNSolution_Disp_3D : public PDNSolution
{
  public:
    PDNSolution_Disp_3D( 
        const APart_Node * const &pNode,
        const FEANode * const &fNode_ptr,
        const int &type );

    ~PDNSolution_Disp_3D();

    void Init_zero_solu( const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr );
};

#endif
