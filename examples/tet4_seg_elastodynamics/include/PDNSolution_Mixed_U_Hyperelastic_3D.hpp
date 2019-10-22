#ifndef PDNSOLUTION_MIXED_U_HYPERELASTIC_3D_HPP
#define PDNSOLUTION_MIXED_U_HYPERELASTIC_3D_HPP
// ==================================================================
// PDNSolution_Mixed_U_Hyperelastic_3D.hpp
// This is the solution class for the fully coupled mixed formulation
// with u-v kinematics.
//
// Author: Ju Liu
// Date: Dec. 11 2016
// ==================================================================
#include "PDNSolution.hpp"
#include "IPLocAssem.hpp"
#include "Math_Tools.hpp"

class PDNSolution_Mixed_U_Hyperelastic_3D : public PDNSolution
{
  public:
    PDNSolution_Mixed_U_Hyperelastic_3D( 
        const class APart_Node * const &pNode,
        const class IPLocAssem * const &locassem_ptr,
        const class FEANode * const &fNode_ptr,
        const int &type );

    virtual ~PDNSolution_Mixed_U_Hyperelastic_3D();

    void Init_test_1(const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

    // Generate random entries.
    void Init_test_2(const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );
  
    // Generate zero disp and pres, velo = sin(x)
    void Init_test_3( const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );


    // case 0: generate full zero vector 
    void Init_zero(const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

    // case 1: generate zero for the disp and pres slots, 
    //         uses specific values for the velo slots.
    //         default velo value is [pi ; 0.0 ; 0.0].
    void Init_zero_u_p_one_v( const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

    // case 2: generate zero for the disp and pres slots, 
    //         velocity in x dir is 10 m/sec * z / 6.0 m
    //         This one is used for testing A. Gil's thick column
    //         oscillation benchmark
    void Init_vx_linear( const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

};

#endif
