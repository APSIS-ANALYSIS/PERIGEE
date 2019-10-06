#ifndef PDNSOLUTION_NS_3D_HPP
#define PDNSOLUTION_NS_3D_HPP
// ==================================================================
// PDNSolution_NS_3D.hpp
// This is the class that we define the initial condition for the 
// three-dimensional incompressible Navier-Stokes equations
//
// Date: March 9 2015
// ==================================================================
#include <cmath>
#include "Sys_Tools.hpp"
#include "IPLocAssem.hpp"
#include "FEANode.hpp"
#include "PDNSolution.hpp"
#include "APart_Node.hpp"
#include "IALocal_BC.hpp"

class PDNSolution_NS_3D : public PDNSolution
{
  public:
    PDNSolution_NS_3D(const class APart_Node * const &pNode,
        const class IALocal_BC * const &lbc,
        const class IPLocAssem * const &locassem_ptr,
        const class FEANode * const &fNode_ptr,
        const int &type );

    virtual ~PDNSolution_NS_3D();

    // Initial condition setting
    // 0. Test case
    void Init_test(const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );
    
    // 1. Test case
    void Init_test_2(const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );
    
    // 2. Three-dimensional driven cavity
    void Init_3D_Driven_Cavity(const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

    // 3. Three-dimensional Poiseuille in a pipe
    void Init_3D_Poiseuille_Pipe( const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

    // 4. Three-dimensional pressure drivenPoiseuille in a pipe
    void Init_3D_Pressure_Pipe( const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

    // 5. Three-dimensional solution vector generated by random number
    void Init_random( const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

    // 6. Three-dimensional pressure driven in a coronary
    void Init_3D_Pressure_Coronary_3patch( 
        const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

    // 7. Three-dimensional dirichlet inflow in coronary
    void Init_3D_Dir_inflow_Coronary_3patch( 
        const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );

    // 8. Three-dimensional pressure driven in a coronary
    void Init_3D_Pressure_Coronary_107patch( 
        const APart_Node * const &pNode_ptr,
        const IPLocAssem * const &locassem_ptr,
        const FEANode * const &fNode_ptr );


};

#endif
