#ifndef ELEMBC_3D_LINEARELASTIC_HPP
#define ELEMBC_3D_LINEARELASTIC_HPP
// ==================================================================
// ElemBC_3D_Linearelastic.hpp
//
// This is the boundary element face specification for 3D linear 
// elastic solid problems.
//
// Date: May 18 2017
// ==================================================================
#include "IElemBC.hpp"

class ElemBC_3D_Linearelastic : public IElemBC
{
  public:
    ElemBC_3D_Linearelastic( const int &nElem_x, const int &nElem_y, 
        const int &nElem_z, const int &bc_type );

    virtual ~ElemBC_3D_Linearelastic();

  private:
    ElemBC_3D_Linearelastic() {};

    // --------------------------------------------------------------
    // BC_test : This is a routine for debugging. bc_type = 0 
    // --------------------------------------------------------------
    void BC_test(const int &nElem_x, const int &nElem_y, const int &nElem_z);
    
    // --------------------------------------------------------------
    // BC_type_1 : No boundary has face integration
    // --------------------------------------------------------------
    void BC_type_1();

    // --------------------------------------------------------------
    // BC_type_2 : For the 3D geometry, all 6 faces need a boundary 
    //             integral. 
    // --------------------------------------------------------------
    void BC_type_2(const int &nElem_x, const int &nElem_y, const int &nElem_z);

    // --------------------------------------------------------------
    // BC_type_3 : For the 3D geometry, top face requires surface 
    //             integration. 
    // --------------------------------------------------------------
    void BC_type_3(const int &nElem_x, const int &nElem_y, const int &nElem_z);

    // --------------------------------------------------------------
    // BC_type_4 : For the 3D geometry, bottom face requires surface 
    //             integration. 
    // --------------------------------------------------------------
    void BC_type_4(const int &nElem_x, const int &nElem_y, const int &nElem_z);

    // --------------------------------------------------------------
    // BC_type_5 : For the 3D geometry, bottom, top, front, back, and 
    //             right face requires surface integration.
    //             I use this in combination with NodalBC type 4 for
    //             cantilever beam problmes. 
    // --------------------------------------------------------------
    void BC_type_5(const int &nElem_x, const int &nElem_y, const int &nElem_z);

    // --------------------------------------------------------------
    // BC_type_6 : For the 3D geometry, bottom, top, left, and right 
    //             face requires surface integration.
    //             I use this in combination with NodalBC type 5 for
    //             vessel wall model using C0 NURBS parametrization. 
    // --------------------------------------------------------------
    void BC_type_6(const int &nElem_x, const int &nElem_y, const int &nElem_z);
    
    // --------------------------------------------------------------
    // BC_type_7 : For the 3D geometry, top, front, back, left, and right 
    //             face requires surface integration.
    //             I use this in combination with NodalBC type 6.
    // --------------------------------------------------------------
    void BC_type_7(const int &nElem_x, const int &nElem_y, const int &nElem_z);

};

#endif
