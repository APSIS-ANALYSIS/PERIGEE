#ifndef NODALBC_3D_LINEARELASTIC_HPP
#define NODALBC_3D_LINEARELASTIC_HPP
// ==================================================================
// NodalBC_3D_Linearelastic.hpp
//
// This is an instantiation of INodalBC for 3D linear elastic solid 
// mech problems.
//
// This class specifies the nodal strong boundary conditions.
//
// Date: May 18 2017
// ==================================================================
#include "INodalBC.hpp"

class NodalBC_3D_Linearelastic : public INodalBC
{
  public:
    NodalBC_3D_Linearelastic( const int &nFunc, const int &nFunc_x, 
        const int &nFunc_y, const int &nFunc_z, const int &bc_type);

    virtual ~NodalBC_3D_Linearelastic();

  private:
    NodalBC_3D_Linearelastic() {};

    // --------------------------------------------------------------
    // BC_type_1 : C0 Dirichlet on all 6 faces
    // --------------------------------------------------------------
    void BC_type_1( const int &nFunc_x, const int &nFunc_y, const int &nFunc_z );


    // --------------------------------------------------------------
    // BC_type_2 : No nodal BC
    // --------------------------------------------------------------
    void BC_type_2();


    // --------------------------------------------------------------
    // BC_type_3: C0 Dirichlet on LEFT, RIGHT, FRONT, BACK, BOTTOM.
    //            We assume the geometry is a cube.
    // --------------------------------------------------------------
    void BC_type_3( const int &nFunc_x, const int &nFunc_y, const int &nFunc_z );


    // --------------------------------------------------------------
    // BC_type_4: C0 Dirichlet on LEFT.
    //            We assume the geometry has a rectangular shape, such as
    //            the cube/beam, etc.
    // --------------------------------------------------------------
    void BC_type_4( const int &nFunc_x, const int &nFunc_y, const int &nFunc_z );


    // --------------------------------------------------------------
    // BC_type_5: C0-Periodic on Front-Back pair.
    //            This is used for the vessel wall model, such as
    //            the /input/geometry_3d_annular_tube.txt, which gives
    //            parametrization of an annular geometry using C0 NURBS
    // --------------------------------------------------------------
    void BC_type_5( const int &nFunc_x, const int &nFunc_y, const int &nFunc_z );


    // --------------------------------------------------------------
    // BC_type_6: C0-Dirichlet on Bottom
    //            We assume the geometry is a cube.
    // --------------------------------------------------------------
    void BC_type_6( const int &nFunc_x, const int &nFunc_y, const int &nFunc_z );


    // --------------------------------------------------------------
    // BC_type_7: C0-Dirichlet on Bottom & Top
    //            We assume the geometry is a cube.
    // --------------------------------------------------------------
    void BC_type_7( const int &nFunc_x, const int &nFunc_y, const int &nFunc_z );

    
    // --------------------------------------------------------------
    // BC_type_8: C0-Dirichlet on Top & BACK
    //            We assume the geometry is a cube.
    // --------------------------------------------------------------
    void BC_type_8( const int &nFunc_x, const int &nFunc_y, const int &nFunc_z );

    
    // --------------------------------------------------------------
    // BC_type_9: C0-Dirichlet on Top & LEFT 
    //            We assume the geometry is a cube.
    // --------------------------------------------------------------
    void BC_type_9( const int &nFunc_x, const int &nFunc_y, const int &nFunc_z );


    // --------------------------------------------------------------
    // BC_type_10: Dirichlet on bottom; on top, all nodes follow one
    //             master nodes.
    //             Used for the z-disp and z-velo in the tensile test
    // --------------------------------------------------------------
    void BC_type_10( const int &nFunc_x, const int &nFunc_y, const int &nFunc_z );

};


#endif
