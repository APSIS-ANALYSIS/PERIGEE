#ifndef MESH_TET_HPP
#define MESH_TET_HPP
// ============================================================================
// Mesh_Tet.hpp
//
// This is the instantiation of the IMesh class for tetrahedral element,
// i.e. unstructured mesh.
//
// Author: Ju Liu
// Date: July 30 2023
// ============================================================================
#include "IMesh.hpp"

class Mesh_Tet : public IMesh
{
  public:
    Mesh_Tet(const int &in_nFunc, const int &in_nElem, const int &in_deg);

    virtual ~Mesh_Tet() = default;

    virtual void print_info() const;

    virtual int get_degree() const {return deg;}

    virtual int get_nFunc() const {return nFunc;}

    virtual int get_nElem() const {return nElem;}

    virtual int get_nLocBas() const {return nLocBas;}

  private:
    const int nFunc, nElem, deg;
   
    // ------------------------------------------------------------------------ 
    // In this class, nLocBas is determined by the element type, and since we
    // are restricted to tet element, the degree determines the value of it.
    // ------------------------------------------------------------------------ 
    int nLocBas;
};

#endif
