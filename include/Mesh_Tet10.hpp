#ifndef MESH_TET10_HPP
#define MESH_TET10_HPP
// ==================================================================
// Mesh_Tet10.hpp
//
// This is the instantiation of the IMesh class for 10-node tetrahedral
// element, i.e. quadratic tets mesh. An unstructured mesh.
//
// Author: Ju Liu
// Date: Jan. 04 2020
// ==================================================================
#include "IMesh.hpp"

class Mesh_Tet10 : public IMesh
{
  public:
    Mesh_Tet10(const int &in_nFunc, const int &in_nElem);

    virtual ~Mesh_Tet10();

    virtual void print_mesh_info() const;

    virtual int get_s_degree() const {return 2;}
    virtual int get_t_degree() const {return 2;}
    virtual int get_u_degree() const {return 2;}

    virtual int get_nFunc() const {return nFunc;}

    virtual int get_nElem() const {return nElem;}

    virtual int get_nLocBas() const {return 10;}

  private:
    const int nFunc, nElem;
};

#endif
