#ifndef MESH_TET4_HPP
#define MESH_TET4_HPP
// ==================================================================
// Mesh_Tet4.hpp
//
// This is the instantiation of the IMesh class for 4-node tetrahedral
// element, i.e. linear tets, mesh. An unstructured mesh.
//
// Date: Jan. 11 2017
// ==================================================================
#include "IMesh.hpp"

class Mesh_Tet4 : public IMesh
{
  public:
    Mesh_Tet4(const int &in_nFunc, const int &in_nElem);

    virtual ~Mesh_Tet4();

    virtual void print_mesh_info() const;

    virtual int get_s_degree() const {return 1;}
    virtual int get_t_degree() const {return 1;}
    virtual int get_u_degree() const {return 1;}

    virtual int get_nFunc() const {return nFunc;}

    virtual int get_nElem() const {return nElem;}

    virtual int get_nLocBas() const {return 4;}

  private:
    const int nFunc, nElem;
};

#endif
