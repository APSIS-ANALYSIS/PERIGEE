#ifndef MESH_FEM_HPP
#define MESH_FEM_HPP
// ==================================================================
// Mesh_FEM.hpp
//
// This is the instantiation of the IMesh class for general element.
//
// Date: Nov. 17 2017
// ==================================================================
#include "IMesh.hpp"

class Mesh_FEM : public IMesh
{
  public:
    // Constructor
    // Assumes the polyminal degree is uniform
    Mesh_FEM(const int &in_nfunc, const int &in_nelem, 
        const int &in_nlocbas, const int &in_deg);

    virtual ~Mesh_FEM() = default;

    virtual void print_info() const;

    virtual int get_s_degree() const {return sdeg;}
    virtual int get_t_degree() const {return tdeg;}
    virtual int get_u_degree() const {return udeg;}

    virtual int get_nFunc() const {return nFunc;}

    virtual int get_nElem() const {return nElem;}

    virtual int get_nLocBas() const {return nLocBas;}

  private:
    const int nFunc, nElem, nLocBas, sdeg, tdeg, udeg;
};

#endif
