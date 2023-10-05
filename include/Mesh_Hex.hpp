#ifndef MESH_HEX_HPP
#define MESH_HEX_HPP

#include "IMesh.hpp"

class Mesh_Hex : public IMesh
{
  public:
    Mesh_Hex(const int &in_nFunc, const int &in_nElem, const int &in_deg);

    virtual ~Mesh_Hex();

    virtual void print_info() const;

    virtual int get_degree() const {return deg;}

    virtual int get_nFunc() const {return nFunc;}

    virtual int get_nElem() const {return nElem;}

    virtual int get_nLocBas() const {return nLocBas;}

  private:
    const int nFunc, nElem, deg;
   
    // ------------------------------------------------------------------------ 
    // In this class, nLocBas is determined by the element type, and since we
    // are restricted to hex element, the degree determines the value of it.
    // ------------------------------------------------------------------------ 
    int nLocBas;
};

#endif