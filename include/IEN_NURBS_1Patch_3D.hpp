#ifndef IEN_NURBS_1PATCH_3D_HPP
#define IEN_NURBS_1PATCH_3D_HPP
// ==================================================================
// IEN_NURBS_1Patch_3D.hpp
// 
// This object gives the IEN array for 1patch NURBS tensor product mesh.
//
// Date: Sept. 24 2013
// ==================================================================
#include "Sys_Tools.hpp"
#include "IIEN.hpp"
#include "IMesh.hpp"

class IEN_NURBS_1Patch_3D : public IIEN
{
  public:
    IEN_NURBS_1Patch_3D( const IMesh * const &mesh );
    virtual ~IEN_NURBS_1Patch_3D();

    virtual int get_IEN( const int &e, const int &l_node) const 
    {return IEN[e][l_node];}

    virtual void print_IEN() const;

  private:
    int ** IEN;
    int nElem, nLocBas;
};

#endif
