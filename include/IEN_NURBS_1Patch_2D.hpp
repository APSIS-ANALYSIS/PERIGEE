#ifndef IEN_NURBS_1PATCH_2D_HPP
#define IEN_NURBS_1PATCH_2D_HPP
// ========================================================
// IEN_NURBS_1Patch_2D.hpp
// Gives the IEN array for single patch tensor product mesh.
//
// Date: Jan 26 2014
// ========================================================

#include "IIEN.hpp"
#include "IMesh.hpp"

class IEN_NURBS_1Patch_2D : public IIEN
{
  public:
    IEN_NURBS_1Patch_2D( const IMesh * const &mesh );

    virtual ~IEN_NURBS_1Patch_2D();

    virtual s_int get_IEN( const s_int &e, const s_int &l_node) const
    {return IEN[e][l_node];}

    virtual void print_IEN() const;

  private:
    s_int ** IEN;
    s_int nElem, nLocBas;
};

#endif
