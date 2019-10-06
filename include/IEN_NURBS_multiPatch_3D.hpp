#ifndef IEN_NURBS_MULTIPATCH_3D_HPP
#define IEN_NURBS_MULTIPATCH_3D_HPP
// ============================================================================
// IEN_NURBS_multiPatch_3D.hpp
// 
// This object gives the IEN array for multi-patch NURBS 3D geometries.
//
// Date: Aug. 31 2015
// ============================================================================

#include "IIEN.hpp"
#include "IMesh.hpp"
#include "IEN_NURBS_1Patch_3D_wPtr.hpp"

class IEN_NURBS_multiPatch_3D : public IIEN
{
  public:
    IEN_NURBS_multiPatch_3D( const IMesh * const &mmesh );
    
    virtual ~IEN_NURBS_multiPatch_3D();

    virtual s_int get_IEN( const s_int &ee, const s_int &l_node ) const
    {return IEN[ee*nLocBas + l_node];}

    virtual void print_IEN() const;

  private:
    s_int * IEN;
    const s_int nElem;
    const int nLocBas;
};

#endif
