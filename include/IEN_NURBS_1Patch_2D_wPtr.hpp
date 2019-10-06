#ifndef IEN_NURBS_1PATCH_2D_WPTR_HPP
#define IEN_NURBS_1PATCH_2D_WPTR_HPP
// ============================================================================
// IEN_NURBS_1Patch_2D_wPtr.hpp
// IEN class for single patch NURBS 2D geometry. This class will store the IEN
// array in a one-dimensional array and give a pointer to that array for
// direct data accessing. This design is intended to reduce the preparation for
// the call of the METIS library.
//
// Date: Oct. 5th 2015
// ============================================================================
#include "IIEN.hpp"
#include "IMesh.hpp"
#include "Sys_Tools.hpp"

class IEN_NURBS_1Patch_2D_wPtr : public IIEN
{
  public:
    IEN_NURBS_1Patch_2D_wPtr( const IMesh * const &mesh );


    // Constructor:
    // It will read in the IEN for the whole mesh and extracts the element with
    // nonzero measure to create its own IEN, which is only for the elements
    // with nonzero measures.
    IEN_NURBS_1Patch_2D_wPtr( IIEN const * const &ien_wz,
        IMesh const * const &mesh );


    virtual ~IEN_NURBS_1Patch_2D_wPtr();

    virtual int get_IEN( const int &e, const int &node ) const;

    virtual void print_IEN() const;

    virtual int * get_IENptr() const;

  private:
    int * IEN;
    int nElem, nLocBas;
};


#endif
