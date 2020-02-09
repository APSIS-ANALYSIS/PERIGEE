#ifndef IEN_NURBS_1PATCH_3D_WITH_POINTER_HPP
#define IEN_NURBS_1PATCH_3D_WITH_POINTER_HPP
// ============================================================================
// NURBS_1Patch_IEN_wPtr.hpp
// Gives the IEN array for 1patch NURBS tensor pruduct mesh
//
// This class returns a pointer to the 1d array of IEN.
// The object is to gives a pointer directly linked to the IEN, which
// is stored as a 1-D array. Hence, IEN can be used directly as eind
// in METIS call.
//
// Date:
// Oct. 10 2013
// ============================================================================

#include "IIEN.hpp"
#include "IMesh.hpp"
#include "Sys_Tools.hpp"

class IEN_NURBS_1Patch_3D_wPtr : public IIEN
{
  public:
    // ------------------------------------------------------------------------ 
    // Constructor for tensor-product NURBS mesh given by the degree,
    // num_of_elem, and num_of_funcs.
    // ------------------------------------------------------------------------ 
    IEN_NURBS_1Patch_3D_wPtr( const IMesh * const &mesh );
   
    // ------------------------------------------------------------------------ 
    // Constructor with zero element removal
    // This constructor will read in the IEN array for the whole mesh and
    // remove zero-measure element from the IEN array quickly using the
    // tensor-product structure of the NURBS mesh. Eventually, the IEN array
    // will only save the info for non-zero elements
    // ------------------------------------------------------------------------ 
    IEN_NURBS_1Patch_3D_wPtr( const IIEN * const &ien_wz,
       const IMesh * const &mesh );

    virtual ~IEN_NURBS_1Patch_3D_wPtr();

    virtual int get_IEN( const int &e, const int &l_node) const;

    virtual void print_IEN() const;

  private:
    int * IEN;
    
    int nElem, nLocBas;
};

#endif
