#ifndef IEN_TEST_HPP
#define IEN_TEST_HPP
// ========================================================
// IEN_Test.hpp
// Test and verify my IEN routine.
//
// Date: Sept. 24th 2013
// ========================================================
#include <cassert>
#include "IMesh.hpp"
#include "IEN_NURBS_1Patch_3D.hpp"

namespace TEST_T
{
  void IEN_NURBS_1Patch_3D_Test( const IIEN * const &ien,
      const IMesh * const &mesh );
  
  void IEN_NURBS_multiPatch_3D_Test( const IMesh * const mmesh, 
     const IIEN * const &in_ien );

  // Compare the contents of two IEN classes to see if they are identical
  void IEN_NURBS_compare( IIEN const * const &ien1, IIEN const * const &ien2,
     const int &nElem, const int &nLocBas );

  // Compare to make sure the nonzero IEN is compatible with the original IEN
  void IEN_NURBS_2D_nz_check( IIEN const * const &ien_wz, 
      IIEN const * const &ien_nz,
      IMesh const * const &mesh );
}

#endif
