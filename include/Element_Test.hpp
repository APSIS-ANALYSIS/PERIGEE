#ifndef ELEMENT_TEST_HPP
#define ELEMENT_TEST_HPP
// ============================================================================
// Element_Test.hpp
// Object:
// Perform various test for the element routine to verify the correctness.
//
// Date: June 27 2015
// ============================================================================
#include "Math_Tools.hpp"
#include "IAGlobal_Mesh_Info.hpp"
#include "ALocal_Elem.hpp"
#include "FEAElement.hpp"

namespace TEST_T
{
  // Calculate the Jacobian matrix times the inverse of the Jacobian matrix
  // and make sure that it is identity matrix.
  void Element_JacxInvJac_check( const std::vector<FEAElement *> &earray,
      const ALocal_Elem * const &locelem,
      const IALocal_meshSize * const &msize );

  // Calculate the Jacobian and inverse of Jacobian and expect the result is
  // close to the Identity matrix with tol = 1.0e-15
  void Element_JacxInvJac_check( FEAElement const * const &ele1 );

  // Compare two element routine gives the same ds_dx matrix
  void Element_invJac_compare(const std::vector<FEAElement *> &earray1,
      const std::vector<FEAElement *> &earray2,
      const ALocal_Elem * const &locelem,
      const IALocal_meshSize * const &msize );

  void Element_invJac_compare( FEAElement const * const &ele1,
      FEAElement const * const &ele2,
      const double &ehx, const double &ehy, const double &ehz,
      const double &tol );

  // compare dx_ds
  void Element_Jac_compare(  FEAElement const * const &ele1,
      FEAElement const * const &ele2, const double &tol );

  // compare two element routien gives the same R, dR_dx, dR_dy, dR_dz
  // d2R_dxx, d2R_dyy, d2R_dzz.
  void Element_R_dR_LapR_compare(
      const std::vector<FEAElement *> &earray1, 
      const std::vector<FEAElement *> &earray2,
      const IAGlobal_Mesh_Info * const &gmiptr,
      const ALocal_Elem * const &locelem,
      const IALocal_meshSize * const &msize );

  // compare two element's get_2D_R_dR_d2R function as well as the 
  // get_detJac function
  void Element_2D_R_dR_d2R_compare(
      FEAElement const * const &ele1, FEAElement const * const &ele2,
      const double &tol );

  // compare two element's R
  void Element_R_compare( FEAElement const * const &ele1,
      FEAElement const * const &ele2, const double &tol ); 

  void Element_grad2d_compare(FEAElement const * const &ele1,
      FEAElement const * const &ele2, const double &tol );

  void Element_grad3d_compare(FEAElement const * const &ele1,
      FEAElement const * const &ele2, const double &tol );

  void Element_d2R_compare( FEAElement const * const &nurbs,
      FEAElement const * const &bs, const double &tol );

  void Element_d2R_3d_compare( FEAElement const * const &nurbs,
      FEAElement const * const &bs, const double &tol );

  void Element_normal_vector_compare( FEAElement const * const &nurbs,
      FEAElement const * const &bs, const double &tol );

  void Element_normal_vector_3d_compare( FEAElement const * const &nurbs,
      FEAElement const * const &bs, const double &tol );
}

#endif
