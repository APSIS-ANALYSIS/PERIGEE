#ifndef BOUNDARYCOND_TEST_HPP
#define BOUNDARYCOND_TEST_HPP
// ============================================================================
// BoundaryCond_Test.hpp
// Object:
// Perform test for the boundary condition specification.
//
// Date: Sept. 19 2015
// ============================================================================
#include <cstdlib>
#include "INodalBC.hpp"
#include "BoundaryCond2D.hpp"
#include "BoundaryCond.hpp"
#include "ElemBC.hpp"

namespace TEST_T
{
  // -------------------------------------------------------------------------- 
  // This routine will check the control points for periodic master and slave
  // nodes.
  // -------------------------------------------------------------------------- 
  void BC_Periodic_coordinate_check( const std::vector<double> &cPts,
      const BoundaryCond * const &bc, const double &tol );

  // -------------------------------------------------------------------------- 
  // Correct the control points for master-slave pairs
  //   Due to the error in refinement algorithms, the control points for
  //   master-slave pairs in the refined mesh may be non-matching. This routine
  //   is used to enforce the control point of the slave node to be identical to
  //   that of the master node.
  // -------------------------------------------------------------------------- 
  void correct_masterslave_ctrlPts( std::vector<double> &cPts,
      const BoundaryCond * const &bc );

  // -------------------------------------------------------------------------- 
  // This routine will compare two BC objects and make sure they give identical
  // BC specification
  // -------------------------------------------------------------------------- 
  void NBC_compare( BoundaryCond2D const * const &bc,
     INodalBC const * const &nbc );

  // -------------------------------------------------------------------------- 
  // This routine will check the ElemBC's node index and coordinates are
  // compatible with the given global volumetric node index
  // -------------------------------------------------------------------------- 
  void EBC_node_compatibility_check( const ElemBC * const &bc, 
      std::vector<double> &ctrlPts );

  // -------------------------------------------------------------------------- 
  // This routine will check that the ElemBC's cell IEN matches a volumetric
  // cell's face.
  // -------------------------------------------------------------------------- 
  void EBC_cell_IEN_check( const ElemBC * const &bc,
      const IIEN * const &ien );
}

#endif
