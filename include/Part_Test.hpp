#ifndef PART_TEST_HPP
#define PART_TEST_HPP
// ==================================================================
// Part_Test.hpp
//
// Object: Perform various logical test for the mesh partition.
//
// Date: Jan. 12 2017
// ==================================================================
#include <cassert>
#include "IPart.hpp"
#include "IIEN.hpp"
#include "IMesh.hpp"
#include "Map_Node_Index.hpp"
#include "INBC_Partition.hpp"
#include "INodalBC.hpp"
#include "EBC_Partition_vtp.hpp"
#include "ALocal_EBC.hpp"

using std::cout;
using std::cerr;
using std::endl;

namespace TEST_T
{
  // Check the LIEN compatibility with the IEN
  void Part_LIEN_Test(
      const IPart * const &part,
      const Map_Node_Index * const &mnindex,
      const IIEN * const &IEN );


  // Check the node numbering compatibility with global node index
  void Part_Node_Test(
      const IPart * const &part,
      const Map_Node_Index * const &mnindex);

  
  // Check the geometry compatibility
  void Part_CtrlPts_Test(
      const IPart * const &part,
      const Map_Node_Index * const &mnindex,
      const std::vector<double> &ctrlPts );


  // Check the compatibility of LID with the local dirichlet and 
  // periodic nodes.
  void Part_NBC_Test( const IPart * const &part,
      const Map_Node_Index * const &mnindex,
      const std::vector<INodalBC *> &nbc,
      const INBC_Partition * const &nbcpart,
      const int &dof );


  // Check the ALocal_EBC, verify that it is identical to the 
  // IEBC_Partition object.
  void EBCPart_AEBC_Test(
      const IEBC_Partition * const &ebc,
      const ALocal_EBC * const &ebc_part );


  // Check the compatibility of the EBC_part with the mesh partition
  // and the original global ebc.
  void Part_EBC_Test( const IEBC_Partition * const &ep,
      const IPart * const &part,
      const ElemBC * const &ebc );

}

#endif
