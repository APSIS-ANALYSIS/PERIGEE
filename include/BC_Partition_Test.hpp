#ifndef BC_PARTITION_TEST
#define BC_PARTITION_TEST
// ==================================================================
// BC_Partition_Test.hpp
// Object: Test the BC_Partition class
//
// Date: Oct. 13 2013
// ==================================================================
#include <cassert>
#include "IPart.hpp"
#include "BoundaryCond.hpp"
#include "Map_Node_Index.hpp"
#include "BC_Partition.hpp"

void BC_Partition_Node_LID_Test( const BC_Partition * const &bcpart,
   const IPart * const &part, const int dof );

void BC_Partition_LID_Node_Test( const BC_Partition * const &bcpart,
    const std::vector<BoundaryCond *> &bc_list, const IPart * const &part,
    const Map_Node_Index * const &mnindex, const int dof );

void BC_Partition_Node_LID_Test_2( const BC_Partition * const &bcpart,
    const std::vector<BoundaryCond *> &bc_list, const IPart * const &part,
    const Map_Node_Index * const &mnindex, const int dof );


bool isIndexInLDN( const s_int index, const int dof, const BC_Partition * const &bcpart );
bool isIndexInLPSN( const s_int index, const int dof, const BC_Partition * const &bcpart );
bool isIndexInLPMN( const s_int index, const int dof, const BC_Partition * const &bcpart );


bool isIndexInDN( const s_int index, const int dof, const BoundaryCond * const &bc,
   const Map_Node_Index * const &mnindex );
bool isIndexInPSN( const s_int index, const int dof, const BoundaryCond * const &bc,
   const Map_Node_Index * const &mnindex, int &pos );
bool isIndexInPMN( const s_int index, const int dof, const BoundaryCond * const &bc,
   const Map_Node_Index * const &mnindex );

#endif
