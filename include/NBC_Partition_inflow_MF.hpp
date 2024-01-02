#ifndef NBC_PARTITION_INFLOW_MF_HPP
#define NBC_PARTITION_INFLOW_MF_HPP
// ============================================================================
// NBC_Partition_inflow_MF.hpp
//
// Record the corresponding matrix location for the dof on the nodes of inflow
// type.
//
// Author: Ju Liu
// Date: Dec. 17 2021
// ============================================================================
#include "NBC_Partition_inflow.hpp"

class NBC_Partition_inflow_MF : public NBC_Partition_inflow
{
  public:
    NBC_Partition_inflow_MF( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const INodalBC * const &nbc,
        const std::vector< std::vector<int> > &grid2id );

    virtual ~NBC_Partition_inflow_MF() = default;

    virtual void write_hdf5( const std::string &FileName ) const;

  protected:
    // size is num_nbc x (3 x Num_LD[ii]), 0 <= ii < num_nbc
    std::vector< std::vector<int> > LDN_MF;
};

#endif
