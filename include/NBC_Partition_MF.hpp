#ifndef NBC_PARTITION_MF_HPP
#define NBC_PARTITION_MF_HPP
// ============================================================================
// NBC_Partition_MF.hpp
//
// Nodal boundary condition implementation for Multi-Field problems.
// 
// For Multi-Field problems, we need to document the assigned row (or column)
// index for the (node, dof) pair.
//
// Date: Dec. 17 2021
// Author: Ju Liu
// ============================================================================
#include "NBC_Partition.hpp"

class NBC_Partition_MF : public NBC_Partition
{
  public:
    NBC_Partition_MF( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const std::vector<INodalBC *> &nbc_list,
        const std::vector< std::vector<int> > &grid2id );

    // If the grid2id mapper is not provided for the constructor, we assume that
    // the row/col index is given by the dof x node_new_numbering + mm,
    // for 0 <= mm < dof
    NBC_Partition_MF( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const std::vector<INodalBC *> &nbc_list );

    virtual ~NBC_Partition_MF();

    virtual void write_hdf5( const std::string &FileName ) const
    { write_hdf5(FileName, "/nbc"); }

    virtual void write_hdf5( const std::string &FileName,
        const std::string &GroupName ) const;

  protected:
    // LID mapped to the MF value, which means the actual row/col index in the
    // matrix problem.
    std::vector<int> LID_MF;

    // LDN, LPSN, LPMN, LocalMaster, LocalMasterSlave mapped to MF value
    // The mapping is from the grid2id mapper, that maps the (node, dof) value
    // to the actual matrix problem's row/col location.
    std::vector<int> LDN_MF, LPSN_MF, LPMN_MF, LocalMaster_MF, LocalMasterSlave_MF;
};

#endif
