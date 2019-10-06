#ifndef NBC_PARTITION_2D_HPP
#define NBC_PARTITION_2D_HPP
// ============================================================================
// NBC_Partition_2D.hpp
//
// Nodal Boundary Condition Partition implementation for two-dimensional
// meshes.
//
// Date: Aug. 19 2015
// ============================================================================
#include "INBC_Partition.hpp"
#include "IPart.hpp"
#include "INodalBC.hpp"

class NBC_Partition_2D : public INBC_Partition
{
  public:
    NBC_Partition_2D( const IPart * const &part,
       const Map_Node_Index * const &mnindex,
       const std::vector<INodalBC *> &nbc_list );

    virtual ~NBC_Partition_2D();

    virtual void write_hdf5(const char * FileName) const;

    virtual void print_info() const;

  private:
    const int cpu_rank;

    // The ID array of the Local nodes.
    // length : nlocghonode * dof
    std::vector<int> LID;

    // Local subdomain's Dirichlet nodes
    // length : Num_LD[0] + ... + Num_LD[dof-1]
    std::vector<int> LDN;

    // Local subdomain's slave nodes and their master nodes
    // length : Num_LPS[0] + ... + Num_LPS[dof-1]
    std::vector<int> LPSN, LPMN;

    // Local subdomain's master nodes and their slaves
    // length : Num_LPM[0] + ... + Num_LPM[dof-1]
    std::vector<int> LocalMaster, LocalMasterSlave;

    // Number of Local Dirichlet Nodes for each dof
    // length : dof
    std::vector<int> Num_LD;

    // Number of Local Periodic Slave nodes
    // length : dof
    std::vector<int> Num_LPS;

    // Number of Local Periodic Master nodes
    // length : dof
    std::vector<int> Num_LPM;
};

#endif
