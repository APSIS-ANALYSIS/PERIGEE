#ifndef NBC_PARTITION_SOLID_HPP
#define NBC_PARTITION_SOLID_HPP
// ============================================================================
// NBC_Partition_Solid.hpp
//
// Solid-specific nodal boundary condition partition with displacement tags.
//
// Date: Mar. 19 2026
// ============================================================================
#include "IPart.hpp"
#include "Map_Node_Index.hpp"
#include "NodalBC_Solid.hpp"

class NBC_Partition_Solid
{
  public:
    NBC_Partition_Solid( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const std::vector<NodalBC_Solid *> &solid_nbc_list_x,
        const std::vector<NodalBC_Solid *> &solid_nbc_list_y,
        const std::vector<NodalBC_Solid *> &solid_nbc_list_z );

    virtual ~NBC_Partition_Solid() = default;

    virtual void write_hdf5( const std::string &FileName ) const
    { write_hdf5(FileName, "/nbc"); }

    virtual void write_hdf5( const std::string &FileName,
        const std::string &GroupName ) const;

  private:
    const int cpu_rank;

    std::vector<int> LID;
    std::vector<int> LDN;
    std::vector<int> LPSN, LPMN;
    std::vector<int> LocalMaster, LocalMasterSlave;
    std::vector<int> Num_LD, Num_LPS, Num_LPM;
    std::vector<int> LDN_is_disp_driven;
};

#endif
