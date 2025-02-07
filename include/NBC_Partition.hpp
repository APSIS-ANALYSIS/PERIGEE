#ifndef NBC_PARTITION_HPP
#define NBC_PARTITION_HPP
// ============================================================================
// NBC_Partition.hpp
//
// Nodal Boundary Condition Partition implementation for 
// three-dimensional meshes.
//
// Date: March 24 2016
// Author: Ju Liu
// ============================================================================
#include "IPart.hpp"
#include "INodalBC.hpp"
#include "Map_Node_Index.hpp"

class NBC_Partition
{
  public:
    // ------------------------------------------------------------------------
    // Generate partition for a single nodal bc file.
    // The purpose of this one is to generate a nodal file for inflow
    // boundary condition for CFD analysis. I get the inflow bc in the
    // input INodalBC structure, this one will write a partitioned
    // nbc file for analysis use. 
    // ------------------------------------------------------------------------
    NBC_Partition( const IPart * const &part,
       const Map_Node_Index * const &mnindex,
       const std::vector<INodalBC *> &nbc_list );

    virtual ~NBC_Partition() = default;

    // ------------------------------------------------------------------------
    // write_hdf5 : write the nodal bc info into the part file, under
    //              the groupname = /nbc.
    //              This function requires that the part_pxxxxx.h5 file has 
    //              been created.
    // \para FileName : the base name for the partition file (default part)
    // ------------------------------------------------------------------------
    virtual void write_hdf5( const std::string &FileName ) const
    { write_hdf5(FileName, "/nbc"); }

    // ------------------------------------------------------------------------
    // write_hdf5 : write the nodal bc info into the part file, under
    //              the given groupname : GroupName.
    //              This function requires that the part_pxxxxx.h5
    //              file has been created.
    // \para FileName : the base name for the partition file (default part)
    // \para GroupName : the group name
    // ------------------------------------------------------------------------
    virtual void write_hdf5( const std::string &FileName, 
        const std::string &GroupName ) const;

    virtual void print_info() const;

    virtual int get_LID( const int &ii ) const {return LID[ii];}

    virtual int get_LDN( const int &ii ) const {return LDN[ii];}

    virtual int get_LPSN( const int &ii ) const {return LPSN[ii];}

    virtual int get_LPMN( const int &ii ) const {return LPMN[ii];}

    virtual int get_LocalMaster( const int &ii ) const {return LocalMaster[ii];}

    virtual int get_LocalMasterSlave( const int &ii ) const 
    {return LocalMasterSlave[ii];}

    virtual int get_Num_LD( const int &ii ) const {return Num_LD[ii];}

    virtual int get_Num_LPS( const int &ii ) const {return Num_LPS[ii];}

    virtual int get_Num_LPM( const int &ii ) const {return Num_LPM[ii];}

  protected:
    const int cpu_rank;

    // The ID array of the Local nodes.
    std::vector<int> LID;

    // Volumetric nodal indices of the local subdomain's Dirichlet nodes
    // Length is sum Num_LD[ii] for 0 <= ii < dof
    std::vector<int> LDN;

    // Volumetric nodal indices of the local subdomain's slave nodes and 
    // their master nodes
    // Length is sum Num_LPS/LPM[ii] for 0 <= ii < dof
    std::vector<int> LPSN, LPMN;

    // Volumetric nodal indices of the local subdomain's master nodes and 
    // their slave nodes
    // Length is sum Num_LPS/LPM[ii] for 0 <= ii < dof
    std::vector<int> LocalMaster, LocalMasterSlave;

    // Number of Local Dirichlet Nodes for each dof
    // Length is dof
    std::vector<int> Num_LD;

    // Number of Local Periodic Slave nodes
    // Length is dof
    std::vector<int> Num_LPS;

    // Number of Local Periodic Master nodes
    // Length is dof
    std::vector<int> Num_LPM;
};

#endif
