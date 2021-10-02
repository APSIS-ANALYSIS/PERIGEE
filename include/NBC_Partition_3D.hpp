#ifndef NBC_PARTITION_3D_HPP
#define NBC_PARTITION_3D_HPP
// ==================================================================
// NBC_Partition_3D.hpp
//
// Nodal Boundary Condition Partition implementation for 
// three-dimensional meshes.
//
// Date: March 24 2016
// Author: Ju Liu
// ==================================================================
#include "INBC_Partition.hpp"
#include "IPart.hpp"
#include "Map_Node_Index.hpp"
#include "INodalBC.hpp"
#include "HDF5_Writer.hpp"

class NBC_Partition_3D : public INBC_Partition
{
  public:
    NBC_Partition_3D( const IPart * const &part,
       const Map_Node_Index * const &mnindex,
       const std::vector<INodalBC *> &nbc_list );

    // --------------------------------------------------------------
    // Generate partition for a single nodal bc file.
    // The purpose of this one is to generate a nodal file for inflow
    // boundary condition for CFD analysis. I get the inflow bc in the
    // input INodalBC structure, this one will write a partitioned
    // nbc file for analysis use. 
    // --------------------------------------------------------------
    NBC_Partition_3D( const IPart * const &part,
       const Map_Node_Index * const &mnindex,
       const INodalBC * const &nbc );

    virtual ~NBC_Partition_3D();

    // --------------------------------------------------------------
    // write_hdf5 : write the nodal bc info into the part file.
    //              this function requires that the part_pxxxxx.h5
    //              file has been created.
    // \para FileName : the base name for the partition file (default part)
    // --------------------------------------------------------------
    virtual void write_hdf5(const char * FileName) const;

    // -------------------------------------------------------------- 
    // write_hdf5 : write the nodal bc info into the part file, under
    //              the given groupname : GroupName.
    //              This function requires that the part_pxxxxx.h5
    //              file has been created.
    // \para FileName : the base name for the partition file (default part)
    // \para GroupName : the group name
    // -------------------------------------------------------------- 
    virtual void write_hdf5(const char * FileName, 
        const char * GroupName) const;

    virtual void print_info() const;

    virtual int get_LID( const int &ii ) const {return LID[ii];}

    virtual int get_LDN( const int &nbc_id, const int &ii ) const
    {return LDN[nbc_id][ii];}

    virtual int get_LPSN( const int &nbc_id, const int &ii ) const
    {return LPSN[nbc_id][ii];}

    virtual int get_LPMN( const int &nbc_id, const int &ii ) const
    {return LPMN[nbc_id][ii];}

    virtual int get_LocalMaster( const int &nbc_id, const int &ii ) const
    {return LocalMaster[nbc_id][ii];}

    virtual int get_LocalMasterSlave( const int &nbc_id, const int &ii ) const 
    {return LocalMasterSlave[nbc_id][ii];}

    virtual int get_Num_LD( const int &nbc_id, const int &ii ) const
    {return Num_LD[nbc_id][ii];}

    virtual int get_Num_LPS( const int &nbc_id, const int &ii ) const
    {return Num_LPS[nbc_id][ii];}

    virtual int get_Num_LPM( const int &nbc_id, const int &ii ) const
    {return Num_LPM[nbc_id][ii];}

  protected:
    const int cpu_rank;

    const int num_nbc;

    // The ID array of the Local nodes.
    std::vector<int> LID;

    // Local subdomain's Dirichlet nodes
    // num_nbc times ( dof x Num_LD[ii][jj] )
    std::vector< std::vector<int> > LDN;

    // Local subdomain's slave nodes and their master nodes
    // num_nbc times ( dof x Num_LPS/LPM[ii][jj] )
    std::vector< std::vector<int> > LPSN, LPMN;

    // Local subdomain's master nodes and their slaves
    // num_nbc times ( dof x Num_LPS/LPM[ii][jj] )
    std::vector< std::vector<int> > LocalMaster, LocalMasterSlave;

    // Number of Local Dirichlet Nodes for each dof
    // num_nbc x dof
    std::vector< std::vector<int> > Num_LD;

    // Number of Local Periodic Slave nodes
    // num_nbc x dof
    std::vector< std::vector<int> > Num_LPS;

    // Number of Local Periodic Master nodes
    // num_nbc x dof
    std::vector< std::vector<int> > Num_LPM;
};

#endif
