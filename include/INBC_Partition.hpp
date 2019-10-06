#ifndef INBC_PARTITION_HPP
#define INBC_PARTITION_HPP
// ============================================================================
// INBC_Partition.hpp
//
// Interface for Nodal Boundary Condition Partition header file.
//
// This is a pure virtual class that provides an interface for nodal boundary
// condition partition. The partitioned nodal info will be cached as a HDF5
// file.
//
// Date: Aug 19 2015
// ============================================================================
#include "Sys_Tools.hpp"

class INBC_Partition
{
  public:
    INBC_Partition(){};

    virtual ~INBC_Partition(){};

    virtual void write_hdf5( const char * FileName ) const = 0;
    
    virtual void write_hdf5( const char * FileName,
       const char * GroupName ) const
    {SYS_T::print_exit("Error: write_hdf5(fname, gname) is not implemented. \n");}

    virtual void write_inflow_hdf5( const char * FileName,
       const std::vector<double> &outVec,
       const double &act_inflow_area ) const
    {SYS_T::print_exit("Error: write_inflow_hdf5 is not implemented. \n");}

    virtual void print_info() const = 0;

    // virtual function to access the private data
    virtual int get_LID( const int &pos ) const
    {SYS_T::print_exit("Error: get_LID is not implemented. \n"); return 0;}

    virtual int get_LDN( const int &pos ) const
    {SYS_T::print_exit("Error: get_LDN is not implemented. \n"); return 0;}

    virtual int get_LPSN( const int &pos ) const
    {SYS_T::print_exit("Error: get_LPSN is not implemented. \n"); return 0;}

    virtual int get_LPMN( const int &pos ) const
    {SYS_T::print_exit("Error: get_LPMN is not implemented. \n"); return 0;}

    virtual int get_LocalMaster( const int &pos ) const
    {SYS_T::print_exit("Error: get_LocalMaster is not implemented. \n"); return 0;}

    virtual int get_LocalMasterSlave( const int &pos ) const
    {SYS_T::print_exit("Error: get_LocalMasterSlave is not implemented. \n"); return 0;}

    virtual int get_Num_LD( const int &pos ) const
    {SYS_T::print_exit("Error: get_Num_LD is not implemented. \n"); return 0;}

    virtual int get_Num_LPS( const int &pos ) const
    {SYS_T::print_exit("Error: get_Num_LPS is not implemented. \n"); return 0;}

    virtual int get_Num_LPM( const int &pos ) const
    {SYS_T::print_exit("Error: get_Num_LPM is not implemented. \n"); return 0;}
};


#endif
