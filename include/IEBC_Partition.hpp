#ifndef IEBC_PARTITION_HPP
#define IEBC_PARTITION_HPP
// ============================================================================
// IEBC_Partition.hpp
//
// Interface for Elemental Boundary Condition Partition Header file.
//
// This is a pure virtual class that provides an interface for elemental
// boundary condition partition. The partitioned elemental info is saved as a
// HDF5 file.
//
// Author: Ju Liu
// Date: Aug. 19 2015
// ============================================================================
#include "Sys_Tools.hpp"

class IEBC_Partition
{
  public:
    IEBC_Partition(){};

    virtual ~IEBC_Partition(){};

    // Write the EBC partitioned info in the root directory of a HDF5 file
    virtual void write_hdf5( const char * FileName ) const = 0;

    // Write the EBC partitioned info in a prescribed directory of a HDF5 file
    virtual void write_hdf5( const char * FileName,
        const char * GroupName ) const
    {SYS_T::print_exit("Error: write_hdf5(fname,gname) is not implemented.\n");}

    virtual void print_info() const = 0;

    virtual int get_num_ebc() const 
    {SYS_T::print_exit("Error: IEBC_Partition::get_num_ebc is not implemented. \n"); return 0;}

    virtual int get_num_local_node(const int &ii) const
    {SYS_T::print_exit("Error: IEBC_Partition::get_num_local_node is not implemented. \n"); return 0;}

    virtual int get_num_local_cell(const int &ii) const
    {SYS_T::print_exit("Error: IEBC_Partition::get_num_local_cell is not implemented. \n"); return 0;}

    virtual int get_cell_nLocBas(const int &ii) const
    {SYS_T::print_exit("Error: IEBC_Partition::get_cell_nLocBas is not implemented. \n"); return 0;}

    virtual double get_local_pt_xyz(const int &ii, const int &jj) const
    {SYS_T::print_exit("Error: IEBC_Partition::get_local_pt_xyz is not implemented. \n"); return 0.0;}

    virtual int get_local_tri_ien(const int &ii, const int &jj) const
    {SYS_T::print_exit("Error: IEBC_Partition::get_local_tri_ien is not implemented. \n"); return 0;}

    virtual int get_local_global_node(const int &ii, const int &jj) const
    {SYS_T::print_exit("Error: IEBC_Partition::get_local_global_node is not implemented. \n"); return 0;}

    virtual int get_local_node_pos(const int &ii, const int &jj) const
    {SYS_T::print_exit("Error: IEBC_Partition::get_local_node_pos is not implemented. \n"); return 0;}

    virtual int get_local_global_cell(const int &ii, const int &jj) const
    {SYS_T::print_exit("Error: IEBC_Partition::get_local_global_cell is not implemented. \n"); return 0;}

};


#endif
