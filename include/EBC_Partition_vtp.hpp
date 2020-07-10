#ifndef EBC_PARTITION_VTP_HPP
#define EBC_PARTITION_VTP_HPP
// ==================================================================
// EBC_Partition_vtp.hpp
//
// Elemental boundary condition partition implementation for data 
// structures readed from .vtp files. This is used mainly for
// unstructured meshes.
//
// Author: Ju Liu
// Date: Jan. 15 2017
// ==================================================================
#include <cassert>
#include "IEBC_Partition.hpp"
#include "IPart.hpp"
#include "ElemBC.hpp"

class EBC_Partition_vtp : public IEBC_Partition
{
  public:
    EBC_Partition_vtp( const IPart * const &part, 
        const Map_Node_Index * const &mnindex,
        const ElemBC * const &ebc );

    virtual ~EBC_Partition_vtp();

    // Write the data to hdf5 file in group /GroupName
    virtual void write_hdf5( const char * const &FileName,
       const char * const &GroupName = "/ebc" ) const;

    virtual void print_info() const;

    virtual int get_num_ebc() const {return num_ebc;}

    virtual int get_num_local_node(const int &ii) const
    {return num_local_node[ii];}

    virtual int get_num_local_cell(const int &ii) const
    {return num_local_cell[ii];}

    virtual int get_cell_nLocBas(const int &ii) const
    {return cell_nLocBas[ii];}

    virtual double get_local_pt_xyz(const int &ii, const int &jj) const
    {return local_pt_xyz[ii][jj];}

    virtual int get_local_tri_ien(const int &ii, const int &jj) const
    {return local_tri_ien[ii][jj];}

    virtual int get_local_global_node(const int &ii, const int &jj) const
    {return local_global_node[ii][jj];}

    virtual int get_local_node_pos(const int &ii, const int &jj) const
    {return local_node_pos[ii][jj];}

    virtual int get_local_global_cell(const int &ii, const int &jj) const
    {return local_global_cell[ii][jj];}

  protected:
    // The rank or id of the subdomain
    const int cpu_rank; 
    
    // The number of groups of bc faces that require bc integral
    const int num_ebc; 

    std::vector<int> num_local_node, num_local_cell, cell_nLocBas;

    // local node's coordinates, num_ebc x (3 x num_local_node[ii]) in size
    std::vector< std::vector<double> > local_pt_xyz;

    // local cell's IEN array, 
    // num_ebc x (cell_nLocBas[ii] x num_local_cell[ii]) in size
    std::vector< std::vector<int> > local_tri_ien;

    // local node's global index, num_ebc x num_local_node[ii] in size
    std::vector< std::vector<int> > local_global_node;
    
    // local node's position in the local_to_global array
    std::vector< std::vector<int> > local_node_pos;

    // local cell's global index num_ebc x num_local_cell[ii] in size
    std::vector< std::vector<int> > local_global_cell;
};

#endif
