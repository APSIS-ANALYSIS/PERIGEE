#ifndef EBC_PARTITION_FEM_HPP
#define EBC_PARTITION_FEM_HPP
// ==================================================================
// EBC_Partition_FEM.hpp
//
// Elemental boundary condition partition implementation mainly for
// unstructured meshes.
// Input: 
// IPart : the local partition information
// Map_Node_Index : Nodal index mapping between old and new indices
// ElemBC : the elemental BC
//
// Author: Ju Liu
// Date: Nov. 22 2017
// ==================================================================
#include <cassert>
#include "IEBC_Partition.hpp"
#include "IPart.hpp"
#include "ElemBC.hpp"

class EBC_Partition_FEM : public IEBC_Partition
{
  public:
    EBC_Partition_FEM( const IPart * const &part, 
        const Map_Node_Index * const &mnindex,
        const ElemBC * const &ebc );

    virtual ~EBC_Partition_FEM();

    virtual void write_hdf5( const char * FileName ) const;
    
    virtual void write_hdf5( const char * FileName,
       const char * GroupName ) const;

    virtual void print_info() const;

    virtual int get_num_ebc() const {return num_ebc;}

    virtual int get_num_local_cell_node(const int &ii) const
    {return num_local_cell_node[ii];}

    virtual int get_num_local_cell(const int &ii) const
    {return num_local_cell[ii];}

    virtual int get_cell_nLocBas(const int &ii) const
    {return cell_nLocBas[ii];}

    virtual double get_local_cell_node_xyz(const int &ii, const int &jj) const
    {return local_cell_node_xyz[ii][jj];}

    virtual int get_local_tri_ien(const int &ii, const int &jj) const
    {return local_tri_ien[ii][jj];}

    virtual int get_local_cell_node_vol_id(const int &ii, const int &jj) const
    {return local_cell_node_vol_id[ii][jj];}

    virtual int get_local_cell_node_pos(const int &ii, const int &jj) const
    {return local_cell_node_pos[ii][jj];}

    virtual int get_local_cell_vol_id(const int &ii, const int &jj) const
    {return local_cell_vol_id[ii][jj];}

  protected:
    const int cpu_rank;
    const int num_ebc;

    std::vector<int> num_local_cell_node, num_local_cell, cell_nLocBas;

    // local cell node's coordinates, num_ebc x (3 x num_local_cell_node[ii]) in size
    std::vector< std::vector<double> > local_cell_node_xyz;

    // local cell's IEN array, 
    // num_ebc x (cell_nLocBas[ii] x num_local_cell[ii]) in size
    std::vector< std::vector<int> > local_tri_ien;

    // local cell node's global index, num_ebc x num_local_cell_node[ii] in size
    std::vector< std::vector<int> > local_cell_node_vol_id;
    
    // local cell node's position in the local_to_global array
    std::vector< std::vector<int> > local_cell_node_pos;

    // local cell's global index num_ebc x num_local_cell[ii] in size
    std::vector< std::vector<int> > local_cell_vol_id;

    // local cell's interior point's coordinates
    // num_ebc x { 3 num_local_cell[ii] }
    std::vector< std::vector<double> > local_intpt;
};

#endif
