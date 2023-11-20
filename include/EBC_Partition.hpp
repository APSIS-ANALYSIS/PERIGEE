#ifndef EBC_PARTITION_HPP
#define EBC_PARTITION_HPP
// ==================================================================
// EBC_Partition.hpp
//
// Elemental boundary condition partition implementation for data 
// structures readed from .vtp files. This is used mainly for
// unstructured meshes.
//
// Author: Ju Liu
// Date: Jan. 15 2017
// ==================================================================
#include "IPart.hpp"
#include "Map_Node_Index.hpp"
#include "ElemBC.hpp"
#include "HDF5_Writer.hpp"

class EBC_Partition
{
  public:
    EBC_Partition( const IPart * const &part, 
        const Map_Node_Index * const &mnindex,
        const ElemBC * const &ebc );

    virtual ~EBC_Partition() = default;

    // Write the data to hdf5 file in group /GroupName
    virtual void write_hdf5( const std::string &FileName,
       const std::string &GroupName ) const;

    virtual void write_hdf5( const std::string &FileName ) const
    { write_hdf5( FileName, "/ebc" ); }

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

    virtual int get_local_cell_ien(const int &ii, const int &jj) const
    {return local_cell_ien[ii][jj];}

    virtual int get_local_cell_node_vol_id(const int &ii, const int &jj) const
    {return local_cell_node_vol_id[ii][jj];}

    virtual int get_cell_local_node_pos(const int &ii, const int &jj) const
    {return local_cell_node_pos[ii][jj];}

    virtual int get_local_cell_vol_id(const int &ii, const int &jj) const
    {return local_cell_vol_id[ii][jj];}

  protected:
    const int cpu_rank; // The rank or id of the subdomain
    
    const int num_ebc; // Number of groups of bc faces that require bc integral

    // length is num_ebc
    std::vector<int> num_local_cell_node, num_local_cell, cell_nLocBas;

    // local cell node's coordinates, num_ebc x (3 x num_local_cell_node[ii]) in size
    std::vector< std::vector<double> > local_cell_node_xyz;

    // local cell's IEN array, 
    // num_ebc x (cell_nLocBas[ii] x num_local_cell[ii]) in size
    std::vector< std::vector<int> > local_cell_ien;

    // local cell node's global index, num_ebc x num_local_cell_node[ii] in size
    // local means local to the CPU's subdomain
    // global means indices in the volumetric mesh
    std::vector< std::vector<int> > local_cell_node_vol_id;
   
    // local cell node's index in the surface wall mesh
    // here, the local cell node is defined as all nodes of the local elements,
    // which means there could be some ghost nodes for the mesh partitioning
    // num_ebc x num_local_node[ii]
    // note: this data is NOT stored into the hdf5 file
    std::vector< std::vector<int> > local_cell_node;

    // local cell node's position in the local_to_global array
    // note: local_to_global array is generated in the Part_xxx based on IPart
    // class, which stores the nodal indices of the local nodes followed by the
    // ghost nodes
    std::vector< std::vector<int> > local_cell_node_pos;

    // local cell's global index num_ebc x num_local_cell[ii] in size
    std::vector< std::vector<int> > local_cell_vol_id;
};

#endif
