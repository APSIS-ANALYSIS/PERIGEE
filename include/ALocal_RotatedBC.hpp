#ifndef ALOCAL_ROTATEDBC_HPP
#define ALOCAL_ROTATEDBC_HPP
// ============================================================================
// ALocal_RotatedBC.hpp
//
// Analysis-use data that is related to the rotated surface nodes.
//
// Date Created: Sep 11 2024
// Author: Yujie Sun
// ============================================================================
#include "HDF5_Reader.hpp"
#include "Math_Tools.hpp"

class ALocal_RotatedBC
{
  public:
    ALocal_RotatedBC( const std::string &fileBaseName, const int &cpu_rank );

    virtual ~ALocal_RotatedBC() = default;

    // ------------------------------------------------------------------------
    // Get global index of a Dirichlet node in the local partition
    // 0 <= node < Num_LD
    // Note: make sure that Num_LD > 0 before calling this get function
    // ------------------------------------------------------------------------
    virtual int get_LDN( const int &node ) const
    {return LDN[node];}

    // ------------------------------------------------------------------------
    // Get the number of Dirichlet nodes in the local partition
    // ------------------------------------------------------------------------
    virtual int get_Num_LD() const {return Num_LD;}

    // ------------------------------------------------------------------------
    // determine whether a given index belongs to the LDN vector
    // ------------------------------------------------------------------------
    virtual bool is_inLDN( const int &ii ) const 
    { return VEC_T::is_invec(LDN, ii); }

    // ------------------------------------------------------------------------
    // get number of nodes belonging to the local partition
    // ------------------------------------------------------------------------
    virtual int get_num_local_node() const
    {return num_local_node;}

    // ------------------------------------------------------------------------
    // get number of cells belonging to the local partition
    // ------------------------------------------------------------------------
    virtual int get_num_local_cell() const
    {return num_local_cell;}

    // ------------------------------------------------------------------------
    // get the number of nodes per element, or nLocBas
    // ------------------------------------------------------------------------
    virtual int get_cell_nLocBas() const
    {return cell_nLocBas;}

    // ------------------------------------------------------------------------
    // access coordinates of a node in the local partition by indexing
    // the local_pt_xyz array
    // 0 <= ii < num_local_node
    // Note: make sure num_local_cell > 0 before using this get function
    // ------------------------------------------------------------------------
    virtual Vector_3 get_local_pt_xyz( const int &ii) const
    {return local_pt_xyz[ii];}

    // ------------------------------------------------------------------------
    // access coordinates of a Local Dirichlet nodes in the local partition by indexing
    // the LDN_pt_xyz array
    // 0 <= ii < Num_LD
    // Note: make sure Num_LD > 0 before using this get function
    virtual Vector_3 get_LDN_pt_xyz( const int &ii ) const
    {return LDN_pt_xyz[ii];}

    // ------------------------------------------------------------------------
    // access an element's IEN array by indexing the local_cell_ien array
    // 0 <= ii < cell_nLocBas x num_local_cell
    // Note: make sure num_local_cell > 0 before using this get function
    // ------------------------------------------------------------------------
    virtual int get_local_cell_ien( const int &ii ) const
    {return local_cell_ien[ii];}

    // ------------------------------------------------------------------------
    // access an element's IEN array by its local element index and surface element
    // node index
    // 0 <= ee < num_local_cell, 0 <= ii < cell_nLocBas
    // Note: make sure num_local_cell > 0 before using this get function
    // ------------------------------------------------------------------------
    virtual int get_local_cell_ien( 
        const int &ee, 
        const int &ii ) const
    { return local_cell_ien[ee * cell_nLocBas + ii]; }

    // ------------------------------------------------------------------------
    // get_ctrlPts_xyz: given surface element index eindex, return the
    // control point coordinates.
    // Users are responsible for allocating/deallocating memory for ctrl_(x/y/z)
    // surface element id: 0 <= eindex < num_local_cell;
    // ctrl_x/y/z : output coordinate arrays, each of length cell_nLocBas.
    // Note: make sure num_local_cell > 0 before using this get function
    // ------------------------------------------------------------------------
    virtual void get_ctrlPts_xyz( const int &eindex,
        double * const &ctrl_x, double * const &ctrl_y,
        double * const &ctrl_z ) const;

    // ------------------------------------------------------------------------
    // get_SIEN: given surface element index eindex, return the IEN.
    // Users are responsible for allocating/deallocating memory for sien
    // eindex : 0 <= eindex < num_local_cell
    // sien : length cell_nLocBas.
    // Note: make sure num_local_cell > 0 before using this get function
    // ------------------------------------------------------------------------
    virtual void get_SIEN( const int &eindex,
        int * const &sien ) const;

    virtual std::vector<int> get_SIEN( const int &eindex ) const;


  private:
    // Number of local Dirichlet nodes. Length num_nbc
    int Num_LD;

    // Global indices of the local Dirichlet nodes (rotated BC)
    // Length is Num_LD for 0<= ii < num_nbc
    std::vector<int> LDN;

    // Number of local nodes & cells, and cell's nLocBas. Length num_nbc
    int num_local_node, num_local_cell, cell_nLocBas;

    // Nodal coordinates of all local nodes
    // Length: num_local_node
    std::vector<Vector_3> local_pt_xyz;

    // Nodal coordinates of local Dirichlet nodes' coordinates
    std::vector<Vector_3> LDN_pt_xyz;

    // Surface IEN array of all local elements
    // Length: cell_nLocBas x num_local_cell
    std::vector<int> local_cell_ien;

    // Indices of all local nodes in the local volumetric partition's
    // local_to_global array
    // Length: num_local_node
    std::vector<int> local_node_pos;
  
    // ------------------------------------------------------------------------ 
    // Disallow default constructor 
    ALocal_RotatedBC() = delete; 

};

#endif
