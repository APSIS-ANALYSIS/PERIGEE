#ifndef ALOCAL_MOVINGBC_HPP
#define ALOCAL_MOVINGBC_HPP
// ============================================================================
// ALocal_MovingBC.hpp
//
// Analysis-use data that is related to the moving surface nodes.
//
// Date Created: Jun 14 2024
// Author: Yujie Sun
// ============================================================================
#include "HDF5_Reader.hpp"
#include "Math_Tools.hpp"

class ALocal_MovingBC
{
  public:
    ALocal_MovingBC( const std::string &fileBaseName, const int &cpu_rank );

    virtual ~ALocal_MovingBC() = default;

    // ------------------------------------------------------------------------
    // Get the number of surfaces with prescribed moving
    // ------------------------------------------------------------------------
    virtual int get_num_nbc() const {return num_nbc;};

    // ------------------------------------------------------------------------
    // Get global index of a Dirichlet node in the local partition
    // 0 <= node < Num_LD[nbc_id]
    // Note: make sure that Num_LD[nbc_id] > 0 before calling this get function
    // ------------------------------------------------------------------------
    virtual int get_LDN( const int &nbc_id, const int &node ) const
    {return LDN[nbc_id][node];}

    // ------------------------------------------------------------------------
    // Get the number of Dirichlet nodes in the local partition
    // ------------------------------------------------------------------------
    virtual int get_Num_LD( const int &nbc_id ) const {return Num_LD[nbc_id];}

    // ------------------------------------------------------------------------
    // determine whether a given index belongs to the LDN[nbc_id] vector
    // ------------------------------------------------------------------------
    virtual bool is_inLDN( const int &nbc_id, const int &ii ) const 
    { return VEC_T::is_invec(LDN[nbc_id], ii); }

    // ------------------------------------------------------------------------
    // get number of nodes belonging to the local partition
    // ------------------------------------------------------------------------
    virtual int get_num_local_node( const int &nbc_id ) const
    {return num_local_node[nbc_id];}

    // ------------------------------------------------------------------------
    // get number of cells belonging to the local partition
    // ------------------------------------------------------------------------
    virtual int get_num_local_cell( const int &nbc_id ) const
    {return num_local_cell[nbc_id];}

    // ------------------------------------------------------------------------
    // get the number of nodes per element, or nLocBas[nbc_id]
    // ------------------------------------------------------------------------
    virtual int get_cell_nLocBas( const int &nbc_id ) const
    {return cell_nLocBas[nbc_id];}

    // ------------------------------------------------------------------------
    // access coordinates of a node in the local partition by indexing
    // the local_pt_xyz array
    // 0 <= ii < num_local_node[nbc_id]
    // Note: make sure num_local_cell[nbc_id] > 0 before using this get function
    // ------------------------------------------------------------------------
    virtual Vector_3 get_local_pt_xyz( const int &nbc_id, const int &ii) const
    {return local_pt_xyz[nbc_id][ii];}

    // ------------------------------------------------------------------------
    // access coordinates of a Local Dirichlet nodes in the local partition by indexing
    // the LDN_pt_xyz array
    // 0 <= ii < Num_LD[nbc_id]
    // Note: make sure Num_LD[nbc_id] > 0 before using this get function
    virtual Vector_3 get_LDN_pt_xyz( const int &nbc_id, const int &ii) const
    {return LDN_pt_xyz[nbc_id][ii];}

    // ------------------------------------------------------------------------
    // access an element's IEN array by indexing the local_cell_ien array
    // 0 <= ii < cell_nLocBas[nbc_id] x num_local_cell[nbc_id]
    // Note: make sure num_local_cell[nbc_id] > 0 before using this get function
    // ------------------------------------------------------------------------
    virtual int get_local_cell_ien( const int &nbc_id, const int &ii ) const
    {return local_cell_ien[nbc_id][ii];}

    // ------------------------------------------------------------------------
    // access an element's IEN array by its local element index and surface element
    // node index
    // 0 <= ee < num_local_cell[nbc_id], 0 <= ii < cell_nLocBas[nbc_id]
    // Note: make sure num_local_cell[nbc_id] > 0 before using this get function
    // ------------------------------------------------------------------------
    virtual int get_local_cell_ien( const int &nbc_id, const int &ee, 
        const int &ii ) const
    { return local_cell_ien[nbc_id][ee * cell_nLocBas[nbc_id] + ii]; }

    // ------------------------------------------------------------------------
    // get_ctrlPts_xyz: given surface element index eindex, return the
    // control point coordinates.
    // Users are responsible for allocating/deallocating memory for ctrl_(x/y/z)
    // surface element id: 0 <= eindex < num_local_cell[nbc_id];
    // ctrl_x/y/z : output coordinate arrays, each of length cell_nLocBas[nbc_id].
    // Note: make sure num_local_cell[nbc_id] > 0 before using this get function
    // ------------------------------------------------------------------------
    virtual void get_ctrlPts_xyz( const int &nbc_id, const int &eindex,
        double * const &ctrl_x, double * const &ctrl_y,
        double * const &ctrl_z ) const;

    // ------------------------------------------------------------------------
    // get_SIEN: given surface element index eindex, return the IEN.
    // Users are responsible for allocating/deallocating memory for sien
    // eindex : 0 <= eindex < num_local_cell[nbc_id]
    // sien : length cell_nLocBas[nbc_id].
    // Note: make sure num_local_cell[nbc_id] > 0 before using this get function
    // ------------------------------------------------------------------------
    virtual void get_SIEN( const int &nbc_id, const int &eindex,
        int * const &sien ) const;

    virtual std::vector<int> get_SIEN( const int &nbc_id, const int &eindex ) const;


  private:
    int num_nbc;

    // Number of local Dirichlet nodes. Length num_nbc
    std::vector<int> Num_LD;

    // Global indices of the local Dirichlet nodes (moving BC)
    // Length is num_nbc x Num_LD[ii] for 0<= ii < num_nbc
    std::vector< std::vector<int> > LDN;

    // Number of local nodes & cells, and cell's nLocBas. Length num_nbc
    std::vector<int> num_local_node, num_local_cell, cell_nLocBas;

    // Nodal coordinates of all local nodes
    // num_nbc times num_local_node[ii]
    std::vector< std::vector<Vector_3> > local_pt_xyz;

    // Nodal coordinates of local Dirichlet nodes' coordinates
    std::vector< std::vector<Vector_3> > LDN_pt_xyz;

    // Surface IEN array of all local elements
    // num_nbc times (cell_nLocBas[ii] x num_local_cell[ii])
    std::vector< std::vector<int> > local_cell_ien;

    // Indices of all local nodes in the local volumetric partition's
    // local_to_global array
    // num_nbc x num_local_node[ii]
    std::vector< std::vector<int> > local_node_pos;
  
    // ------------------------------------------------------------------------ 
    // Disallow default constructor 
    ALocal_MovingBC() = delete; 

};

#endif