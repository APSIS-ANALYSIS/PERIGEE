#ifndef ELEMBC_2D_HPP
#define ELEMBC_2D_HPP
// ==================================================================
// ElemBC_2D.hpp
//
// This is an instantiation of ElemBC for 2D problems.
// 
// The boundary is defined as a 2-1 = 1 dimensional manifold. In this
// class, each element boundary is identified with an index ebc_id,
// 0 <= ebc_id < num_ebc.
//
// num_node[ebc_id] : number of nodes on this boundary
// num_cell[ebc_id] : number of cells on this boundary
// cell_nLocBas[ebc_id] : number of local basis functions in cells
//
// pt_xyz[ebc_id] : the list of x-y-z coordinates of the nodes on 
//                  this boundary, with length 3 x num_node[ebc_id]
//
// cell_interior_pt[ebc_id] : the list of x-y-z coordindates of a point 
//                  in the interior of the surface element w.r.t. the
//                  2D domain. This will help define the outward normal
//                  vector
//
// tri_ien[ebc_id] : the list of local IEN array.
//
// global_node[ebc_id] : the list of the local nodes' global indices
//
// global_cell[ebc_id] : the list of the surface 2D element indices
//                       that the 1D cell belongs to. 
//
// Date created: Nov. 17 2017
// ==================================================================
#include "Vec_Tools.hpp"
#include "ElemBC.hpp"

class ElemBC_2D : public ElemBC
{
  public:
    ElemBC_2D( const std::vector<int> &edge_idx,
        const std::vector<int> &in_num_cell,
        const std::vector<int> &in_num_node,
        const std::vector<int> &in_nlocbas,
        const std::vector< std::vector<double> > &in_pt_coor,
        const std::vector< std::vector<int> > &in_ien_loc,
        const std::vector< std::vector<int> > &in_bcpt,
        const std::vector< std::vector<int> > &edge2elem,
        const IIEN * const &VIEN,
        const std::vector<double> &vctrlPts );

    virtual ~ElemBC_2D();

    virtual int get_num_ebc() const {return num_ebc;}

    // 0 <= ebc_id < num_ebc
    virtual int get_num_node(const int &ebc_id) const
    {return num_node[ebc_id];}

    virtual int get_num_cell(const int &ebc_id) const
    {return num_cell[ebc_id];}

    virtual int get_cell_nLocBas(const int &ebc_id) const
    {return cell_nLocBas[ebc_id];}

    // 0 <= node < num_node[ebc_id]
    // dir = 0, 1, 2
    virtual double get_pt_xyz(const int &ebc_id, const int &node,
        const int &dir) const
    {return pt_xyz[ebc_id][3*node+dir];}

    // 0 <= cell < num_cell[ebc_id]
    // 0 <= lnode < cell_nLocBas[ebc_id]
    virtual int get_ien(const int &ebc_id, const int &cell,
        const int &lnode) const
    {return tri_ien[ebc_id][cell_nLocBas[ebc_id]*cell + lnode];}

    // 0 <= node_index < num_node[ebc_id]
    virtual int get_global_node(const int &ebc_id, const int &node_index) const
    {return global_node[ebc_id][node_index];}

    // 0 <= cell_index < num_cell[ebc_id]
    virtual int get_global_cell(const int &ebc_id, const int &cell_index) const
    {return global_cell[ebc_id][cell_index];}

    // 0 <= cell_index < num_cell[ebc_id]
    virtual double get_intpt_xyz(const int &ebc_id, const int &cell_index,
        const int &dir) const
    {return cell_interior_pt[ebc_id][3*cell_index+dir];}

    virtual void print_info() const;

  protected:
    ElemBC_2D() {};

    int num_ebc;

    int * num_node;     // length num_ebc
    int * num_cell;     // length num_ebc
    int * cell_nLocBas; // length num_ebc

    // {num_ebc} times {3 x num_node[ii]} in size
    std::vector< std::vector<double> > pt_xyz;

    // {num_ebc} times {3 x num_cell[ii]} in size
    // stores the interior point of the boundary cell w.r.t. the surface element
    // this will help define the outward normal
    std::vector< std::vector<double> > cell_interior_pt;

    // {num_ebc} times {cell_nLocBas[ii] x num_cell[ii]} in size
    std::vector< std::vector<int> > tri_ien;

    // {num_ebc} times {num_node[ii]} in size
    std::vector< std::vector<int> > global_node;

    // {num_ebc} times {num_cell[ii]} in size
    std::vector< std::vector<int> > global_cell;
};

#endif
