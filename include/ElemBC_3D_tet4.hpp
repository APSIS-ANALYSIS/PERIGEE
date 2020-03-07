#ifndef ELEMBC_3D_TET4_HPP
#define ELEMBC_3D_TET4_HPP
// ==================================================================
// ElemBC_3D_tet4.hpp
//
// This is an instantiation of ElemBC for 3D problems described by 
// linear tetrahedral elements. For linear tets, the boundary is
// described by triangle elements. This information for the triagnles
// can be obtained directly from the mesher, such as tetgen; and one
// can read these info from vtp files.
//
// The data structure is given as follows. The input is a list of .vtp
// file names; the length of the vtpfileList is the number of different
// boundary surfaces that requires boundary integration. This number 
// is also recorded as num_ebc;
//
// For each portion of the elemental boundary with index ii, 
// we record the number of nodes, 
//           the number of cells, and 
//           the nLocBas (essentially the size of the connectivity array)
// as num_node[ii], num_cell[ii], cell_nLocBas[ii].
//
// pt_xyz[ii] gives the nodal point's coordinates;
//
// tri_ien[ii] gives the cell's IEN array in a 1D array;
// 
// global_node[ii] gives the global nodal indices;
//
// global_cell[ii] gives the global cell index.
//
// Date: Jan. 10 2017
// ==================================================================
#include "ElemBC.hpp"
#include "Tet_Tools.hpp"

class ElemBC_3D_tet4 : public ElemBC
{
  public:
    ElemBC_3D_tet4( const std::vector<std::string> &vtpfileList );

    virtual ~ElemBC_3D_tet4();

    virtual int get_num_ebc() const {return num_ebc;}

    virtual int get_num_node(const int &ebc_id) const 
    {return num_node[ebc_id];}

    virtual int get_num_cell(const int &ebc_id) const
    {return num_cell[ebc_id];}

    virtual int get_cell_nLocBas(const int &ebc_id) const
    {return cell_nLocBas[ebc_id];}

    virtual double get_pt_xyz(const int &ebc_id, const int &node,
        const int &dir) const
    {return pt_xyz[ebc_id][3*node+dir];}

    virtual int get_ien(const int &ebc_id, const int &cell, 
        const int &lnode) const
    {return tri_ien[ebc_id][3*cell + lnode];}

    virtual int get_global_node(const int &ebc_id, const int &node_index) const
    {return global_node[ebc_id][node_index];}

    virtual void get_global_node(const int &ebc_id, 
        std::vector<int> &out) const
    {out = global_node[ebc_id];}
    
    virtual int get_global_cell(const int &ebc_id, const int &cell_index) const
    {return global_cell[ebc_id][cell_index];}

    virtual void print_info() const;

    // The surface element (triangle)'s IEN needs to be reset so that
    // the triangle's local IEN gives point ordering such that the 
    // vector 01 cross vector 02 gives the outward normal vector 
    // The tetrahedral face number is given by the opposing vertex's
    // number. The tet face 0 for linear tets, includes the node 1,2,3.
    // vector 1->2 cross vector 1->3 gives the outward normal vector.
    // A complete Tet-face node numbering is:
    //  Tet-Face-0 : Node 1 2 3
    //  Tet-Face-1 : Node 0 3 2
    //  Tet-Face-2 : Node 0 1 3
    //  Tet-Face-3 : Node 0 2 1
    // This reset function will use the VIEN to find the triangle's face
    // id, and correct the tri_ien if the numbering does not match the
    // above pattern.
    virtual void resetTriIEN_outwardnormal( const IIEN * const &VIEN );

    virtual void get_normal_vec( const int &ebc_id, double &out_nx,
        double &out_ny, double &out_nz ) const
    {SYS_T::commPrint("Warning: get_normal_vec is not implemented. \n");}

    virtual void get_intNA( const int &ebc_id, std::vector<double> &fintNA ) const
    {SYS_T::commPrint("Warning: get_intNA is not implemented. \n");}

  protected:
    // prohibit the default constructor
    ElemBC_3D_tet4() {};

    int num_ebc;
    int * num_node;     // length num_ebc
    int * num_cell;     // length num_ebc
    int * cell_nLocBas; // length num_ebc

    // num_ebc times 3 x num_node[ii] in size
    std::vector< std::vector<double> > pt_xyz;

    // num_ebc times cell_nLocBas[ii] x num_cell[ii] in size
    std::vector< std::vector<int> > tri_ien;

    // num_ebc times num_node[ii] in size
    std::vector< std::vector<int> > global_node;

    // num_ebc times num_cell[ii] in size
    std::vector< std::vector<int> > global_cell;
};

#endif
