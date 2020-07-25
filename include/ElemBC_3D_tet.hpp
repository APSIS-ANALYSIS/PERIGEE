#ifndef ELEMBC_3D_TET_HPP
#define ELEMBC_3D_TET_HPP
// ==================================================================
// ElemBC_3D_tet.hpp
//
// This is an instantiation of ElemBC for 3D problems. It records the
// elemental bc information for separate surface files.
//
// Date: Feb. 6 2020
// ==================================================================
#include "ElemBC.hpp"
#include "Tet_Tools.hpp"

class ElemBC_3D_tet : public ElemBC
{
  public:
    ElemBC_3D_tet( const std::vector<std::string> &vtkfileList,
       const int &elemtype=501 );

    virtual ~ElemBC_3D_tet();

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
    {return tri_ien[ebc_id][ cell_nLocBas[ebc_id] * cell + lnode ];}

    virtual int get_global_node(const int &ebc_id, const int &node_index) const
    {return global_node[ebc_id][node_index];}

    virtual void get_global_node(const int &ebc_id, std::vector<int> &out) const
    {out = global_node[ebc_id];}

    virtual int get_global_cell(const int &ebc_id, const int &cell_index) const
    {return global_cell[ebc_id][cell_index];}

    virtual void print_info() const;

    // For linear element (type 501), the face node numbering is
    //   Tet-Face-0 : Node 1 2 3
    //   Tet-Face-1 : Node 0 3 2
    //   Tet-Face-2 : Node 0 1 3
    //   Tet-Face-3 : Node 0 2 1
    // For quadratic element (type 502), the face node numbering is
    //   Tet-Face-0 : Node 1 2 3 5 9 8
    //   Tet-Face-1 : Node 0 3 2 7 9 6
    //   Tet-Face-2 : Node 0 1 3 4 8 7
    //   Tet-Face-3 : Node 0 2 1 6 5 4
    virtual void resetTriIEN_outwardnormal( const IIEN * const &VIEN );

    // Access the data in ElemBC_3D_tet_outflow, outward normal vector
    virtual void get_normal_vec( const int &ebc_id, double &out_nx,
        double &out_ny, double &out_nz ) const
    {SYS_T::commPrint("Warning: get_normal_vec is not implemented. \n");}

    // Access the data in ElemBC_3D_tet_outflow, basis surface integration
    virtual void get_intNA( const int &ebc_id, std::vector<double> &fintNA ) const
    {SYS_T::commPrint("Warning: get_intNA is not implemented. \n");}

    // Access the data in ElemBC_3D_tet_wall, wall thickness used in CMM
    virtual void get_wall_thickness( const int &ebc_id, std::vector<double> &th ) const
    {SYS_T::commPrint("Warning: get_wall_thickness is not implemented. \n");}

    // Access the data in ElemBC_3D_tet_wall, wall youngs modulus used in CMM
    virtual void get_wall_youngsmod( const int &ebc_id, std::vector<double> &E ) const
    {SYS_T::commPrint("Warning: get_wall_youngsmod is not implemented. \n");}

    // write variable wall properties from ElemBC_3D_tet_wall for CMM
    virtual void write_wall_prop( const int &ebc_id ) const
    {SYS_T::commPrint("Warning: write_wall_prop is not implemented. \n");}

  protected:
    const int elem_type;
     
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
