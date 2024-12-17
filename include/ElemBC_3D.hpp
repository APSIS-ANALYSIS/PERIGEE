#ifndef ELEMBC_3D_HPP
#define ELEMBC_3D_HPP
// ============================================================================
// ElemBC_3D.hpp
//
// This is an instantiation of ElemBC for 3D problems. It records the elemental 
// bc information for separate surface files.
//
// Date: Sep. 7 2023
// ============================================================================
#include "ElemBC.hpp"

class ElemBC_3D : public ElemBC
{
  public:
    // Default constructor: prescribe an elembc class with no surface mesh data.
    ElemBC_3D( const FEType &elemtype );

    ElemBC_3D( const std::vector<std::string> &vtkfileList, const FEType &elemtype );

    virtual ~ElemBC_3D() = default;

    virtual int get_num_ebc() const {return num_ebc;}

    virtual int get_num_node(const int &ebc_id) const {return num_node[ebc_id];}

    virtual int get_num_cell(const int &ebc_id) const {return num_cell[ebc_id];}

    virtual int get_cell_nLocBas(const int &ebc_id) const {return cell_nLocBas[ebc_id];}

    virtual double get_pt_xyz(const int &ebc_id, const int &node, const int &dir) const
    {return pt_xyz[ebc_id][3*node+dir];}

    virtual int get_ien(const int &ebc_id, const int &cell, const int &lnode) const
    {return sur_ien[ebc_id][ cell_nLocBas[ebc_id] * cell + lnode ];}

    virtual int get_global_node(const int &ebc_id, const int &node_index) const
    {return global_node[ebc_id][node_index];}

    virtual void get_global_node(const int &ebc_id, std::vector<int> &out) const
    {out = global_node[ebc_id];}

    virtual std::vector<int> get_global_node(const int &ebc_id) const
    {return global_node[ebc_id];}

    virtual int get_global_cell(const int &ebc_id, const int &cell_index) const
    {return global_cell[ebc_id][cell_index];}

    virtual void print_info() const;

    // For linear element (type Tet4), the face node numbering is
    //   Tet-Face-0 : Node 1 2 3
    //   Tet-Face-1 : Node 0 3 2
    //   Tet-Face-2 : Node 0 1 3
    //   Tet-Face-3 : Node 0 2 1
    // For quadratic element (type Tet10), the face node numbering is
    //   Tet-Face-0 : Node 1 2 3 5 9 8
    //   Tet-Face-1 : Node 0 3 2 7 9 6
    //   Tet-Face-2 : Node 0 1 3 4 8 7
    //   Tet-Face-3 : Node 0 2 1 6 5 4
    // For trilinear element (type Hex8), the face node numbering is
    //   Hex-Face-0 : Node 0 3 2 1
    //   Hex-Face-1 : Node 4 5 6 7
    //   Hex-Face-2 : Node 0 1 5 4
    //   Hex-Face-3 : Node 1 2 6 5
    //   Hex-Face-4 : Node 2 3 7 6
    //   Hex-Face-5 : Node 0 4 7 3
    // For triquadratic element (type Hex27), the face node numbering is
    //   Hex-Face-0 : Node 0 3 2 1 11 10 9 8 24
    //   Hex-Face-1 : Node 4 5 6 7 12 13 14 15 25
    //   Hex-Face-2 : Node 0 1 5 4 8 17 12 16 22
    //   Hex-Face-3 : Node 1 2 6 5 9 18 13 17 21
    //   Hex-Face-4 : Node 2 3 7 6 10 19 14 18 23
    //   Hex-Face-5 : Node 0 4 7 3 16 15 19 11 20
    virtual void resetSurIEN_outwardnormal( const IIEN * const &VIEN );

    // Access the data in ElemBC_3D_outflow, outward normal vector
    virtual Vector_3 get_normal_vec( const int &ebc_id ) const
    {SYS_T::commPrint("Warning: get_normal_vec is not implemented. \n"); return Vector_3();}

    // Access the data in ElemBC_3D_outflow, basis surface integration
    virtual std::vector<double> get_intNA( const int &ebc_id ) const
    {SYS_T::commPrint("Warning: get_intNA is not implemented.\n"); return {};}

    // Access the data in ElemBC_3D_wall_turbulence, wall model type
    virtual int get_wall_model_type() const
    {SYS_T::commPrint("Warning: get_wall_model_type is not implemented. \n"); return -1;}

    // Access the data in ElemBC_3D_wall_turbulence, face id of volume element
    virtual int get_faceID( const int &cell_index ) const
    {SYS_T::commPrint("Warning: get_faceID is not implemented. \n"); return {};}

  protected:
    const FEType elem_type;
    const int num_ebc;
    
    std::vector<int> num_node;     // length num_ebc
    std::vector<int> num_cell;     // length num_ebc
    std::vector<int> cell_nLocBas; // length num_ebc

    // num_ebc times 3 x num_node[ii] in size
    std::vector< std::vector<double> > pt_xyz;

    // num_ebc times cell_nLocBas[ii] x num_cell[ii] in size
    std::vector< std::vector<int> > sur_ien;

    // num_ebc times num_node[ii] in size
    std::vector< std::vector<int> > global_node;

    // num_ebc times num_cell[ii] in size
    std::vector< std::vector<int> > global_cell;

    // ------------------------------------------------------------------------
    // Disallow default constructor
    ElemBC_3D() = delete;

    // Reset function for the IEN array of different element types.
    void reset501IEN_outwardnormal( const IIEN * const &VIEN );

    void reset502IEN_outwardnormal( const IIEN * const &VIEN );

    void reset601IEN_outwardnormal( const IIEN * const &VIEN );

    void reset602IEN_outwardnormal( const IIEN * const &VIEN );
};

#endif
