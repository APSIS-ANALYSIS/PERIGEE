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
#include "FEType.hpp"

class ElemBC_3D : public ElemBC
{
  public:
    // Default constructor: prescribe an elembc class with no surface mesh data.
    ElemBC_3D( const FEType &in_elemtype );

    ElemBC_3D( const std::vector<std::string> &vtkfileList, const FEType &in_elemtype );

    ~ElemBC_3D() override = default;

    int get_num_ebc() const override {return num_ebc;}

    int get_num_node(int ebc_id) const override {return num_node[ebc_id];}

    int get_num_cell(int ebc_id) const override {return num_cell[ebc_id];}

    int get_cell_nLocBas(int ebc_id) const override {return cell_nLocBas[ebc_id];}

    double get_pt_xyz(int ebc_id, int node, int dir) const override
    {return pt_xyz[ebc_id][3*node+dir];}

    int get_ien(int ebc_id, int cell, int lnode) const override
    {return sur_ien[ebc_id][ cell_nLocBas[ebc_id] * cell + lnode ];}

    int get_global_node(int ebc_id, int node_index) const override
    {return global_node[ebc_id][node_index];}

    void get_global_node(int ebc_id, std::vector<int> &out) const override
    {out = global_node[ebc_id];}

    std::vector<int> get_global_node(int ebc_id) const override
    {return global_node[ebc_id];}

    int get_global_cell(int ebc_id, int cell_index) const override
    {return global_cell[ebc_id][cell_index];}

    void print_info() const override;

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
    void resetSurIEN_outwardnormal( const IIEN * const &VIEN ) override;

    // Access the data in ElemBC_3D_outflow, outward normal vector
    Vector_3 get_normal_vec( int ebc_id ) const override
    {SYS_T::commPrint("Warning: get_normal_vec is not implemented. \n"); return Vector_3();}

    // Access the data in ElemBC_3D_outflow, basis surface integration
    std::vector<double> get_intNA( int ebc_id ) const override
    {SYS_T::commPrint("Warning: get_intNA is not implemented.\n"); return {};}

    // Access the data in ElemBC_3D_wall_turbulence, wall model type
    int get_wall_model_type() const override
    {SYS_T::commPrint("Warning: get_wall_model_type is not implemented. \n"); return -1;}

    // Access the data in ElemBC_3D_wall_turbulence, face id of volume element
    int get_faceID( int cell_index ) const override
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
