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

  protected:
    ElemBC_3D_tet() {}; // Disallow default constructor
    
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
