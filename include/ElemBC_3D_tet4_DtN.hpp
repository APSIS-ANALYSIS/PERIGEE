#ifndef ELEMBC_3D_TET4_DTN_HPP
#define ELEMBC_3D_TET4_DTN_HPP
// ==================================================================
// ElemBC_3D_tet4_DtN.hpp
//
// A derived class based on the standard ElemBC_3D_Tet4.hpp. This
// class has a vector for each ebc with length equaling to the number
// of nodes for that ebc surface. This vector will store
//                int_{Gamma} N_A n_i dA.
// This vector will be used as a rank-1 modification to the consistent
// tangent matrix to account for the lumped parameter model coupling.
//
// Date: Nov. 6 2017
// Author: Ju Liu
// ==================================================================
#include "ElemBC_3D_tet4.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "HDF5_Writer.hpp"

class ElemBC_3D_tet4_DtN : public ElemBC_3D_tet4
{
  public:
    // Constructor: vtpfileList : the list of vtp files describing the
    //                            boundary surfaces that require a surface
    //                            integral
    //              dtn_bc_list:  the list of id's of the surfaces that
    //                            are Dirichlet-to-Neumann boundary.
    //                            vtpfileList[ dtn_bc_list[ii] ]
    //                            is the surface that is a DtN surface.
    //                            dtn_bc_list.size() gives the number of
    //                            DtN surfaces.
    ElemBC_3D_tet4_DtN( const std::vector<std::string> &vtpfileList,
       const std::vector<int> &dtn_bc_list );

    virtual ~ElemBC_3D_tet4_DtN();

    virtual void print_info() const;

  private:
    // prohibit the default constructor
    ElemBC_3D_tet4_DtN() {};

    // Number of the Dirichlet-to-Neumann boundary surfaces.
    int num_dtn_surface;

    // num_ebc times num_node[ii] in size, stores the int_{Gamma} N_A dA 
    // for each surface.
    std::vector< std::vector<double> > intNA;
};

#endif
