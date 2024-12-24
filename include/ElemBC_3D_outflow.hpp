#ifndef ELEMBC_3D_OUTFLOW_HPP
#define ELEMBC_3D_OUTFLOW_HPP
// ============================================================================
// ElemBC_3D_outflow.hpp
//
// A derived class from the ElemBC_3D.hpp
//
// This class has a vector with length num_node for each ebc id that stores 
// int_{Gamma} N_A dS
//
// Date: Sep. 28 2023
// ============================================================================
#include "ElemBC_3D.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FEAElement_Quad9_3D_der0.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Quad.hpp"

class ElemBC_3D_outflow : public ElemBC_3D
{
  public:
    // Constructor for outlet boundaries that will be associated
    // with reduced order models in the analysis code. 
    // \para vtkfileList: the list of the vtp files for the outlet surfaces
    // \para outlet_normal_vec: the outward normal vectors for the surfaces
    ElemBC_3D_outflow( const std::vector<std::string> &vtkfileList,
       const std::vector< Vector_3 > &outlet_normal_vec,
       const FEType &in_elemtype );

    virtual ~ElemBC_3D_outflow() = default;

    virtual void print_info() const;

    virtual Vector_3 get_normal_vec( const int &ebc_id ) const
    { return outNormal[ebc_id]; }

    virtual std::vector<double> get_intNA( const int &ebc_id ) const
    { return intNA[ebc_id]; }

  private:
    // Disallow the default constructor
    ElemBC_3D_outflow() = delete;

    // It stores the surface integral of each nodal basis function on the 
    // outlet surfaces.
    // Dimension is num_ebc x num_node[ii].
    std::vector< std::vector<double> > intNA;

    // It explicitly stores the outlet face normal vector.
    // Dimension num_ebc
    std::vector< Vector_3 > outNormal;
};

#endif
