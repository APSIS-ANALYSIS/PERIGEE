#ifndef ELEMBC_3D_TET4_OUTFLOW_HPP
#define ELEMBC_3D_TET4_OUTFLOW_HPP
// ==================================================================
// ElemBC_3D_tet4_outflow.hpp
//
// A derived class based on the standard ElemBC_3D_tet4.hpp.
//
// This class has a vector with length num_node for each ebc id storing
//                    int_{Gamma} N_A dS.
// This vector will be used as a rank-one modification to the consistent
// tangent matrix to account for various lumped parameter models for
// the 0D-3D coupling.
//
// Date: June 5 2018
// Author: Ju Liu
// ==================================================================
#include "ElemBC_3D_tet4.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "HDF5_Writer.hpp"

class ElemBC_3D_tet4_outflow : public ElemBC_3D_tet4
{
  public:
    ElemBC_3D_tet4_outflow( const std::vector<std::string> &vtpfileList,
       const std::vector< std::vector<double> > &outlet_normal_vec );

    virtual ~ElemBC_3D_tet4_outflow();

    virtual void print_info() const;

    virtual void get_normal_vec( const int &ebc_id, double &out_nx,
        double &out_ny, double &out_nz ) const
    {
      out_nx = outNormal[ebc_id][0];
      out_ny = outNormal[ebc_id][1];
      out_nz = outNormal[ebc_id][2];
    }

    virtual void get_intNA( const int &ebc_id, 
        std::vector<double> &fintNA ) const
    { fintNA = intNA[ebc_id]; }

  private:
    // prohibit the default constructor
    ElemBC_3D_tet4_outflow() {};

    // intNA has dimension num_ebc x num_node[ii].
    // It stores the surface integral of each nodal basis function
    // on the corresponding surface.
    std::vector< std::vector<double> > intNA;

    // outNormal has dimension num_ebc x 3
    // It explicitly stores the outlet face normal vector.
    std::vector< std::vector<double> > outNormal;
};

#endif
