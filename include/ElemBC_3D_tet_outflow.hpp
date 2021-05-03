#ifndef ELEMBC_3D_TET_OUTFLOW_HPP
#define ELEMBC_3D_TET_OUTFLOW_HPP
// ==================================================================
// ElemBC_3D_tet_outflow.hpp
//
// A derived class from the ElemBC_3D_tet.hpp
//
// This class has a vector with length num_node for each ebc id that
// stores int_{Gamma} N_A dS
//
// Date: Feb. 6 2020
// ==================================================================
#include "ElemBC_3D_tet.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "QuadPts_Gauss_Triangle.hpp"

class ElemBC_3D_tet_outflow : public ElemBC_3D_tet
{
  public:
    // Constructor for outlet boundaries that will be associated
    // with reduced order models in the analysis code. 
    // \para vtkfileList: the list of the vtp files for the outlet surfaces
    // \para outlet_normal_vec: the outward normal vectors for the surfaces
    ElemBC_3D_tet_outflow( const std::vector<std::string> &vtkfileList,
       const std::vector< std::vector<double> > &outlet_normal_vec,
       const int &elemtype = 501 );

    virtual ~ElemBC_3D_tet_outflow();

    virtual void print_info() const;

    virtual void get_normal_vec( const int &ebc_id, double &out_nx,
        double &out_ny, double &out_nz ) const
    {
      out_nx = outNormal[ebc_id][0];
      out_ny = outNormal[ebc_id][1];
      out_nz = outNormal[ebc_id][2];
    }

    virtual std::vector<double> get_intNA( const int &ebc_id ) const
    { return intNA[ebc_id]; }

  private:
    std::vector< std::vector<double> > intNA;
    std::vector< std::vector<double> > outNormal;
};

#endif
