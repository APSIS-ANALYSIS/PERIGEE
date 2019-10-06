#ifndef VISDATAPREP_3D_NS_HPP
#define VISDATAPREP_3D_NS_HPP
// ==================================================================
// VisDataPrep_3DNS.hpp
// This is the data preparation routine for 3D Navier-Stokes equations.
//
// Date: March 13 2015
// ==================================================================

#include "IVisDataPrep.hpp"

class VisDataPrep_3DNS : public IVisDataPrep
{
  public:
    VisDataPrep_3DNS();

    virtual ~VisDataPrep_3DNS();

    virtual void get_pointArray(
        const std::string solution_file_name,
        const std::string analysis_node_mapping_file,
        const std::string post_node_mapping_file,
        const APart_Node * const &nNode_ptr,
        const IAGlobal_Mesh_Info * const &gInfo_ptr,
        const int &input_dof,
        double ** &pointArrays ) const;
};

#endif
