#ifndef VISDATAPREP_2DNLHEAT_HPP
#define VISDATAPREP_2DNLHEAT_HPP
// ==================================================================
// VisDataPrep_2DNLHeat.hpp
// This is the data preparation for 2D Nonlinear Heat problem.
//
// Date: April 20 2014
// ==================================================================

#include "IVisDataPrep.hpp"

class VisDataPrep_2DNLHeat : public IVisDataPrep
{
  public:
    VisDataPrep_2DNLHeat();
    virtual ~VisDataPrep_2DNLHeat();

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
