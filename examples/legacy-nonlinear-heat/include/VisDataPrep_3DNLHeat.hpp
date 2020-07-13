#ifndef VISDATAPREP_3DNLHEAT_HPP
#define VISDATAPREP_3DNLHEAT_HPP
// ==================================================================
// VisDataPrep_3DNLHeat.hpp
// This is the data preparation for 3D Nonlinear Heat problem.
//
// Dec. 17 2013  
// ==================================================================

#include "IVisDataPrep.hpp"

class VisDataPrep_3DNLHeat : public IVisDataPrep
{
  public:
    VisDataPrep_3DNLHeat();
    virtual ~VisDataPrep_3DNLHeat();

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
