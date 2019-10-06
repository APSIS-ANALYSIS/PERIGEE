#ifndef VISDATAPREP_2DLINEARELASTIC_HPP
#define VISDATAPREP_2DLINEARELASTIC_HPP
// ==================================================================
// VisDataPrep_2DLinearElastic.hpp
// This is the data preparation for 2D Linear Elastic problem.
//
// Date: Sept 12 2016
// ==================================================================

#include "IVisDataPrep.hpp"

class VisDataPrep_2DLinearElastic : public IVisDataPrep
{
  public:
    VisDataPrep_2DLinearElastic();

    virtual ~VisDataPrep_2DLinearElastic();

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
