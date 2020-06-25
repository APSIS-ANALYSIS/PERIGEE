#ifndef VISDATAPREP_LAP_HPP
#define VISDATAPREP_LAP_HPP
// ==================================================================
// VisDataPrep_Lap.hpp
//
// This is the data preparation for visualization the laplace eqn.
//
// Date: June 5 2020
// ==================================================================
#include "IVisDataPrep.hpp"

class VisDataPrep_Lap : public IVisDataPrep
{
  public:
    VisDataPrep_Lap();

    virtual ~VisDataPrep_Lap();

    virtual int get_ptarray_size() const {return 1;}

    virtual int get_ptarray_comp_length( const int &ii ) const
    { return pt_array_len[ii]; }

    virtual void get_pointArray(
        const std::string solution_file_name,
        const std::string analysis_node_mapping_file,
        const std::string post_node_mapping_file,
        const APart_Node * const &nNode_ptr,
        const IAGlobal_Mesh_Info * const &gInfo_ptr,
        const int &input_dof,
        double ** &solArrays ) const;

  private:
    std::vector<int> pt_array_len;
};

#endif
