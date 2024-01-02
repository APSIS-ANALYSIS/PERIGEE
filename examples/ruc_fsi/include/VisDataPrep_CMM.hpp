#ifndef VISDATAPREP_CMM_HPP
#define VISDATAPREP_CMM_HPP
// ==================================================================
// VisDataPrep_CMM.hpp
//
// This is the data preparation for visualizing NS problems.
//
// Author: Ju Liu
// Date Created: Aug 5 2017
// ==================================================================
#include "IVisDataPrep.hpp"

class VisDataPrep_CMM : public IVisDataPrep
{
  public:
    VisDataPrep_CMM();

    virtual ~VisDataPrep_CMM() = default;

    // Return the number of physical fields to be read from solution
    // vector
    virtual int get_ptarray_size() const {return 3;}
   
    // Return the number of components for each physical field
    virtual int get_ptarray_comp_length( const int &ii ) const
    { return pt_array_len[ii]; }

    virtual void get_pointArray(
        const std::vector<std::string> solution_file_names,
        const std::string analysis_node_mapping_file,
        const std::string post_node_mapping_file,
        const APart_Node * const &nNode_ptr,
        const int &in_nfunc,
        double ** &solArrays ) const;

  private:
    std::vector<int> pt_array_len;
};

#endif
