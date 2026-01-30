#ifndef VISDATAPREP_TRANSPORT_HPP
#define VISDATAPREP_TRANSPORT_HPP
// ==================================================================
// VisDataPrep_Transport.hpp
//
// This is the data preparation for visualizing transport problems.
//
// Date Created: Nov. 13 2023
// ==================================================================
#include "IVisDataPrep.hpp"

class VisDataPrep_Transport : public IVisDataPrep
{
  public:
    VisDataPrep_Transport();

    virtual ~VisDataPrep_Transport() = default;

    // Return the number of physical fields to be read from solution
    // vector
    virtual int get_ptarray_size() const {return 1;}
   
    // Return the number of components for each physical field
    virtual int get_ptarray_comp_length( const int &ii ) const
    { return pt_array_len[ii]; }

    virtual void get_pointArray(
        const std::string solution_file_name,
        const std::vector<int> &analysis_node_mapping,
        const std::vector<int> &post_node_mapping,
        const APart_Node * const &nNode_ptr,
        const int &input_nfunc,
        const int &input_dof,
        double ** &solArrays ) const;

  private:
    std::vector<int> pt_array_len;
};

#endif
