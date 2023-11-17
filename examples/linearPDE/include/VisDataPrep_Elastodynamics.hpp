#ifndef VISDATAPREP_ELASTODYNAMICS_HPP
#define VISDATAPREP_ELASTODYNAMICS_HPP
// ==================================================================
// VisDataPrep_Elastodynamics.hpp
//
// This is the data preparation for visualizing elastodynamics 
// problems.
//
// Date Created: Nov. 5 2023
// ==================================================================
#include "IVisDataPrep.hpp"

class VisDataPrep_Elastodynamics : public IVisDataPrep
{
  public:
    VisDataPrep_Elastodynamics();

    virtual ~VisDataPrep_Elastodynamics() = default;

    // Return the number of physical fields to be read from solution
    // vector
    virtual int get_ptarray_size() const {return 1;}
   
    // Return the number of components for each physical field
    virtual int get_ptarray_comp_length( const int &ii ) const
    { return pt_array_len[ii]; }

    virtual void get_pointArray(
        const std::string solution_file_name,
        const std::string analysis_node_mapping_file,
        const std::string post_node_mapping_file,
        const APart_Node * const &nNode_ptr,
        const int &input_nfunc,
        const int &input_dof,
        double ** &solArrays ) const;

  private:
    std::vector<int> pt_array_len;
};

#endif
