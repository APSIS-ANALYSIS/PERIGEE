#ifndef VISDATAPREP_SMOOTH_VOL_HPP
#define VISDATAPREP_SMOOTH_VOL_HPP
// ==================================================================
// VisDataPrep_Smooth_Vol.hpp
//
// This is the data preparation for visualizing smoothed solutions
//
// Date Created: Jan. 22 2024
// ==================================================================
#include "IVisDataPrep.hpp"

class VisDataPrep_Smooth_Vol : public IVisDataPrep
{
  public:
    VisDataPrep_Smooth_Vol();

    virtual ~VisDataPrep_Smooth_Vol() = default;

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