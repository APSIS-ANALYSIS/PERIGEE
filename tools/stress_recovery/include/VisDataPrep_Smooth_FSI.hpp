#ifndef VISDATAPREP_SMOOTH_FSI_HPP
#define VISDATAPREP_SMOOTH_FSI_HPP
// ==================================================================
// VisDataPrep_Smooth_FSI.hpp
//
// This is the data preparation for visualizing smoothed solutions
//
// Date Created: Jan. 22 2024
// ==================================================================
#include "IVisDataPrep.hpp"

class VisDataPrep_Smooth_FSI : public IVisDataPrep
{
  public:
    VisDataPrep_Smooth_FSI();

    virtual ~VisDataPrep_Smooth_FSI() = default;

    // Return the number of physical fields to be read from solution
    // vector
    virtual int get_ptarray_size() const {return 3;}
   
    // Return the number of components for each physical field
    virtual int get_ptarray_comp_length( const int &ii ) const
    { return pt_array_len[ii]; }

    virtual void get_pointArray(
        const std::string &cauchy_solution_file_name,
        const std::string &disp_sol_file_name,
        const std::string &pres_solution_file_name,
        const std::string &an_v_node_mapping_file,
        const std::string &an_p_node_mapping_file,
        const std::string &pn_v_node_mapping_file,
        const std::string &pn_p_node_mapping_file,
        const APart_Node * const &pNode_v,
        const APart_Node * const &pNode_p,
        const int &input_nfunc_v,
        const int &input_nfunc_p,
        double ** &pointArrays ) const;

  private:
    std::vector<int> pt_array_len;
};

#endif