#ifndef VISDATAPREP_HYPERELASTIC_HPP
#define VISDATAPREP_HYPERELASTIC_HPP
// ============================================================================
// VisDataPrep_Hyperelastic.hpp
//
// This is the data preparation routine for HYPERELASTIC equations.
//
// Date: Jan 18 2022
// ============================================================================
#include "IVisDataPrep.hpp"

class VisDataPrep_Hyperelastic : public IVisDataPrep
{
  public:
    VisDataPrep_Hyperelastic( const bool &is_ref );

    virtual ~VisDataPrep_Hyperelastic() = default;

    virtual int get_ptarray_size() const {return 3;}

    virtual int get_ptarray_comp_length( const int &ii ) const
    { return pt_array_len[ii]; }

    virtual void get_pointArray(
        const std::string &disp_solution_file_name,
        const std::string &pres_solution_file_name,
        const std::string &velo_solution_file_name,
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

