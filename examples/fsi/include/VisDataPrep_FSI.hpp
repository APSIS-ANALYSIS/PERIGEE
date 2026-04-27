#ifndef VISDATAPREP_FSI_HPP
#define VISDATAPREP_FSI_HPP
// ============================================================================
// VisDataPrep_FSI.hpp
// 
// This is the data preparation for 3D FSI.
//
// Date: Jan 17 2022
// ============================================================================
#include "IVisDataPrep.hpp"

class VisDataPrep_FSI : public IVisDataPrep
{
  public:
    VisDataPrep_FSI();

    virtual ~VisDataPrep_FSI() = default;

    virtual int get_ptarray_size() const {return 3;}

    virtual int get_ptarray_comp_length( const int &ii ) const
    { return pt_array_len[ii]; }

    virtual void get_pointArray(
        const std::string &disp_solution_file_name,
        const std::string &pres_solution_file_name,
        const std::string &velo_solution_file_name,
        const std::vector<int> &an_v_node_mapping,
        const std::vector<int> &an_p_node_mapping,
        const std::vector<int> &pn_v_node_mapping,
        const std::vector<int> &pn_p_node_mapping,
        const APart_Node * const &pNode_v,
        const APart_Node * const &pNode_p,
        double ** &pointArrays ) const;

  private:
    std::vector<int> pt_array_len;
};

#endif
