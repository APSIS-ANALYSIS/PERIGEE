#ifndef VISDATAPREP_ALE_NS_HPP
#define VISDATAPREP_ALE_NS_HPP
// ============================================================================
// VisDataPrep_ALE_NS.hpp
//
// This is the data preparation routine for 3D ALE-NS equations.
//
// Date: Jan 17 2022
// ============================================================================
#include "IVisDataPrep.hpp"

class VisDataPrep_ALE_NS : public IVisDataPrep
{
  public:
    VisDataPrep_ALE_NS();

    virtual ~VisDataPrep_ALE_NS() = default;

    virtual int get_ptarray_size() const {return 3;}

    virtual int get_ptarray_comp_length( const int &ii ) const
    { return pt_array_len[ii]; }

    virtual void get_pointArray(
        const std::string &disp_sol_file_name,
        const std::string &pres_sol_file_name,
        const std::string &velo_sol_file_name,
        const std::string &an_v_mapping_file,
        const std::string &an_p_mapping_file,
        const std::string &pn_v_mapping_file,
        const std::string &pn_p_mapping_file,
        const APart_Node * const &pNode_v,
        const APart_Node * const &pNode_p,
        const int &input_nfunc_v,
        const int &input_nfunc_p,
        double ** &pointArrays ) const;

  private:
    std::vector<int> pt_array_len;
};

#endif
