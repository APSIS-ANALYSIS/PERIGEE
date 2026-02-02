#ifndef VISDATAPREP_SOLID_HPP
#define VISDATAPREP_SOLID_HPP
// ============================================================================
// VisDataPrep_Solid.hpp
//
// Data preparation routine for hyperelastic solid visualization.
//
// Date: Feb. 01 2026
// ============================================================================
#include "IVisDataPrep.hpp"

class VisDataPrep_Solid : public IVisDataPrep
{
  public:
    VisDataPrep_Solid();

    virtual ~VisDataPrep_Solid() = default;

    virtual int get_ptarray_size() const {return 3;}

    virtual int get_ptarray_comp_length( const int &ii ) const
    { return pt_array_len[ii]; }

    virtual void get_pointArray(
        const std::vector<std::string> solution_file_names,
        const std::string analysis_node_mapping_file,
        const std::string post_node_mapping_file,
        const APart_Node * const &nNode_ptr,
        const int &input_nfunc,
        double ** &pointArrays ) const;

  private:
    std::vector<int> pt_array_len;
};

#endif
