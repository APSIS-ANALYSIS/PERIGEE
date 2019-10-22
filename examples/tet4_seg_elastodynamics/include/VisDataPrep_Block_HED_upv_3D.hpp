#ifndef VISDATAPREP_BLOCK_HED_UPV_3D_HPP
#define VISDATAPREP_BLOCK_HED_UPV_3D_HPP
// ==================================================================
// VisDataPrep_Block_HED_upv_3D.hpp
//
// This is the data preparation for 3D HyperElastoDynamics (HED) with
// block solver, and input u-p-v disp, pres, velo solution vectors.
//
// Author: Ju Liu
// Date created: Feb. 27 2018
// ==================================================================
#include "IVisDataPrep.hpp"

class VisDataPrep_Block_HED_upv_3D : public IVisDataPrep
{
  public:
    VisDataPrep_Block_HED_upv_3D( const bool &is_ref );


    virtual ~VisDataPrep_Block_HED_upv_3D();


    virtual int get_ptarray_size() const {return 3;}


    virtual int get_ptarray_comp_length( const int &ii ) const
    { return pt_array_len[ii]; }


    // solution_file_names store the names for the disp, pres, and 
    // velo.
    virtual void get_pointArray(
        const std::vector<std::string> solution_file_names,
        const std::string analysis_node_mapping_file,
        const std::string post_node_mapping_file,
        const APart_Node * const &nNode_ptr,
        const int &in_nfunc,
        double ** &pointArrays ) const;


  private:
    std::vector<int> pt_array_len;
};

#endif
