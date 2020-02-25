#ifndef VISDATAPREP_NS_3D_HPP
#define VISDATAPREP_NS_3D_HPP
// ==================================================================
// VisDataPrep_NS_3D.hpp
//
// This is the data preparation for 3D Navier-Stokes (using tet 
// element.)
//
// Date Created: June 27 2017
// ==================================================================
#include "IVisDataPrep.hpp"

class VisDataPrep_NS_3D : public IVisDataPrep
{
  public:
    VisDataPrep_NS_3D();

    virtual ~VisDataPrep_NS_3D();

    virtual int get_ptarray_size() const {return 2;}

    virtual int get_ptarray_comp_length( const int &ii ) const
    { return pt_array_len[ii]; }

    virtual void get_pointArray(
        const std::string solution_file_name,
        const std::string analysis_node_mapping_file,
        const std::string post_node_mapping_file,
        const APart_Node * const &nNode_ptr,
        const IAGlobal_Mesh_Info * const &gInfo_ptr,
        const int &input_dof,
        double ** &pointArrays ) const;
  
  private:
    std::vector<int> pt_array_len;
};

#endif
