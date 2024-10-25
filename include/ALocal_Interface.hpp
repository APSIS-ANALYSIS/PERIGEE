#ifndef ALOCAL_INTERFACE_HPP
#define ALOCAL_INTERFACE_HPP
// ============================================================================
// ALocal_Interface.hpp
//
// FEM-analysis-use Local subdomain's interface.
// Our model may contains more than one interface, therefore this class is 
// similar with ALocal_EBC, but it uses the volume elements' info.
//
// The volume elements attached to an interface would be from the fixed volume
// and the rotated volume. We distribute the fixed volume elements to each
// partition, but store the rotated volume elements' information to every
// partition h5 file which involves the interface. 

// When calculate the integral on the interface, we go through the fixed elements,
// use the quadrature points of the fixed face, and search for the opposite points
// on the rotated face.
//
// Author: Xuanming Huang
// Date Created: Jun. 24  2024
// ============================================================================
#include <vector>
#include "PDNSolution.hpp"

class ALocal_Interface
{
  public:
    ALocal_Interface( const std::string &fileBaseName, const int &cpu_rank);

    ~ALocal_Interface() = default;

    int get_num_itf() const {return num_itf;}

    int get_nLocBas() const {return nLocBas;}

    int get_num_tag(const int &ii) const {return num_tag[ii];}

    int get_num_fixed_ele(const int &ii) const {return num_fixed_ele[ii];}

    virtual int get_num_fixed_ele(const int &ii) const
    {return num_local_fixed_ele[ii];}

    virtual int get_max_num_fixed_ele(const int &ii) const
    {return max_num_local_fixed_ele[ii];}

    virtual int get_fixed_ele(const int &ii, const int &ee) const
    {return local_fixed_ele[ii][ee];} // return "jj" for other functions

    int get_fixed_face_id(const int &ii, const int &jj) const
    {return fixed_ele_face_id[ii][jj];}

    int get_fixed_ele_tag(const int &ii, const int &jj) const
    {return fixed_ele_tag[ii][jj];}

    virtual int get_num_tagged_fixed_ele(const int &ii, const int &tag) const
    {return num_tagged_fixed_ele[ii][tag];}

    virtual int get_tagged_fixed_ele(const int &ii, const int &tag, const int &ee) const
    {return tagged_fixed_ele[ii][tag][ee];} // return "jj" for other functions

    virtual int get_fixed_lien(const int &ii, const int &kk) const
    {return fixed_lien[ii][kk];}

    int get_fixed_LID(const int &ii, const int &dof_index, const int &node) const
    {return fixed_LID[ii][dof_index * num_fixed_node[ii] + node];}

    virtual void get_fixed_ele_ctrlPts(const int &ii, const int &jj,
      double * const volctrl_x,  double * const volctrl_y,  double * const volctrl_z) const;

    virtual int get_max_num_rotated_ele(const int &ii) const
    {return max_num_local_rotated_ele[ii];}

    virtual int get_num_rotated_ele(const int &ii) const
    {return num_local_rotated_ele[ii];}

    virtual int get_rotated_ele(const int &ii, const int &ee) const
    {return local_rotated_ele[ii][ee];} // return "jj" for other functions

    virtual int get_rotated_face_id(const int &ii, const int &jj) const
    {return rotated_ele_face_id[ii][jj];}

    virtual int get_rotated_ele_tag(const int &ii, const int &jj) const
    {return rotated_ele_tag[ii][jj];}

    virtual int get_num_tagged_rotated_ele(const int &ii, const int &tag) const
    {return num_tagged_rotated_ele[ii][tag];}

    virtual int get_tagged_rotated_ele(const int &ii, const int &tag, const int &ee) const
    {return tagged_rotated_ele[ii][tag][ee];} // return "jj" for other functions

    virtual int get_rotated_lien(const int &ii, const int &kk) const
    {return rotated_lien[ii][kk];}

    int get_rotated_LID(const int &ii, const int &dof_index, const int &node) const
    {return rotated_LID[ii][dof_index * num_rotated_node[ii] + node];}

    virtual void get_rotated_ele_ctrlPts(const int &ii, const int &jj,
        double * const volctrl_x,  double * const volctrl_y,  double * const volctrl_z) const;

  protected:
    // the number of interfaces
    int num_itf;

    // the number of local basis function
    int nLocBas;

    int nqp_sur;

    // the number of local fixed volume elements of each interface
    // size: num_itf
    std::vector<int> num_local_fixed_ele;

    // max numbers of the local fixed volume elements among all parts
    // size: num_itf
    std::vector<int> max_num_local_fixed_ele;

    // the number of local rotated volume elements of each interface
    // size: num_itf
    std::vector<int> num_local_rotated_ele;
    
    // max numbers of the local rotated volume elements among all parts
    // size: num_itf
    std::vector<int> max_num_local_rotated_ele;

    // the number of the nodes from the fixed volume elements
    // size: num_itf
    std::vector<int> num_fixed_node;

    // the number of the nodes from the rotated volume elements
    // size: num_itf
    std::vector<int> num_rotated_node;

    // the number of interval tag of each pair of interfaces
    // size: num_itf
    std::vector<int> num_tag;

    // stores the local fixed element indices in this part
    std::vector<std::vector<int>> local_fixed_ele;

    // stores the face id of fixed volume element
    // size: num_itf x num_all_fixed_ele[ii]
    std::vector<std::vector<int>> fixed_ele_face_id;

    // stores the volume element's IEN array of the fixed "layer"
    // size: num_itf x (nlocbas x num_all_fixed_ele[ii])
    std::vector<std::vector<int>> fixed_lien;

    // stores the interval tag of fixed volume element
    // size: num_itf x num_all_fixed_ele[ii]
    std::vector<std::vector<int>> fixed_ele_tag;

    std::vector<std::vector<std::vector<int>>> tagged_fixed_ele;

    std::vector<std::vector<int>> num_tagged_fixed_ele;

    // stores the initial coordinates of the nodes from the fixed volume elements
    // size: num_itf x (3 x num_fixed_node[ii])
    std::vector<std::vector<double>> fixed_pt_xyz;

    // the (mapped) global node id corresponding to the fixed_pt_xyz
    // size: num_itf x num_fixed_node[ii]
    // just for debug
    std::vector<std::vector<int>> fixed_node_id;

    // ID array
    // size: num_itf x (num_fixed_node[ii] x dof)
    std::vector<std::vector<int>> fixed_LID;
    
    // stores the local rotated element indices in this part
    std::vector<std::vector<int>> local_rotated_ele;

    // stores the face id of all the rotated volume element
    // size: num_itf x num_all_rotated_ele[ii]
    std::vector<std::vector<int>> rotated_ele_face_id;

    // stores the volume element's IEN array of the rotated "layer"
    // size: num_itf x (nlocbas x num_all_rotated_ele[ii])
    std::vector<std::vector<int>> rotated_lien;

    // stores the interval tag of rotated volume element
    // size: num_itf x num_all_rotated_ele[ii]
    std::vector<std::vector<int>> rotated_ele_tag;

    std::vector<std::vector<std::vector<int>>> tagged_rotated_ele;

    std::vector<std::vector<int>> num_tagged_rotated_ele;

    // stores the initial coordinates of the nodes from the rotated volume elements
    // size: num_itf x (3 x num_rotated_node[ii])
    std::vector<std::vector<double>> rotated_pt_xyz;

    // the (mapped) global node id corresponding to the init_rotated_node_xyz
    // size: num_itf x num_rotated_node[ii]
    // just for debug
    std::vector<std::vector<int>> rotated_node_id;

    // ID array
    // size: num_itf x (num_rotated_node[ii] x dof)
    std::vector<std::vector<int>> rotated_LID;

    ALocal_Interface() = delete;
};
#endif
