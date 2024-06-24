#ifndef INTERFACE_PAIR_HPP
#define INTERFACE_PAIR_HPP
// ============================================================================
// Interface_pair.hpp
// 
// A combination of two surface mesh information of an interface.
// Assume there is a fixed surface and a rotated surface.
// The numerical quadrature points will be generated on the fixed surface,
// and the counterpart will be searched on the rotated surface.
// For the efficiency of searching, the following data is critical:
// 
// 1. Partition tag of the fixed element.
//   The fixed elements are uniformly distributed to each cpu,
//   independent from the global partition the volumes.
// 2. Interval tag of the fixed and rotated element.
//   They will be grouped in to several intervals according to their initial
//   position. The rotated elements that have the same interval tag with a fixed
//   one will obtain a priority in the searching.
//
// The other data constructions are similar with ElemBC_turbulence_wall_model.
// 
// Author: Xuanming Huag
// Date: Jun 24 2024
// ============================================================================

#include "Sys_Tools.hpp"
#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "VTK_Tools.hpp"
#include "HDF5_Tools.hpp"
#include "IIEN.hpp"
#include "Tet_Tools.hpp"
#include "Hex_Tools.hpp"

class Interface_pair
{
  public:
    // Constructor for type 0 interface-pair, with a direction(0, 1, 2) input
    Interface_pair( const std::string &fixed_vtkfile,
                    const std::string &rotated_vtkfile,
                    const std::string &fixed_h5file,
                    const int &total_num_fixed_elem,
                    const int &total_num_fixed_pt,
                    const std::vector<double> &all_vol_ctrlPts,
                    const IIEN * const &VIEN,
                    const int &elemtype_in,
                    const std::vector<double> &intervals_in,
                    const int &direction_in );

    // Constructor for type 1 interface-pair, with the centroid of the surface input
    Interface_pair( const std::string &fixed_vtkfile,
                    const std::string &rotated_vtkfile,
                    const std::string &fixed_h5file,
                    const int &total_num_fixed_elem,
                    const int &total_num_fixed_pt,
                    const std::vector<double> &all_vol_ctrlPts,
                    const IIEN * const &VIEN,
                    const int &elemtype_in,
                    const std::vector<double> &intervals_in,
                    const Vector_3 &centroid_in );

    virtual int get_num_fixed_ele() const
    {return  num_fixed_ele;}

    virtual int get_fixed_part_tag(const int &cell_index) const
    {return fixed_part_tag[cell_index];}

    virtual int get_fixed_faceID(const int &cell_index) const
    {return fixed_face_id[cell_index];}

    virtual std::vector<int> get_FL_vien() const
    {return fixed_layer_vien;}

    virtual std::vector<int> get_FLN_GID() const
    {return fixed_layer_global_node;}

    virtual std::vector<double> get_FLN_xyz() const
    {return fixed_layer_pt_xyz;}

    virtual std::vector<int> get_FIT() const
    {return fixed_interval_tag;}

    virtual std::vector<int> get_rotated_faceID() const
    {return rotated_face_id;}

    virtual std::vector<int> get_RL_vien() const
    {return rotated_layer_vien;}

    virtual std::vector<int> get_RLN_GID() const
    {return rotated_layer_global_node;}

    virtual std::vector<double> get_RLN_xyz() const
    {return rotated_layer_pt_xyz;}

    virtual std::vector<int> get_RIT() const
    {return rotated_interval_tag;}

    virtual ~Interface_pair(){};

  private:
    // 0: Lofted along an axis
    // 1: Top or bottom
    const int interface_type;

    // the number of local basis function in a surface/volume element
    int s_nLocBas;
    int v_nLocBas;

    // the axial direction of type 0 interface-pair
    const int T0_axial_direction;

    // the centroid of the type 1 interface-pair
    const Vector_3 T1_surface_centroid;

    // the number of the fixed layer elements
    int num_fixed_ele;

    // the partition tag of the fixed layer elements
    std::vector<int> fixed_part_tag;

    // the face id of the fixed layer elements
    std::vector<int> fixed_face_id;

    // the ien array of the fixed layer
    std::vector<int> fixed_layer_vien;

    // the GlobalNodeID of the fixed layer nodes
    std::vector<int> fixed_layer_global_node;

    // the xyz-coordinate of nodes, corresponding to the fixed_layer_global_node
    std::vector<double> fixed_layer_pt_xyz;

    // the interval tag of the fixed layer elements
    std::vector<int> fixed_interval_tag;

    // the number of the rotated layer elements
    int num_rotated_ele;

    // the face id of the rotated layer elements
    std::vector<int> rotated_face_id;

    // the ien array of the rotated layer
    std::vector<int> rotated_layer_vien;

    // the GlobalNodeID of the rotated layer nodes
    std::vector<int> rotated_layer_global_node;

    // the xyz-coordinate of nodes, corresponding to the rotated_layer_global_node
    std::vector<double> rotated_layer_pt_xyz;

    // the interval tag of the rotated layer elements
    std::vector<int> rotated_interval_tag;

    virtual void Initialize(const std::string &fixed_vtkfile,
      const std::string &rotated_vtkfile,
      const std::string &fixed_h5file,
      const int &total_num_fixed_elem,
      const int &total_num_fixed_pt,
      const std::vector<double> &all_vol_ctrlPts,
      const IIEN * const &VIEN,
      const int &elemtype_in,
      const std::vector<double> &intervals_in);

    // Check the validity of the input interval
    virtual void Check_interval(const std::vector<double> &intervals);

    // Attach the interval tags to all the elements
    virtual void Tag(const std::vector<double> &intervals,
      const int &num_fixed_ele,
      const std::vector<double> &fixed_sur_pt_xyz,
      const std::vector<int> &fixed_sur_ien,
      const int &num_rotated_ele,
      const std::vector<double> &rotated_sur_pt_xyz,
      const std::vector<int> &rotated_sur_ien);

    // Deside the interval tag of a single element
    virtual int Group(const Vector_3 &ele_centroid, const std::vector<double> &intervals);

    Interface_pair() = delete;

};

#endif