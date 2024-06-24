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
#include "HDF5_Reader.hpp"
#include "Math_Tools.hpp"

class ALocal_Interface
{
    public:
        ALocal_Interface( const std::string &fileBaseName, const int &cpu_rank );

        virtual ~ALocal_Interface() = default;

        virtual void print_info() const;

        virtual int get_num_itf() const
        {return num_itf;}

        virtual int get_num_fixed_ele(const int &ii) const
        {return num_fixed_ele[ii];}

        virtual int get_num_tag(const int &ii) const
        {return num_tag[ii];}

        virtual int get_num_rotated_ele(const int &ii, const int &tag) const
        {return num_rotated_ele[ii][tag];}

        virtual int get_num_rotated_node(const int &ii) const
        {return num_rotated_node[ii];}

        virtual int get_fixed_face_id(const int &ii, const int &jj) const
        {return fixed_ele_face_id[ii][jj];}

        virtual int get_fixed_ele_tag(const int &ii, const int &jj) const
        {return fixed_ele_tag[ii][jj];}

        virtual int get_fixed_layer_ien(const int &ii, const int &jj) const
        {return fixed_layer_ien[ii][jj];}

        virtual double get_fixed_node_xyz(const int &ii, const int &jj) const
        {return fixed_node_xyz[ii][jj];}

        virtual int get_fixed_node_id(const int &ii, const int &jj) const
        {return fixed_node_id[ii][jj];}

        virtual int get_rotated_layer_ien(const int &ii, const int &tag, const int &jj) const
        {return rotated_layer_ien[ii][tag][jj];}

        virtual int get_rotated_face_id(const int &ii, const int &tag, const int &jj) const
        {return rotated_layer_face_id[ii][tag][jj];}

        virtual double get_init_rotated_node_xyz(const int &ii, const int &jj) const
        {return init_rotated_node_xyz[ii][jj];}

        virtual int get_rotated_node_id(const int &ii, const int &jj) const
        {return rotated_node_id[ii][jj];}

        // Get the current rotated nodes' xyz with given rotation rule and time
        virtual Vector_3 get_curr_xyz(const int &ii, const int &node, const double &tt) const;

        virtual void get_fixed_ele_ctrlPts(const int &ii, const int &ee,
            double * const volctrl_x,  double * const volctrl_y,  double * const volctrl_z) const;

        virtual void get_rotated_ele_ctrlPts(const int &ii, const int &tag,const int &ee, const double &tt,
            double * const volctrl_x,  double * const volctrl_y,  double * const volctrl_z) const;

    protected:
        // the number of interfaces
        int num_itf;

        // the number of local basis function
        int nLocBas;

        // the number of fixed volume elements in this part
        // size: num_itf
        std::vector<int> num_fixed_ele;

        // the number of rotated volume elements of each interface
        // size: num_itf x num_tag[ii]
        std::vector<std::vector<int>> num_rotated_ele;

        // the number of the nodes from the fixed volume elements
        // size: num_itf
        std::vector<int> num_fixed_node;

        // the number of the nodes from the rotated volume elements
        // size: num_itf
        std::vector<int> num_rotated_node;

        // the number of interval tag of each pair of interfaces
        // size: num_itf
        std::vector<int> num_tag;

        // stores the face id of fixed volume element
        // size: num_itf x num_fixed_ele[ii]
        std::vector<std::vector<int>> fixed_ele_face_id;

        std::vector<std::vector<int>> fixed_ele_tag;

        // stores the volume element's IEN array of the fixed "layer"
        // size: num_itf x (nlocbas x num_fixed_ele[ii])
        std::vector<std::vector<int>> fixed_layer_ien;

        // stores the initial coordinates of the nodes from the fixed volume elements
        // size: num_itf x (3 x num_fixed_node[ii])
        std::vector<std::vector<double>> fixed_node_xyz;

        // the (mapped) global node id corresponding to the fixed_node_xyz
        // size: num_itf x num_fixed_node[ii]
        std::vector<std::vector<int>> fixed_node_id;
        
        // stores the volume element's IEN array of the rotated "layer"
        // size: num_itf x (nlocbas x num_rotated_ele[ii])
        std::vector<std::vector<std::vector<int>>> rotated_layer_ien;
    
        // stores the face id of all the rotated volume element
        // size: num_itf x num_rotated_ele[ii]
        std::vector<std::vector<std::vector<int>>> rotated_layer_face_id;

        // stores the initial coordinates of the nodes from the rotated volume elements
        // size: num_itf x (3 x num_rotated_node[ii])
        std::vector<std::vector<double>> init_rotated_node_xyz;

        // the (mapped) global node id corresponding to the init_rotated_node_xyz
        // size: num_itf x num_rotated_node[ii]
        std::vector<std::vector<int>> rotated_node_id;

        ALocal_Interface() = delete;
};
#endif
