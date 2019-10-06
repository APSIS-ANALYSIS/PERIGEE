#ifndef BOUNDARY_COND_2D_HPP
#define BOUNDARY_COND_2D_HPP
// ==================================================================
// BoundaryCond2D.hpp
// Object:
// Define the Dirichlet, Periodic slave/master pairs in the boundary;
// Define the element on the boundary with specification of face 
// directions.
//
// Output:
// 1. number of Dirichlet nodes and these nodes' global indices;
// 2. number of Periodic nodes and the master-slave global indices;
// 3. ID array of the essential boundary nodes;
// 4. the global indices of elements that have faces on the top, 
//    bottom, left, right faces.
//
//         top
//     -------------
//   l |           | r
//   e |           | i
//   f |           | g
//   t |           | h
//     |           |
//     -------------  
//         bottom
//
// Date: Feb 17 2014
// ==================================================================
#include <cassert>
#include "IMesh.hpp"
#include "Vec_Tools.hpp"

class BoundaryCond2D
{
  public:
    // constructor
    BoundaryCond2D( const IMesh * const &mesh, const int &bc_type );
    
    // destructor
    virtual ~BoundaryCond2D();

    // Read permission to private data
    int get_dir_nodes(unsigned int ii) const {return dir_nodes[ii];}
    int get_per_slave_nodes(unsigned int ii) const {return per_slave_nodes[ii];}
    int get_per_master_nodes(unsigned int ii) const {return per_master_nodes[ii];}

    unsigned int get_num_dir_nodes() const {return num_dir_nodes;}
    unsigned int get_num_per_nodes() const {return num_per_nodes;}

    int get_ID(unsigned int ii) const {return ID[ii];}
    unsigned int get_num_ID() const {return ID.size();}

    int get_left_elem(unsigned int ii) const {return left_elem[ii];}
    unsigned int get_num_left_elem() const {return num_left_elem;}
    int get_right_elem(unsigned int ii) const {return right_elem[ii];}
    unsigned int get_num_right_elem() const {return num_right_elem;}
    int get_top_elem(unsigned int ii) const {return top_elem[ii];}
    unsigned int get_num_top_elem() const {return num_top_elem;}
    int get_bottom_elem(unsigned int ii) const {return bottom_elem[ii];}
    unsigned int get_num_bottom_elem() const {return num_bottom_elem;}

    // Info print
    void print_info() const;
    void print_edge_cornder_nodes(const IMesh * const &mesh) const;

  private:
    // -------------- DATA ------------------
    // Dirichlet nodes
    std::vector<int> dir_nodes;
    unsigned int num_dir_nodes;

    // Periodic master - slave nodes
    std::vector<int> per_slave_nodes, per_master_nodes;
    unsigned int num_per_nodes;

    // ID array vector
    std::vector<int> ID;

    // Boundary elements in each face
    std::vector<int> left_elem, right_elem, top_elem, bottom_elem;
    unsigned int num_left_elem, num_right_elem, num_top_elem, num_bottom_elem;

    // -------------- FUNCTION ---------------
    // 1. Initialization function for the above data structure
    void clear_nodes();
    void clear_elems();

    // 2. Auxilliary functions for specifying boundary nodes indices
    // --------------------------------------------------------------
    // 2.1 Generate_BCNodes_A: output the global indices of boundary
    //     nodes in the left right top bottom edge, based on the
    //     IMesh class. Four corners of the mesh are also specified.
    //     This generation assumes the tensor product structure of 
    //     the IMesh. 
    //     e.g.
    //     12 13 14 15  |  left: 4, 8
    //      8  9 10 11  |  right: 3, 7
    //      4  5  6  7  |  top: 13 14
    //      0  1  2  3  |  bottom: 1 2
    //                  |  corner: 0 3 12 15
    // --------------------------------------------------------------
    void Generate_BCNodes_A( const IMesh * const &mesh,
        std::vector<int> &left, std::vector<int> &right, std::vector<int> &top, 
        std::vector<int> &bottom, std::vector<int> &corner ) const; 


    // --------------------------------------------------------------
    // 2.2 Generate_BCNodes_B: output the global indices of boundary
    //     nodes up to 2nd level in the left right top bottom. 16 nodes
    //     in the four corners are specified separately. For example
    //
    //    20 21 22 23 24  | lef1: 10; lef2 11
    //    15 16 17 18 19  | rig1: 14; rig2: 13
    //    10 11 12 13 14  | top1: 22; top2 17
    //     5  6  7  8  9  | bot1: 2;  bot2: 7
    //     0  1  2  3  4  | corner_1: 0 4 20 24
    //                    | corner_2: 1 3 21 23
    //                    | corner_3: 5 15 9 19
    //                    | corner_4: 6 8 16 18
    // --------------------------------------------------------------
    void Generate_BCNodes_B( const IMesh * const &mesh,
        std::vector<int> &lef1, std::vector<int> &rig1, std::vector<int> &top1, std::vector<int> &bot1,
        std::vector<int> &lef2, std::vector<int> &rig2, std::vector<int> &top2, std::vector<int> &bot2,
        std::vector<int> &corner1, std::vector<int> &corner2,
        std::vector<int> &corner3, std::vector<int> &corner4 ) const; 


    // --------------------------------------------------------------
    // 2.3 Generate_BCElems_A: output the boundary element with 
    //     differental facial orientations. For an element mesh 
    //     6 7 8  | Bot: 0 1 2
    //     3 4 5  | Top: 6 7 8
    //     0 1 2  | Lef: 0 3 6
    //            | Rig: 2 5 8
    // --------------------------------------------------------------
    void Generate_BCElems_A( const IMesh * const &mesh,
       std::vector<int> &lef, std::vector<int> &rig, 
       std::vector<int> &bot, std::vector<int> &top ) const;


    // 3. Define the boundary nodes and elements
    // --------------------------------------------------------------
    // BC_type_1: all nodes are dirichlet for C0 function, 
    //            no boundary elements
    // --------------------------------------------------------------
    void BC_type_1( const IMesh * const &mesh );

    // --------------------------------------------------------------
    // BC_type_2: \nabla \c \cdot n = 0 strong enforced for C1 element, 
    //            no boundary integral elements
    // --------------------------------------------------------------
    void BC_type_2( const IMesh * const &mesh );

    // --------------------------------------------------------------
    // BC_type_3: strong (i.e. nodal) imposition of C0 periodic bc, 
    //            no boundary integral elements
    // --------------------------------------------------------------
    void BC_type_3( const IMesh * const &mesh );

    // --------------------------------------------------------------
    // BC_type_4: strong (i.e. nodal) imposition of C1 periodic bc, 
    //            no boundary integral elements
    // --------------------------------------------------------------
    void BC_type_4( const IMesh * const &mesh );

    // --------------------------------------------------------------
    // BC_type_5: weak imposition of boundary elements, i.e, all 
    //            element with boundary faces are specified as lef,
    //            rig, top, and bot elements.
    // --------------------------------------------------------------
    void BC_type_5( const IMesh * const &mesh );

    // --------------------------------------------------------------
    // BC_type_6: Weak imposition of nabal c \dot n = 0, i.e., 
    //            adiabatic bc for heat eqn, or "do nothing bc".
    //            no boundary integral elements  
    // --------------------------------------------------------------
    void BC_type_6( const IMesh * const &mesh );

    // --------------------------------------------------------------
    // BC_type_7: strong dirichlet bc at top and bottom, including
    //            four corner nodes;  
    //            adiabatic bc or "do nothing bc" on left and right 
    //            surface
    //            no boundary integral elements  
    // --------------------------------------------------------------
    void BC_type_7( const IMesh * const &mesh );

    // --------------------------------------------------------------
    // BC_type_8: strong dirichlet bc at left and right, excluding
    //            four corner nodes;  
    //            adiabatic bc or "do nothing bc" on top and bottom 
    //            surface
    //            no boundary integral elements  
    // --------------------------------------------------------------
    void BC_type_8( const IMesh * const &mesh );

    // --------------------------------------------------------------
    // BC_type_9: strong dirichlet bc at bottom, including 2 corner nodes
    //            adiabatic bc or "do nothing bc" on Left, Right,
    //            and Top surface with 2 corners
    // --------------------------------------------------------------
    void BC_type_9( const IMesh * const &mesh );

    // 4. Generate ID array according to the dir and per nodes for strong
    //    imposition: ID[dir_nodes] = -1, ID[slave_nodes] = master_nodes
    //    else, ID[ii] = ii;
    void Create_ID(const IMesh * const &mesh);
};

#endif
