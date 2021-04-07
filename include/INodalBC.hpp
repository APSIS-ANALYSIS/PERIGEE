#ifndef INODALBC_HPP
#define INODALBC_HPP
// ============================================================================
// INodalBC.hpp
// 
// Object:
// This is a pure virtual class that provides an interface for nodal boundary
// condition specification. It provides the global node index for strong
// Dirichlet node; it provides the global node index for strong periodic master
// slave relations. 
//
// Date: Aug 16 2015
// ============================================================================
#include "Sys_Tools.hpp"
#include "Vec_Tools.hpp"

class INodalBC
{
  public:
    INodalBC();

    virtual ~INodalBC();

    // ------------------------------------------------------------------------
    // get_dir_nodes returns the ii-th dirichlet node's global nodal index.
    // The parameter ii runs 0 <= ii < get_num_dir_nodes().
    // ------------------------------------------------------------------------
    virtual unsigned int get_dir_nodes(unsigned int ii) const 
    {return dir_nodes[ii];}

    
    // ------------------------------------------------------------------------
    // get_per_slave_nodes returns the ii-th master-slave pair's slave node's
    // global index. The parameter ii runs as 0 <= ii < get_num_per_nodes().
    // ------------------------------------------------------------------------
    virtual unsigned int get_per_slave_nodes(unsigned int ii) const 
    {return per_slave_nodes[ii];}

    
    // ------------------------------------------------------------------------
    // get_per_slave_nodes returns the ii-th master-slave pair's master node's
    // global index. The parameter ii runs as 0 <= ii < get_num_per_nodes().
    // ------------------------------------------------------------------------
    virtual unsigned int get_per_master_nodes(unsigned int ii) const 
    {return per_master_nodes[ii];}

    
    // ------------------------------------------------------------------------
    // get_num_dir_nodes() gives the total number of Dirichlet nodes that are 
    // strongly enforced.
    // ------------------------------------------------------------------------
    virtual unsigned int get_num_dir_nodes() const {return num_dir_nodes;}


    // ------------------------------------------------------------------------
    // get_num_per_nodes() gives the total number of Periodic(slave) nodes that 
    // are strongly enforced.
    // ------------------------------------------------------------------------
    virtual unsigned int get_num_per_nodes() const {return num_per_nodes;}


    // ------------------------------------------------------------------------
    // get_ID returns the ID-index for the ii-th node. For Dirichlet nodes, its
    // ID-index is -1 (and -1-index will be automatically ignored in FEM
    // assembly); for peroidic slave node, its ID-index is its master's index (
    // in FEM assmelby, the contribution from slave nodes is attributed to the
    // master nodes, and in the slave row, the slave column is 1, the master
    // column is -1, the RHS is 0.) For the rest interior or BC nodes, its ID
    // index is just its global index. 
    // ------------------------------------------------------------------------
    virtual int get_ID(unsigned int ii) const {return ID[ii];}
 

    // ------------------------------------------------------------------------
    // get_num_ID reutrns the length of the ID array.
    // ------------------------------------------------------------------------
    virtual unsigned int get_num_ID() const {return ID.size();}

    
    // ------------------------------------------------------------------------
    // print_info() will print the content of dir_nodes, per_slave_nodes, and
    // per_master_nodes on screen.
    // ------------------------------------------------------------------------
    virtual void print_info() const;

    
    // --------------------------------------------------------------
    // get_para_1() passes additional parameters from the specific
    //              instantiations.
    // --------------------------------------------------------------
    virtual double get_para_1() const
    {
      SYS_T::print_fatal("Error: INodalBC::get_para_1 is not implemented.\n");
      return 0.0;
    }

    
    // --------------------------------------------------------------
    // get_para_2() passes additional parameters from the specific
    //              instantiations.
    // --------------------------------------------------------------
    virtual double get_para_2(const int &ii) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_para_2 is not implemented.\n");
      return 0.0;
    }
  
    // --------------------------------------------------------------
    // get_para_3() passes additional parameters from the specific
    //              instantiations.
    // --------------------------------------------------------------
    virtual int get_para_3() const
    {
      SYS_T::print_fatal("Error: INodalBC::get_para_3 is not implemented.\n");
      return 0.0;
    }

    // --------------------------------------------------------------
    // get_para_4() passes additional parameters from the specific
    //              instantiations.
    // --------------------------------------------------------------
    virtual double get_para_4(const int &ii) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_para_4 is not implemented.\n");
      return 0.0;
    }
    
    // --------------------------------------------------------------
    // get_para_5() passes additional parameters from the specific
    //              instantiations.
    // --------------------------------------------------------------
    virtual double get_para_5(const int &ii) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_para_5 is not implemented.\n");
      return 0.0;
    }
    
    // --------------------------------------------------------------
    // get_para_6() passes additional parameters from the specific
    //              instantiations.
    // --------------------------------------------------------------
    virtual double get_para_6() const
    {
      SYS_T::print_fatal("Error: INodalBC::get_para_6 is not implemented.\n");
      return 0.0;
    }

    // --------------------------------------------------------------
    // get_intNA returns the basis N_A integral on surface
    // --------------------------------------------------------------
    virtual void get_intNA( std::vector<double> &fintNA ) const
    {SYS_T::commPrint("Warning: get_intNA is not implemented. \n");}


    // --------------------------------------------------------------
    // get_num_node returns the number of nodes on surface
    // --------------------------------------------------------------
    virtual int get_num_node() const
    {
      SYS_T::commPrint("Warning: get_num_node is not implemented. \n");
      return -1;
    }
  
    // --------------------------------------------------------------
    // get_num_cell returns the number of cells on surface
    // --------------------------------------------------------------
    virtual int get_num_cell() const
    {
      SYS_T::commPrint("Warning: get_num_cell is not implemented. \n");
      return -1;
    }
  
    // --------------------------------------------------------------
    // get_num_nLocBas returns the number of local basis on surface
    // --------------------------------------------------------------
    virtual int get_nLocBas() const
    {
      SYS_T::commPrint("Warning: get_nLocBas is not implemented. \n");
      return -1;
    }
  
    // --------------------------------------------------------------
    // get_ien returns the surface mesh ien array
    // --------------------------------------------------------------
    virtual int get_ien(const int &cell, const int &lnode) const
    {
      SYS_T::commPrint("Warning: get_ien is not implemented. \n");
      return -1;
    }

    // --------------------------------------------------------------
    // get_pt_xyz returns the points' coordinates
    // --------------------------------------------------------------
    virtual double get_pt_xyz(const int &node, const int &dir) const
    {
      SYS_T::commPrint("Warning: get_pt_xyz is not implemented. \n");
      return -1.0;
    }

    // --------------------------------------------------------------
    // get_global_node returns the node's global volumetric index
    // --------------------------------------------------------------
    virtual int get_global_node(const int &node_idx) const
    {
      SYS_T::commPrint("Warning: get_global_node is not implemented. \n");
      return -1;
    }

    // --------------------------------------------------------------
    // get_global_cell returns the cell's global volumetric index
    // --------------------------------------------------------------
    virtual int get_global_cell(const int &cell_idx) const
    {
      SYS_T::commPrint("Warning: get_global_cell is not implemented. \n");
      return -1;
    } 

    // --------------------------------------------------------------
    // get_cap_id returns the cap id [0, num_caps] of each dir node
    // --------------------------------------------------------------
    virtual void get_cap_id( std::vector<int> &capid ) const
    {SYS_T::commPrint("Warning: get_cap_id is not implemented. \n");}

    // --------------------------------------------------------------
    // get_dominant_comp returns the dominant comp index of each cap's
    // unit normal vector
    // --------------------------------------------------------------
    virtual void get_dominant_comp( std::vector<int> &dom_comp ) const 
    {SYS_T::commPrint("Warning: get_dominant_comp is not implemented. \n");}

    // --------------------------------------------------------------
    // get_outnormal returns each cap's unit normal vector
    // --------------------------------------------------------------
    virtual void get_outnormal( std::vector<double> &outvec ) const
    {SYS_T::commPrint("Warning: get_outnormal is not implemented. \n");}
  
  protected:
    std::vector<unsigned int> dir_nodes;
    unsigned int num_dir_nodes;

    std::vector<unsigned int> per_slave_nodes, per_master_nodes;
    unsigned int num_per_nodes;

    std::vector<int> ID;

    // ------------------------------------------------------------------------
    // Create_ID() will generate the ID array based on dir_nodes and
    // per_xxx_nodes.
    // ------------------------------------------------------------------------
    virtual void Create_ID( const unsigned int &num_node );


    // ------------------------------------------------------------------------
    // Generate_BCNode_2D_A
    // Output the global indices of boundary nodes in the left, right, top,
    // bottom, and corner for a tensor-product two-dimensional mesh.
    //     12 13 14 15  |  left: 4, 8
    //      8  9 10 11  |  right: 3, 7
    //      4  5  6  7  |  top: 13 14
    //      0  1  2  3  |  bottom: 1 2
    //                  |  corner: 0 3 12 15
    // ------------------------------------------------------------------------
    void Generate_BCNode_2D_A( const int &nFunc_x, const int &nFunc_y,
        std::vector<unsigned int> &left, std::vector<unsigned int> &right,
        std::vector<unsigned int> &top, std::vector<unsigned int> &bottom,
        std::vector<unsigned int> &corner ) const; 

    // ------------------------------------------------------------------------
    // Generate_BCNode_2D_B
    // Output the global indices of boundary nodes up to the 2nd layer in the
    // four directions. 16 nodes in the 4 corners are specified. Example
    //    20 21 22 23 24  | lef1: 10; lef2 11
    //    15 16 17 18 19  | rig1: 14; rig2: 13
    //    10 11 12 13 14  | top1: 22; top2 17
    //     5  6  7  8  9  | bot1: 2;  bot2: 7
    //     0  1  2  3  4  | corner_1: 0 4 20 24
    //                    | corner_2: 1 3 21 23
    //                    | corner_3: 5 9 15 19
    //                    | corner_4: 6 8 16 18
    // ------------------------------------------------------------------------
    void Generate_BCNode_2D_B( const int &nFunc_x, const int &nFunc_y,
        std::vector<unsigned int> &lef1, std::vector<unsigned int> &rig1,
        std::vector<unsigned int> &top1, std::vector<unsigned int> &bot1,
        std::vector<unsigned int> &lef2, std::vector<unsigned int> &rig2,
        std::vector<unsigned int> &top2, std::vector<unsigned int> &bot2,
        std::vector<unsigned int> &corner1, std::vector<unsigned int> &corner2,
        std::vector<unsigned int> &corner3, std::vector<unsigned int> &corner4 ) const;



    // -------------------------------------------------------------------------
    // Generate_BCNode_3D_A
    // This may help generate simpler bc specification routine. But this function 
    // is not necessary for generating boundary condition class.
    // In [0,1]^3 cube, front: [1,y,z]; back: [0,y,z]; 
    //                  left: [x,0,z]; right: [x,1,z];
    //                  top:   [x,y,1]; bottom: [x,y,0].  
    // -------------------------------------------------------------------------
    void Generate_BCNode_3D_A( const int &nFunc_x, const int &nFunc_y, 
        const int &nFunc_z,
        std::vector<int> &front, std::vector<int> &back, std::vector<int> &left,
        std::vector<int> &right, std::vector<int> &top, std::vector<int> &bottom ) const;


    // -------------------------------------------------------------------------
    // Generate_BCNode_3D_B
    // Gives more specific specification of boundary nodes:
    // face nodes in front, back, left, right, top, and bottom.
    // edge nodes in edge 01 02 13 23 45 46 57 67 15 37 04 26
    // corner nodes at 0 1 2 3 4 5 6 7
    // -------------------------------------------------------------------------
    void Generate_BCNode_3D_B( const int &nFunc_x, const int &nFunc_y, 
        const int &nFunc_z,
        std::vector<int> &front, std::vector<int> &back, std::vector<int> &left,
        std::vector<int> &right, std::vector<int> &top, std::vector<int> &bottom,
        std::vector<int> &edge01, std::vector<int> &edge02, std::vector<int> &edge13,
        std::vector<int> &edge23, std::vector<int> &edge45, std::vector<int> &edge46,
        std::vector<int> &edge57, std::vector<int> &edge67, std::vector<int> &edge15,
        std::vector<int> &edge37, std::vector<int> &edge04, std::vector<int> &edge26,
        std::vector<int> &corner ) const;


    // -------------------------------------------------------------------------
    // Generate_BCNode_3D_C
    // Gives the specification of interior cube's boundary
    // nodes. Such specification is useful for strong enforcement of C1 type
    // boundary conditions. The date structure is defined in such a way:
    // 1. The six faces are defined similar to the way we defined in _A;
    // 2. The 12 edges, excluding 8 corner points are listed similar to _B;
    // 3. The 8 corner nodes are defined.
    // Refer to the function BC_type_13 for an example of the use of this _C
    // node generation function.
    // -------------------------------------------------------------------------
    void Generate_BCNode_3D_C( const int &nFunc_x, const int &nFunc_y,
        const int &nFunc_z,
        std::vector<int> &front, std::vector<int> &back, std::vector<int> &left,
        std::vector<int> &right, std::vector<int> &top, std::vector<int> &bottom,
        std::vector<int> &edge01, std::vector<int> &edge02, std::vector<int> &edge13,
        std::vector<int> &edge23, std::vector<int> &edge45, std::vector<int> &edge46,
        std::vector<int> &edge57, std::vector<int> &edge67, std::vector<int> &edge15,
        std::vector<int> &edge37, std::vector<int> &edge04, std::vector<int> &edge26,
        std::vector<int> &corner ) const;


    // ------------------------------------------------------------------------
    // Generate_BCNode_3D_D is a fast implementation of the Generate_BCNodes_B
    // function.
    // 0. In this implementation, we added a parameter f_start for multipatch
    // cases. For single patch or the 0-th patch, f_start = 0. For multipatch,
    // f_start = mesh->get_nFunc_start().
    // 
    // 1. One difference is that we will pass the data out by pointers which 
    // allocate dynamic arrays. The users are RESPONSIBLE for delete all these 
    // arrays to free the heap space.
    //
    // 2. To make data local, we put front node and back node together, left
    // node and right node togeter, top and bottom face nodes together. We will
    // pass out a parater to identify the starting index for the second data.
    // back_front: 
    //       0   <= index < nbf   gives the back face nodes;
    //       nbf <= index < 2*nbf gives the front face nodes;
    // left_right:
    //       0   <= index < nlr   gives the left face nodes;
    //       nlr <= index < 2*nlr gives the right face nodes;
    // bottom_top:
    //       0   <= index < nbt   gives the bottom face nodes;
    //       nbt <= index < 2*nbt gives the top face nodes.
    // edge01 stores the edge 01 23 45 67 nodes:
    //       0   <= ii < nex   gives edge 01;
    //       nex <= ii < 2*nex gives edge 23;
    //     2*nex <= ii < 3*nex gives edge 45;
    //     3*nex <= ii < 4*nex gives edge 67.
    // edge02 stores the edge 02 13 46 57 nodes:
    //       0   <= ii < ney   gives edge 02;
    //       ney <= ii < 2*ney gives edge 13;
    //     2*ney <= ii < 3*ney gives edge 46;
    //     3*ney <= ii < 4*ney gives edge 57.
    // edge04 stores the edge 04 15 26 37 nodes:
    //       0   <= ii < nez   gives edge 04;
    //       nez <= ii < 2*nez gives edge 15;
    //     2*nez <= ii < 3*nez gives edge 26;
    //     3*nez <= ii < 4*nez gives edge 37.
    // corner pointer stores 8 nodes from 0 to 7.
    // ------------------------------------------------------------------------
    void Generate_BCNode_3D_D( const int &nFunc_x, const int &nFunc_y,
        const int &nFunc_z, const int &f_start,
        int * &back_front, int &nbf,
        int * &left_right, int &nlr,
        int * &bottom_top, int &nbt,
        int * &edge01, int &nex,
        int * &edge02, int &ney,
        int * &edge04, int &nez,
        int * &corner ) const;

};

#endif
