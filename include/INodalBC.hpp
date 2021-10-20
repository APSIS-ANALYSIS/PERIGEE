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
#include "Vector_3.hpp"

class INodalBC
{
  public:
    INodalBC();

    virtual ~INodalBC();

    // ------------------------------------------------------------------------
    // get_dir_nodes returns the ii-th dirichlet node's global nodal index.
    // The parameter ii runs 0 <= ii < get_num_dir_nodes(nbc_id).
    // ------------------------------------------------------------------------
    virtual unsigned int get_dir_nodes(unsigned int &ii) const = 0;
    
    // ------------------------------------------------------------------------
    // get_per_slave_nodes returns the ii-th master-slave pair's slave node's
    // global index. The parameter ii runs as 0 <= ii < get_num_per_nodes(nbc_id).
    // ------------------------------------------------------------------------
    virtual unsigned int get_per_slave_nodes(unsigned int &ii) const = 0;
    
    // ------------------------------------------------------------------------
    // get_per_slave_nodes returns the ii-th master-slave pair's master node's
    // global index. The parameter ii runs as 0 <= ii < get_num_per_nodes(nbc_id).
    // ------------------------------------------------------------------------
    virtual unsigned int get_per_master_nodes(unsigned int &ii) const = 0;
    
    // ------------------------------------------------------------------------
    // get_num_dir_nodes() gives the total number of Dirichlet nodes that are 
    // strongly enforced.
    // ------------------------------------------------------------------------
    virtual unsigned int get_num_dir_nodes() const = 0;
    
    // ------------------------------------------------------------------------
    // get_num_per_nodes() gives the total number of Periodic(slave) nodes that 
    // are strongly enforced.
    // ------------------------------------------------------------------------
    virtual unsigned int get_num_per_nodes() const = 0;
    
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
    virtual double get_para_1(const int &nbc_id) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_para_1 is not implemented.\n");
      return 0.0;
    }

    // --------------------------------------------------------------
    // get_para_2() passes additional parameters from the specific
    //              instantiations.
    // --------------------------------------------------------------
    virtual Vector_3 get_para_2(const int &nbc_id) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_para_2 is not implemented.\n");
      return Vector_3();
    }
  
    // --------------------------------------------------------------
    // get_para_3() passes additional parameters from the specific
    //              instantiations.
    // --------------------------------------------------------------
    virtual int get_para_3(const int &nbc_id) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_para_3 is not implemented.\n");
      return -1;
    }

    // --------------------------------------------------------------
    // get_para_4() passes additional parameters from the specific
    //              instantiations.
    // --------------------------------------------------------------
    virtual Vector_3 get_para_4(const int &nbc_id) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_para_4 is not implemented.\n");
      return Vector_3();
    }
    
    // --------------------------------------------------------------
    // get_para_5() passes additional parameters from the specific
    //              instantiations.
    // --------------------------------------------------------------
    virtual double get_para_5(const int &nbc_id, const int &ii) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_para_5 is not implemented.\n");
      return 0.0;
    }
    
    // --------------------------------------------------------------
    // get_para_6() passes additional parameters from the specific
    //              instantiations.
    // --------------------------------------------------------------
    virtual double get_para_6(const int &nbc_id) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_para_6 is not implemented.\n");
      return 0.0;
    }

    // --------------------------------------------------------------
    // get_num_nbc() returns the number of nodal bc sets 
    // --------------------------------------------------------------
    virtual int get_num_nbc() const
    {
      SYS_T::commPrint("Warning: get_num_nbc is not implemented. \n");
      return -1;
    }

    // --------------------------------------------------------------
    // get_intNA returns the basis N_A integral on surface
    // --------------------------------------------------------------
    virtual std::vector<double> get_intNA(const int &nbc_id) const
    {SYS_T::commPrint("Warning: get_intNA is not implemented. \n"); return {};}

    // --------------------------------------------------------------
    // get_num_node returns the number of nodes on surface
    // --------------------------------------------------------------
    virtual int get_num_node(const int &nbc_id) const
    {
      SYS_T::commPrint("Warning: get_num_node is not implemented. \n");
      return -1;
    }
  
    // --------------------------------------------------------------
    // get_num_cell returns the number of cells on surface
    // --------------------------------------------------------------
    virtual int get_num_cell(const int &nbc_id) const
    {
      SYS_T::commPrint("Warning: get_num_cell is not implemented. \n");
      return -1;
    }
  
    // --------------------------------------------------------------
    // get_num_nLocBas returns the number of local basis on surface
    // --------------------------------------------------------------
    virtual int get_nLocBas(const int &nbc_id) const
    {
      SYS_T::commPrint("Warning: get_nLocBas is not implemented. \n");
      return -1;
    }
  
    // --------------------------------------------------------------
    // get_ien returns the surface mesh ien array
    // --------------------------------------------------------------
    virtual int get_ien(const int &nbc_id, const int &cell, const int &lnode) const
    {
      SYS_T::commPrint("Warning: get_ien is not implemented. \n");
      return -1;
    }

    // --------------------------------------------------------------
    // get_pt_xyz returns the points' coordinates
    // --------------------------------------------------------------
    virtual double get_pt_xyz(const int &nbc_id, const int &node, const int &dir) const
    {
      SYS_T::commPrint("Warning: get_pt_xyz is not implemented. \n");
      return -1.0;
    }

    // --------------------------------------------------------------
    // get_global_node returns the node's global volumetric index
    // --------------------------------------------------------------
    virtual int get_global_node(const int &nbc_id, const int &node_idx) const
    {
      SYS_T::commPrint("Warning: get_global_node is not implemented. \n");
      return -1;
    }

    // --------------------------------------------------------------
    // get_global_cell returns the cell's global volumetric index
    // --------------------------------------------------------------
    virtual int get_global_cell(const int &nbc_id, const int &cell_idx) const
    {
      SYS_T::commPrint("Warning: get_global_cell is not implemented. \n");
      return -1;
    } 

    // --------------------------------------------------------------
    // get_num_caps returns the number of cap surfaces
    // --------------------------------------------------------------
    virtual int get_num_caps() const
    {
      SYS_T::commPrint("Warning: get_num_caps in not implemented\n");
      return -1;
    }
    
    // --------------------------------------------------------------
    // get_ring_bc_type returns essential bc type for ring nodes
    // --------------------------------------------------------------
    virtual int get_ring_bc_type() const
    {
      SYS_T::commPrint("Warning: get_ring_bc_type in not implemented\n");
      return -1;
    }

    // --------------------------------------------------------------
    // get_cap_id returns the cap id [0, num_caps] of each dir node
    // --------------------------------------------------------------
    virtual std::vector<int> get_cap_id() const
    {SYS_T::commPrint("Warning: get_cap_id is not implemented.\n"); return {};}

    // --------------------------------------------------------------
    // get_outnormal returns each cap's unit normal vector
    // --------------------------------------------------------------
    virtual std::vector<double> get_outnormal() const
    {SYS_T::commPrint("Warning: get_outnormal is not implemented.\n"); return {};}
  
    // --------------------------------------------------------------
    // get_rotation_matrix returns each cap's 3x3 rotation matrix for
    //                     skew boundary conditions
    // --------------------------------------------------------------
    virtual std::vector<double> get_rotation_matrix() const
    {SYS_T::commPrint("Warning: get_rotation_matrix is not implemented. \n"); return {};}
  
  protected:

    // num_nbc x num_dir_nodes[ii] in dimension
    //std::vector< std::vector<unsigned int> > dir_nodes;
    //std::vector<unsigned int> num_dir_nodes;

    // num_nbc x num_per_nodes[ii] in dimension
    //std::vector< std::vector<unsigned int> > per_slave_nodes, per_master_nodes;
    //std::vector<unsigned int> num_per_nodes;

    std::vector<int> ID;

    // ------------------------------------------------------------------------
    // Create_ID() will generate the ID array based on get_dir_nodes and
    // get_per_xxx_nodes.
    // ------------------------------------------------------------------------
    virtual void Create_ID( const unsigned int &num_node );
};

#endif
