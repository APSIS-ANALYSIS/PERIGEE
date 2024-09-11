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
#include "IIEN.hpp"

class INodalBC
{
  public:
    INodalBC() = default;

    virtual ~INodalBC() = default;

    // ------------------------------------------------------------------------
    // get_dir_nodes returns the ii-th dirichlet node's global nodal index.
    // The parameter ii runs 0 <= ii < get_num_dir_nodes(nbc_id).
    // ------------------------------------------------------------------------
    virtual unsigned int get_dir_nodes(const unsigned int &ii) const = 0;
    
    // ------------------------------------------------------------------------
    // get_per_slave_nodes returns the ii-th master-slave pair's slave node's
    // global index. The parameter ii runs as 0 <= ii < get_num_per_nodes(nbc_id).
    // ------------------------------------------------------------------------
    virtual unsigned int get_per_slave_nodes(const unsigned int &ii) const = 0;
    
    // ------------------------------------------------------------------------
    // get_per_slave_nodes returns the ii-th master-slave pair's master node's
    // global index. The parameter ii runs as 0 <= ii < get_num_per_nodes(nbc_id).
    // ------------------------------------------------------------------------
    virtual unsigned int get_per_master_nodes(const unsigned int &ii) const = 0;
    
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
    virtual int get_ID(const unsigned int &ii) const {return ID[ii];}
 
    // ------------------------------------------------------------------------
    // print_info() will print the content of dir_nodes, per_slave_nodes, and
    // per_master_nodes on screen.
    // ------------------------------------------------------------------------
    virtual void print_info() const
    {
      std::cout<<std::endl;
      std::cout<<"======== BC info ======="<<std::endl;
      if( get_num_dir_nodes() > 0 )
      {
        std::cout<<"Dirichlet nodes: "<<std::endl;
        for(unsigned int ii=0; ii<get_num_dir_nodes(); ++ii)
          std::cout<<get_dir_nodes(ii)<<'\t';
        std::cout<<std::endl;
      }

      if( get_num_per_nodes() > 0 )
      {
        std::cout<<"Periodic master - slave nodes: "<<std::endl;
        for(unsigned int ii=0; ii<get_num_per_nodes(); ++ii)
          std::cout<<get_per_master_nodes(ii)<<'\t'<<get_per_slave_nodes(ii)<<std::endl;
      }

      std::cout<<std::endl<<"ID array: "<<std::endl;
      for(unsigned int ii=0; ii<ID.size(); ++ii)
        std::cout<<ID[ii]<<'\t';
      std::cout<<'\n';
      std::cout<<std::endl<<"========================"<<std::endl;
    }

    // ------------------------------------------------------------------------
    // get_inf_active_area() return inflow cap surface active area
    // ------------------------------------------------------------------------
    virtual double get_inf_active_area(const int &nbc_id) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_inf_active_area is not implemented.\n");
      return 0.0;
    }

    // ------------------------------------------------------------------------
    // get_dir_nodes_on_inlet return the Dirichlet node index on each inlet
    // surface.
    // ------------------------------------------------------------------------
    virtual unsigned int get_dir_nodes_on_inlet( const int &nbc_id, const unsigned int &ii ) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_dir_nodes_on_inlet is not implemented.\n");
      return 0;
    }

    // ------------------------------------------------------------------------
    // get_num_dir_nodes_on_inlet return the number of Dirichlet nodes on each 
    // inlet surface.
    // ------------------------------------------------------------------------
    virtual unsigned int get_num_dir_nodes_on_inlet( const int &nbc_id ) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_num_dir_nodes_on_inlet is not implemented.\n");
      return 0;
    }

    // ------------------------------------------------------------------------
    // get_outnormal() return inflow cap surface outward normal vector
    // ------------------------------------------------------------------------
    virtual Vector_3 get_outnormal(const int &nbc_id) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_outnormal is not implemented.\n");
      return Vector_3();
    }

    // ------------------------------------------------------------------------
    // get_num_out_bc_pts() return the number of outline boundary points
    // ------------------------------------------------------------------------
    virtual int get_num_out_bc_pts(const int &nbc_id) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_num_out_bc_pts is not implemented.\n");
      return -1;
    }

    // ------------------------------------------------------------------------
    // get_centroid() return inflow cap surface centroid coordinates
    // ------------------------------------------------------------------------
    virtual Vector_3 get_centroid(const int &nbc_id) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_centroid is not implemented.\n");
      return Vector_3();
    }

    // ------------------------------------------------------------------------
    // get_outline_pts() return the outline points. ii ranges from 0 to 3 x
    // num_out_bc_pts[nbc_id]
    // ------------------------------------------------------------------------
    virtual double get_outline_pts(const int &nbc_id, const int &ii) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_outline_pts is not implemented.\n");
      return 0.0;
    }

    // ------------------------------------------------------------------------
    // get_face_area() return inflow cap surface face area
    // ------------------------------------------------------------------------
    virtual double get_face_area(const int &nbc_id) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_face_area is not implemented.\n");
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

    // --------------------------------------------------------------
    // Reset the triangle element's surface IEN so that the outward normal
    // vector is defined.
    // --------------------------------------------------------------
    virtual void resetSurIEN_outwardnormal( const IIEN * const &VIEN )
    {SYS_T::print_fatal("Warning: resetSurIEN_outwardnormal is not implemented. \n");}

    // --------------------------------------------------------------
    // get the dirichlet-type nodal index on different nbc_id surfaces
    // --------------------------------------------------------------
    virtual unsigned int get_dir_nodes_on_rotated_surface( const unsigned int &ii ) const
    {
      SYS_T::print_fatal("Error: INodalBC::get_dir_nodes_on_rotated_surface is not implemented. \n");
      return 0;
    }

    virtual unsigned int get_num_dir_nodes_on_rotated_surface() const
    {
      SYS_T::print_fatal("Error: INodalBC::get_num_dir_nodes_on_rotated_surface is not implemented. \n");
      return 0;
    }   

    // Access to num_node
    virtual int get_num_node() const
    {
      SYS_T::commPrint("Warning: get_num_node is not implemented. \n");
      return -1;
    }

    // Access to num_cell
    virtual int get_num_cell() const
    {
      SYS_T::commPrint("Warning: get_num_cell is not implemented. \n");
      return -1;
    }

    // Access to (surface) nLocBas
    virtual int get_nLocBas() const
    {
      SYS_T::commPrint("Warning: get_nLocBas is not implemented. \n");
      return -1;
    }

    // Access to (surface) ien
    virtual int get_ien(const int &cell, const int &lnode) const
    {
      SYS_T::commPrint("Warning: get_nLocBas is not implemented. \n");
      return -1;
    }

    // Access to point coordinates, 
    // node = 0, ..., num_node-1.
    // dir  = 0, 1, 2.
    virtual double get_pt_xyz(const int &node, const int &dir) const
    {
      SYS_T::commPrint("Warning: get_nLocBas is not implemented. \n");
      return -1.0;
    }

    // Access to volumetric nodal index
    virtual int get_global_node(const int &node_idx) const
    {
      SYS_T::commPrint("Warning: get_nLocBas is not implemented. \n");
      return -1;
    }

    // Access to volumetric cell index
    virtual int get_global_cell(const int &cell_idx) const
    {
      SYS_T::commPrint("Warning: get_nLocBas is not implemented. \n");
      return -1;
    }

  protected:
    std::vector<int> ID {};

    // ------------------------------------------------------------------------
    // Create_ID() will generate the ID array based on get_dir_nodes and
    // get_per_xxx_nodes.
    // ------------------------------------------------------------------------
    virtual void Create_ID( const unsigned int &num_node )
    {
      ID.resize(num_node); VEC_T::shrink2fit(ID);

      for(unsigned int ii = 0; ii<ID.size(); ++ii) ID[ii] = ii;

      for(unsigned int ii = 0; ii<get_num_dir_nodes(); ++ii)
        ID[ get_dir_nodes(ii) ] = -1;

      for(unsigned int ii = 0; ii<get_num_per_nodes(); ++ii)
        ID[ get_per_slave_nodes(ii) ] = get_per_master_nodes(ii);
    }
};

#endif
