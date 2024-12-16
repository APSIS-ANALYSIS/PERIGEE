#ifndef NODALBC_3D_INFLOW_HPP
#define NODALBC_3D_INFLOW_HPP
// ============================================================================
// NodalBC_3D_inflow.hpp
// 
// This is an instantiation of INodalBC for 3D Inflow type boundary conditions.
//
// For Inflow boundary conditions, 1. there is no periodic type boundary 
// condition; 2. the nodes that belong to the wall are excluded from the 
// dir_nodes list.
//
// Author: Ju Liu
// Date: Aug. 6 2017
// ============================================================================
#include "INodalBC.hpp"
#include "Tet_Tools.hpp"
#include "Hex_Tools.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FEAElement_Quad9_3D_der0.hpp"

class NodalBC_3D_inflow : public INodalBC
{
  public:
    // ------------------------------------------------------------------------
    // Generate the inflow bc given by a list of inffiles.
    // ------------------------------------------------------------------------
    NodalBC_3D_inflow( const std::vector<std::string> &inffileList,
        const std::string &wallfile,
        const int &nFunc,
        const std::vector<Vector_3> &in_outnormal,
        const FEType &elemtype );

    virtual ~NodalBC_3D_inflow() = default;

    virtual unsigned int get_dir_nodes(const unsigned int &ii) const
    {return dir_nodes[ii];}

    virtual unsigned int get_per_slave_nodes(const unsigned int &ii) const
    {
      SYS_T::print_fatal("Error: NodalBC_3D_inflow::get_per_slave_nodes: periodic nodes are not defined.\n");
      return 0;
    }

    virtual unsigned int get_per_master_nodes(const unsigned int &ii) const
    {
      SYS_T::print_fatal("Error: NodalBC_3D_inflow::get_per_master_nodes: periodic nodes are not defined.\n");
      return 0;
    }

    virtual unsigned int get_num_dir_nodes() const {return num_dir_nodes;}

    virtual unsigned int get_num_per_nodes() const {return 0;}

    // get the dirichlet-type nodal index on different nbc_id surfaces
    virtual unsigned int get_dir_nodes_on_inlet( const int &nbc_id, 
        const unsigned int &ii ) const { return dir_nodes_on_inlet[nbc_id][ii]; }

    virtual unsigned int get_num_dir_nodes_on_inlet( const int &nbc_id ) const
    { return num_dir_nodes_on_inlet[nbc_id]; }

    virtual double get_inf_active_area(const int &nbc_id) const {return inf_active_area[nbc_id];}

    virtual Vector_3 get_outnormal(const int &nbc_id) const {return outnormal[nbc_id];}

    // Access to the number of outline boundary points.
    virtual int get_num_out_bc_pts(const int &nbc_id) const {return num_out_bc_pts[nbc_id];}

    // Access to the centroid coordinates.
    virtual Vector_3 get_centroid(const int &nbc_id) const {return centroid[nbc_id];}

    // Access to the outline points. ii ranges from 0 to 3 x num_out_bc_pts[nbc_id];
    virtual double get_outline_pts(const int &nbc_id, const int &ii) const
    {return outline_pts[nbc_id][ii];}

    // Access to the face area
    virtual double get_face_area(const int &nbc_id) const {return face_area[nbc_id];}

    // Access to the integral of NA
    virtual std::vector<double> get_intNA(const int &nbc_id) const {return intNA[nbc_id];}

    // Access to num_nbc
    virtual int get_num_nbc() const {return num_nbc;}

    // Access to num_node
    virtual int get_num_node(const int &nbc_id) const {return num_node[nbc_id];}

    // Access to num_cell
    virtual int get_num_cell(const int &nbc_id) const {return num_cell[nbc_id];}

    // Access to (surface) nLocBas
    virtual int get_nLocBas(const int &nbc_id) const {return nLocBas[nbc_id];}

    // Access to (surface) ien
    virtual int get_ien(const int &nbc_id, const int &cell, const int &lnode) const
    {return sur_ien[nbc_id][ nLocBas[nbc_id] * cell + lnode ];}

    // Access to point coordinates, 
    // node = 0, ..., num_node-1.
    // dir  = 0, 1, 2.
    virtual double get_pt_xyz(const int &nbc_id, const int &node, const int &dir) const
    {return pt_xyz[nbc_id][3*node+dir];}

    // Access to volumetric nodal index
    virtual int get_global_node(const int &nbc_id, const int &node_idx) const
    {return global_node[nbc_id][node_idx];}

    // Access to volumetric cell index
    virtual int get_global_cell(const int &nbc_id, const int &cell_idx) const
    {return global_cell[nbc_id][cell_idx];}

    // For linear element (type Tet4), the face node numbering is
    //   Tet-Face-0 : Node 1 2 3
    //   Tet-Face-1 : Node 0 3 2
    //   Tet-Face-2 : Node 0 1 3
    //   Tet-Face-3 : Node 0 2 1
    // For quadratic element (type Tet10), the face node numbering is
    //   Tet-Face-0 : Node 1 2 3 5 9 8
    //   Tet-Face-1 : Node 0 3 2 7 9 6
    //   Tet-Face-2 : Node 0 1 3 4 8 7
    //   Tet-Face-3 : Node 0 2 1 6 5 4
    // For trilinear element (type Hex8), the face node numbering is
    //   Hex-Face-0 : Node 0 3 2 1
    //   Hex-Face-1 : Node 4 5 6 7
    //   Hex-Face-2 : Node 0 1 5 4
    //   Hex-Face-3 : Node 1 2 6 5
    //   Hex-Face-4 : Node 2 3 7 6
    //   Hex-Face-5 : Node 0 4 7 3
    // For triquadratic element (type Hex27), the face node numbering is
    //   Hex-Face-0 : Node 0 3 2 1 11 10 9 8 24
    //   Hex-Face-1 : Node 4 5 6 7 12 13 14 15 25
    //   Hex-Face-2 : Node 0 1 5 4 8 17 12 16 22
    //   Hex-Face-3 : Node 1 2 6 5 9 18 13 17 21
    //   Hex-Face-4 : Node 2 3 7 6 10 19 14 18 23
    //   Hex-Face-5 : Node 0 4 7 3 16 15 19 11 20
    virtual void resetSurIEN_outwardnormal( const IIEN * const &VIEN );

  private:
    // Disallow default constructor
    NodalBC_3D_inflow() = delete;
    
    // The dirichlet nodes on each inlet surface
    // length num_nbc x num_dir_nodes_on_inlet[ii]
    std::vector< std::vector<unsigned int> > dir_nodes_on_inlet;
    std::vector<unsigned int> num_dir_nodes_on_inlet;

    // All dirichlet nodes stored in a row and the total number of dirichlet
    // nodes
    std::vector<unsigned int> dir_nodes;
    unsigned int num_dir_nodes;

    // number of inlet surfaces and element type
    const int num_nbc;
    const FEType elem_type;
    
    // This is the area calculated by setting the wall nodes to be zero.
    // It is designed to compute the area to give a plug flow profile with
    // a prescribed flow rate. Length num_nbc.
    std::vector<double> inf_active_area;
    
    // This is the actual area of the face. Length num_nbc.
    std::vector<double> face_area;

    // Arithmetic mean of inlet surface points. Length num_nbc.
    std::vector< Vector_3 > centroid;

    // Unit outward normal vector for the surface. Length num_nbc.
    std::vector< Vector_3 > outnormal;

    // Number of boundary points, which are defined as points on the
    // boundary of each inlet surface, i.e., they also belong to
    // the wall mesh. Length num_nbc.
    std::vector<int> num_out_bc_pts;

    // x-y-z coordinates of the outline points.
    // num_nbc times ( 3 x num_out_bc_pts[ii] ) in size.
    std::vector< std::vector<double> > outline_pts;

    // It stores the surface integral of each nodal basis function on
    // the inlet surface.
    // num_nbc x num_node[ii]
    std::vector< std::vector<double> > intNA;

    // number of nodes and cells on each surface. Length num_nbc.
    // Note: num_node[ii] does not equal num_dir_nodes[ii] in this class.
    // nLocBas[ii] is 3, 6, 4 or 9.
    std::vector<int> num_node, num_cell, nLocBas;

    // IEN for each surface.
    // num_nbc times ( nLocBas[ii] x num_cell[ii] ) in size.
    std::vector< std::vector<int> > sur_ien;

    // coordinates of all nodes on each surface
    // num_nbc times ( 3 x num_node[ii] ) in size.
    std::vector< std::vector<double> > pt_xyz;

    // Surface nodes' volumetric mesh ID.
    // num_nbc times num_node[ii] in size.
    std::vector< std::vector<int> > global_node;

    // Surface cell's volumetric mesh ID
    // num_nbc times num_cell[ii] in size.
    std::vector< std::vector<int> > global_cell;

    // ------------------------------------------------------------------------
    // initialization function designed to simplify different constructor
    // implementations.
    virtual void init( const std::vector<std::string> &inffileList,
        const std::string &wallfile,
        const int &nFunc,
        const std::vector<Vector_3> &in_outnormal,
        const FEType &elemtype );

    // Reset function for the IEN array of different element types.
    void reset501IEN_outwardnormal( const IIEN * const &VIEN );

    void reset502IEN_outwardnormal( const IIEN * const &VIEN );

    void reset601IEN_outwardnormal( const IIEN * const &VIEN );

    void reset602IEN_outwardnormal( const IIEN * const &VIEN );
};

#endif
