#ifndef NODALBC_3D_INFLOW_HPP
#define NODALBC_3D_INFLOW_HPP
// ==================================================================
// NodalBC_3D_inflow.hpp
// 
// This is an instantiation of INodalBC for 3D Inflow type boundary
// conditions.
//
// For Inflow boundary conditions, 1. there is no periodic type
// boundary condition; 2. the nodes that belong to the wall are 
// excluded from the dir_nodes list.
//
// Author: Ju Liu
// Date: Aug. 6 2017
// ==================================================================
#include "INodalBC.hpp"
#include "Tet_Tools.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"

class NodalBC_3D_inflow : public INodalBC
{
  public:
    // --------------------------------------------------------------
    // Generate an empty inflow type boundary condition class
    // --------------------------------------------------------------
    NodalBC_3D_inflow(const int &nFunc);
    
    // --------------------------------------------------------------
    // Generate the inflow bc given by the inffile.
    // --------------------------------------------------------------
    NodalBC_3D_inflow( const std::string &inffile,
        const std::string &wallfile,
        const int &nFunc,
        const Vector_3 &in_outnormal,
        const int &elemtype = 501 );

    // --------------------------------------------------------------
    // Generate the inflow bc given by a list of inffiles.
    // --------------------------------------------------------------
    NodalBC_3D_inflow( const std::vector<std::string> &inffileList,
        const std::string &wallfile,
        const int &nFunc,
        const std::vector<Vector_3> &in_outnormal,
        const int &elemtype = 501 );

    virtual ~NodalBC_3D_inflow() {};

    virtual unsigned int get_dir_nodes(unsigned int &ii) const
    {return dir_nodes[ii];}

    virtual unsigned int get_per_slave_nodes(unsigned int &ii) const
    {return per_slave_nodes[ii];}

    virtual unsigned int get_per_master_nodes(unsigned int &ii) const
    {return per_master_nodes[ii];}

    virtual unsigned int get_num_dir_nodes() const {return num_dir_nodes;}

    virtual unsigned int get_num_per_nodes() const {return num_per_nodes;}

    virtual double get_inf_active_area(const int &nbc_id) const {return inf_active_area[nbc_id];}

    virtual Vector_3 get_para_2(const int &nbc_id) const {return outnormal[nbc_id];}

    // Access to the number of outline boundary points.
    virtual int get_para_3(const int &nbc_id) const {return num_out_bc_pts[nbc_id];}

    // Access to the centroid coordinates.
    virtual Vector_3 get_para_4(const int &nbc_id) const {return centroid[nbc_id];}

    // Access to the outline points. ii ranges from 0 to 3 x num_out_bc_pts[nbc_id];
    virtual double get_para_5(const int &nbc_id, const int &ii) const
    {return outline_pts[nbc_id][ii];}

    // Access to the face area
    virtual double get_para_6(const int &nbc_id) const {return face_area[nbc_id];}

    // Access to the integral of NA
    virtual std::vector<double> get_intNA(const int &nbc_id) const
    {return intNA[nbc_id];}

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
    {return tri_ien[nbc_id][ nLocBas[nbc_id] * cell + lnode ];}

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

  private:
    NodalBC_3D_inflow() : num_nbc(0) {};

    std::vector<unsigned int> dir_nodes;
    unsigned int num_dir_nodes;

    std::vector<unsigned int> per_slave_nodes, per_master_nodes;
    unsigned int num_per_nodes;

    // number of inlet surfaces
    const int num_nbc;
    
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
    // nLocBas[ii] is either 3 or 6.
    std::vector<int> num_node, num_cell, nLocBas;

    // IEN for each surface.
    // num_nbc times ( nLocBas[ii] x num_cell[ii] ) in size.
    std::vector< std::vector<int> > tri_ien;

    // coordinates of all nodes on each surface
    // num_nbc times ( 3 x num_node[ii] ) in size.
    std::vector< std::vector<double> > pt_xyz;

    // Surface nodes' volumetric mesh ID.
    // num_nbc times num_node[ii] in size.
    std::vector< std::vector<int> > global_node;

    // Surface cell's volumetric mesh ID
    // num_nbc times num_cell[ii] in size.
    std::vector< std::vector<int> > global_cell;
};

#endif
