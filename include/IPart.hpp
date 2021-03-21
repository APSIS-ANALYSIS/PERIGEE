#ifndef IPART_HPP
#define IPART_HPP
// ==================================================================
// IPart.hpp
// Object:
// An interface for partitioned mesh information. This is a pure 
// abstract interface. Detailed implementation is in the derived classes.
//
// Details:
// 1. elem_loc[nlocalele]: the array storing the global indices of the 
//    elements that belong to this processor.
// 2. node_loc[nlocalnode]: the array storing the global indices of the
//    nodes that belong to this processor.
// 3. ghost_node[nghostnode]: the array that stores the global indices of the
//    ghost nodes that belong to this processor. For i in ghost_nodes,
//    i is in IEN[e][*], e \in elem_loc, but i is not in node_loc.
// 4. nghostnode, ntotnode, nolcalnode, nbadnode: basic int variables
//    that stores the number of different nodes in this processor.
//    nghostnode + nlocalnode = ntotalnode + nbadnode
//    if nbadnode > 0, then the partition is not good.
// 5. LIEN[e][i]: the Local IEN array. e ranges from 0 to nlocalele, 
//    while i ranges from 0 to nlocbas. It gives the local indices of 
//    the corresponding global node.
// 6. local_to_global[nlocalnode]: This is a mapping that maps the local
//    node indices to its global indices. An important relation is that
//     local_to_global[ LIEN[e][i] ] = IEN[ elem_loc[e] ][ i ].
// 7. cpu_rank, cpu_size: MPI parameters telling the number of total
//    processors and the id of this processor.
// 8. bool isdual, int dual_edge_ncommon, bool isMETIS: parameters that
//    deifine the partition.
// 9. nElem, nElem_x, nElem_y, nElem_z, nFunc, nFunc_x, nFunc_y, nFunc_z
//    hx_max hy_max, hz_max, hx_min hy_min hz_min nlocbas: global mesh
//    parameters.
// 10. ctrlPts_x(y,z,w)_loc: local copy of control variables. These 
//     coordniates should have been projected down. 
// 11. hx hy hz: global mesh size of local element in xyz directions.
//
// Date: Sept. 30 2013
// ==================================================================
#include "Sys_Tools.hpp"
#include "Vec_Tools.hpp"
#include "IGlobal_Part.hpp"
#include "Map_Node_Index.hpp"
#include "IIEN.hpp"

class IPart
{
  public:
    IPart(){};
    
    virtual ~IPart(){};

    // --------------------------------------------------------------
    // Generate partition based on the IGlobal_Part class, and the
    // nodal indices are updated based on the Map_Node_Index class.
    // Input:
    // nElem : the number of cells in the global domain
    // nFunc : the number of nodes in the global domain
    // nLocBas : the number of basis functions for the cell
    // cpu_rank : the rank of the cpu for this subdomain 
    // gpart : the global partitioning information
    // mnindex : the map between old and new indices
    // IEN : the global IEN array
    // ctrlPts : the global control points coordinates
    // isPrintinfo : flag for information display. set as true if one
    //                wants to view detailed information for the partitioning
    // Output:
    // elem_loc : the cell indices that belong to this partition
    // nlocalele : the number of local cells (=elem_loc.size)
    // node_loc : the (new) nodal indices that belong to this partition
    // node_loc_original : the old nodal indices
    // node_ghost : the ghost nodes' indices
    // local_to_global : the map from local subdomain nodal indices to the
    //                   global indices
    // nlocalnode : the number of local nodes
    // nghostnode : the number of ghost nodes
    // ntotalnode : nlocalnode + nghostnode - nbadnode
    // nbadnode : the nubmer of bad node, which is the node the local 
    //            subdomain does not need, but belongs to node_loc.
    // nlocghonode : nlocalnode + nghostnode
    // ctrlPts_x_loc : the local and ghost nodes' x coordinate
    // ctrlPts_y_loc : the local and ghost nodes' y coordinate
    // ctrlPts_z_loc : the local and ghost nodes' z coordinate
    // LIEN : the Local IEN array, which is based on the local cell indices
    //        and the local node indices.
    // --------------------------------------------------------------
    virtual void GenPart( const int &nElem,
        const int &nFunc, const int &nLocBas,
        const int &cpu_rank,
        const IGlobal_Part * const &gpart,
        const Map_Node_Index * const &mnindex,
        const IIEN * const &IEN,
        const std::vector<double> &ctrlPts,
        const bool &isPrintinfo,
        std::vector<int> &elem_loc,
        int &nlocalele,
        std::vector<int> &node_loc,
        std::vector<int> &node_loc_original,
        std::vector<int> &node_ghost,
        std::vector<int> &local_to_global,
        int &nlocalnode,
        int &nghostnode,
        int &ntotalnode,
        int &nbadnode,
        int &nlocghonode,
        std::vector<double> &ctrlPts_x_loc,
        std::vector<double> &ctrlPts_y_loc,
        std::vector<double> &ctrlPts_z_loc,
        std::vector<std::vector<int> > &LIEN ) const;

    // Pure virtual function: Write the partition into HDF5 files
    virtual void write( const char * inputFileName ) const = 0;

    // 1. function access element partition.
    virtual int get_elem_loc(int pos) const
    {SYS_T::print_exit("Error: get_elem_loc is not implemented. \n"); return 0;}

    virtual int get_nlocalele() const
    {SYS_T::print_exit("Error: get_nlocalele is not implemented. \n"); return 0;}

    // 2. function access node partition
    virtual int get_node_loc(int pos) const
    {SYS_T::print_exit("Error: get_node_loc is not implemented. \n"); return 0;}

    virtual int get_node_loc_original(int pos) const
    {SYS_T::print_exit("Error: get_node_loc_original is not implemented. \n"); return 0;}

    virtual int get_node_ghost(int pos) const
    {SYS_T::print_exit("Error: get_node_ghost is not implemented. \n"); return 0;}

    virtual int get_local_to_global(int pos) const
    {SYS_T::print_exit("Error: get_local_to_global is not implemented. \n"); return 0;}

    virtual int get_nlocalnode() const
    {SYS_T::print_exit("Error: get_nlocalnode is not implemented. \n"); return 0;}

    virtual int get_nghostnode() const
    {SYS_T::print_exit("Error: get_nghostnode is not implemented. \n"); return 0;}

    virtual int get_ntotalnode() const
    {SYS_T::print_exit("Error: get_ntotalnode is not implemented. \n"); return 0;}

    virtual int get_nbadnode() const
    {SYS_T::print_exit("Error: get_nbadnode is not implemented. \n"); return 0;}

    virtual int get_nlocghonode() const
    {SYS_T::print_exit("Error: get_nlocghonode is not implemented. \n"); return 0;}

    virtual bool isElemInPart(int gloindex) const
    {SYS_T::print_exit("Error: isElemInPart is not implemented. \n"); return false;}

    virtual bool isNodeInPart(int gloindex) const
    {SYS_T::print_exit("Error: isNodeInPart is not implemented. \n"); return false;}

    virtual int get_elemLocIndex(const int &gloindex) const
    {SYS_T::print_exit("Error: get_elemLocIndex is not implemented. \n"); return 0;}


    virtual int get_nodeLocGhoIndex(const int &gloindex) const
    {SYS_T::print_exit("Error: get_nodeLocGhoIndex is not implemented. \n"); return 0;}

    // 3. function that return the partition method parameters
    virtual int get_cpu_rank() const
    {SYS_T::print_exit("Error: get_cpu_rank is not implemented. \n"); return 0;}

    virtual int get_cpu_size() const
    {SYS_T::print_exit("Error: get_cpu_size is not implemented. \n"); return 0;}

    virtual bool get_part_isDual() const
    {SYS_T::print_exit("Error: get_part_isDual is not implemented. \n"); return false;}

    virtual int get_dual_edge_ncommon() const
    {SYS_T::print_exit("Error: get_dual_edge_ncommon is not implemented. \n"); return 0;}

    virtual bool get_isMETIS() const
    {SYS_T::print_exit("Error: get_isMETIS is not implemented. \n"); return false;}

    // 3. Global mesh information
    virtual int get_nElem() const
    {SYS_T::print_exit("Error: get_nElem is not implemented. \n"); return 0;}

    virtual int get_nElem_x() const
    {SYS_T::print_exit("Error: get_nElem_x is not implemented. \n"); return 0;}

    virtual int get_nElem_y() const
    {SYS_T::print_exit("Error: get_nElem_y is not implemented. \n"); return 0;}

    virtual int get_nElem_z() const
    {SYS_T::print_exit("Error: get_nElem_z is not implemented. \n"); return 0;}

    virtual int get_nFunc() const
    {SYS_T::print_exit("Error: get_nFunc is not implemented. \n"); return 0;}

    virtual int get_nFunc_x() const
    {SYS_T::print_exit("Error: get_nFunc_x is not implemented. \n"); return 0;}

    virtual int get_nFunc_y() const
    {SYS_T::print_exit("Error: get_nFunc_y is not implemented. \n"); return 0;}

    virtual int get_nFunc_z() const
    {SYS_T::print_exit("Error: get_nFunc_z is not implemented. \n"); return 0;}

    virtual int get_sDegree() const
    {SYS_T::print_exit("Error: get_sDegree is not implemented. \n"); return 0;}

    virtual int get_tDegree() const
    {SYS_T::print_exit("Error: get_tDegree is not implemented. \n"); return 0;}

    virtual int get_uDegree() const
    {SYS_T::print_exit("Error: get_uDegree is not implemented. \n"); return 0;}

    virtual int get_nLocBas() const
    {SYS_T::print_exit("Error: get_nLocBas is not implemented. \n"); return 0;}

    virtual double get_hx_max() const
    {SYS_T::print_exit("Error: get_hx_max is not implemented. \n"); return 0.0;}

    virtual double get_hy_max() const
    {SYS_T::print_exit("Error: get_hy_max is not implemented. \n"); return 0.0;}

    virtual double get_hz_max() const
    {SYS_T::print_exit("Error: get_hz_max is not implemented. \n"); return 0.0;}

    virtual double get_hx_min() const
    {SYS_T::print_exit("Error: get_hx_min is not implemented. \n"); return 0.0;}

    virtual double get_hy_min() const
    {SYS_T::print_exit("Error: get_hy_min is not implemented. \n"); return 0.0;}

    virtual double get_hz_min() const
    {SYS_T::print_exit("Error: get_hz_min is not implemented. \n"); return 0.0;}

    virtual double get_hx(int ee) const
    {SYS_T::print_exit("Error: get_hx is not implemented. \n"); return 0.0;}

    virtual double get_hy(int ee) const
    {SYS_T::print_exit("Error: get_hy is not implemented. \n"); return 0.0;}

    virtual double get_hz(int ee) const
    {SYS_T::print_exit("Error: get_hz is not implemented. \n"); return 0.0;}

    // 4. LIEN array
    virtual int get_LIEN(int ee, int ii) const
    {SYS_T::print_exit("Error: get_LIEN is not implemented. \n"); return 0;}

    // 5. Control points in this processor
    virtual double get_ctrlPts_x_loc(int pos) const
    {SYS_T::print_exit("Error: get_ctrlPts_x is not implemented. \n"); return 0.0;}

    virtual double get_ctrlPts_y_loc(int pos) const
    {SYS_T::print_exit("Error: get_ctrlPts_y is not implemented. \n"); return 0.0;}

    virtual double get_ctrlPts_z_loc(int pos) const
    {SYS_T::print_exit("Error: get_ctrlPts_z is not implemented. \n"); return 0.0;}

    virtual double get_ctrlPts_w_loc(int pos) const
    {SYS_T::print_exit("Error: get_ctrlPts_w is not implemented. \n"); return 0.0;}


    // 6. print on screen the partition results 
    virtual void print_info() const
    {SYS_T::print_exit("Error: print_info is not implemented. \n");}

    virtual void print_part_ele() const
    {SYS_T::print_exit("Error: print_part_ele is not implemented. \n");}

    virtual void print_part_node() const
    {SYS_T::print_exit("Error: print_part_node is not implemented. \n");}

    virtual void print_part_ghost_node() const
    {SYS_T::print_exit("Error: print_part_ghost_node is not implemented. \n");}

    virtual void print_part_local_to_global() const
    {SYS_T::print_exit("Error: print_part_local_to_global is not implemented. \n");}

    virtual void print_part_LIEN() const
    {SYS_T::print_exit("Error: print_part_LIEN is not implemented. \n");}

    virtual void print_part_loadbalance_edgecut() const
    {SYS_T::print_exit("Error: print_part_loadbalance_edgecut is not implemented. \n");}

    virtual void print_part_bezier_ext_x(int elem) const
    {SYS_T::print_exit("Error: print_part_bezier_ext_x is not implemented. \n");}

    virtual void print_part_bezier_ext_y(int elem) const
    {SYS_T::print_exit("Error: print_part_bezier_ext_y is not implemented. \n");}

    virtual void print_part_bezier_ext_z(int elem) const
    {SYS_T::print_exit("Error: print_part_bezier_ext_z is not implemented. \n");}

};

#endif
