#ifndef IPART_HPP
#define IPART_HPP
// ============================================================================
// IPart.hpp
// A pure abstract interface for mesh partitioning.
//
// Details:
// 1. elem_loc[nlocalele]: the array stores the global indices of the  elements
//    that belong to this processor/subdomain.
// 2. node_loc[nlocalnode]: the array stores the global indices of the nodes
//    that belong to this processor/subdomain.
// 3. ghost_node[nghostnode]: the array that stores the global indices of the
//    ghost nodes that belong to this processor. For node in ghost_nodes,
//    node is in {IEN[e][*]}, for e in elem_loc, but node is NOT in node_loc.
// 4. nghostnode, ntotnode, nolcalnode, nbadnode: basic integere variables
//    that stores the number of different nodes in this processor. They satisfy
//    the following relation:
//            nghostnode + nlocalnode = ntotalnode + nbadnode
//    if nbadnode > 0, then the partition is typically not accepted.
// 5. LIEN[e][i]: the Local IEN array. e ranges from 0 to nlocalele, 
//    while i ranges from 0 to nlocbas. It gives the local (local to the
//    subdomain) indices of the corresponding global node.
// 6. local_to_global[nlocalnode]: This is a mapping that maps the local
//    node indices to its global indices. It is defined as follows,
//            local_to_global[ LIEN[e][i] ] = IEN[ elem_loc[e] ][ i ].
// 7. cpu_rank, cpu_size: MPI parameters telling the number of total processors
//    and the id of this processor.
// 8. dual_edge_ncommon: parameters that deifine the partition.
// 9. nElem, nFunc, nlocbas: global mesh parameters.
// 10. ctrlPts_x(y,z,w)_loc: local copy of control variables. 
//
// Date: Sept. 30 2013
// ============================================================================
#include "Sys_Tools.hpp"

class IPart
{
  public:
    IPart() = default;
    
    virtual ~IPart() = default;

    // Pure virtual function: Write the partition into HDF5 files
    virtual void write( const std::string &inputFileName ) const = 0;

    // 1. function access element partition.
    virtual int get_elem_loc(const int &pos) const
    {SYS_T::print_fatal("Error: get_elem_loc is not implemented. \n"); return 0;}

    virtual int get_nlocalele() const
    {SYS_T::print_fatal("Error: get_nlocalele is not implemented. \n"); return 0;}

    // 2. function access node partition
    virtual int get_node_loc(const int &pos) const
    {SYS_T::print_fatal("Error: get_node_loc is not implemented. \n"); return 0;}

    virtual int get_node_loc_original(const int &pos) const
    {SYS_T::print_fatal("Error: get_node_loc_original is not implemented. \n"); return 0;}

    virtual int get_node_ghost(const int &pos) const
    {SYS_T::print_fatal("Error: get_node_ghost is not implemented. \n"); return 0;}

    virtual int get_local_to_global(const int &pos) const
    {SYS_T::print_fatal("Error: get_local_to_global is not implemented. \n"); return 0;}

    virtual int get_nlocalnode() const
    {SYS_T::print_fatal("Error: get_nlocalnode is not implemented. \n"); return 0;}

    virtual int get_nghostnode() const
    {SYS_T::print_fatal("Error: get_nghostnode is not implemented. \n"); return 0;}

    virtual int get_ntotalnode() const
    {SYS_T::print_fatal("Error: get_ntotalnode is not implemented. \n"); return 0;}

    virtual int get_nbadnode() const
    {SYS_T::print_fatal("Error: get_nbadnode is not implemented. \n"); return 0;}

    virtual int get_nlocghonode() const
    {SYS_T::print_fatal("Error: get_nlocghonode is not implemented. \n"); return 0;}

    virtual bool isElemInPart(const int &gloindex) const
    {SYS_T::print_fatal("Error: isElemInPart is not implemented. \n"); return false;}

    virtual bool isNodeInPart(const int &gloindex) const
    {SYS_T::print_fatal("Error: isNodeInPart is not implemented. \n"); return false;}

    virtual int get_elemLocIndex(const int &gloindex) const
    {SYS_T::print_fatal("Error: get_elemLocIndex is not implemented. \n"); return 0;}

    virtual int get_nodeLocGhoIndex(const int &gloindex) const
    {SYS_T::print_fatal("Error: get_nodeLocGhoIndex is not implemented. \n"); return 0;}

    // 3. function that return the partition method parameters
    virtual int get_cpu_rank() const
    {SYS_T::print_fatal("Error: get_cpu_rank is not implemented. \n"); return 0;}

    virtual int get_cpu_size() const
    {SYS_T::print_fatal("Error: get_cpu_size is not implemented. \n"); return 0;}

    // 3. Global mesh information
    virtual int get_nElem() const
    {SYS_T::print_fatal("Error: get_nElem is not implemented. \n"); return 0;}

    virtual int get_nFunc() const
    {SYS_T::print_fatal("Error: get_nFunc is not implemented. \n"); return 0;}

    virtual int get_sDegree() const
    {SYS_T::print_fatal("Error: get_sDegree is not implemented. \n"); return 0;}

    virtual int get_tDegree() const
    {SYS_T::print_fatal("Error: get_tDegree is not implemented. \n"); return 0;}

    virtual int get_uDegree() const
    {SYS_T::print_fatal("Error: get_uDegree is not implemented. \n"); return 0;}

    virtual int get_nLocBas() const
    {SYS_T::print_fatal("Error: get_nLocBas is not implemented. \n"); return 0;}

    // 4. LIEN array
    virtual int get_LIEN(const int &ee, const int &ii) const
    {SYS_T::print_fatal("Error: get_LIEN is not implemented. \n"); return 0;}

    // 5. Control points in this processor
    virtual double get_ctrlPts_x_loc(const int &pos) const
    {SYS_T::print_fatal("Error: get_ctrlPts_x is not implemented. \n"); return 0.0;}

    virtual double get_ctrlPts_y_loc(const int &pos) const
    {SYS_T::print_fatal("Error: get_ctrlPts_y is not implemented. \n"); return 0.0;}

    virtual double get_ctrlPts_z_loc(const int &pos) const
    {SYS_T::print_fatal("Error: get_ctrlPts_z is not implemented. \n"); return 0.0;}

    virtual double get_ctrlPts_w_loc(const int &pos) const
    {SYS_T::print_fatal("Error: get_ctrlPts_w is not implemented. \n"); return 0.0;}

    // 6. print on screen the partition results 
    virtual void print_info() const
    {SYS_T::print_fatal("Error: print_info is not implemented. \n");}

    virtual void print_part_ele() const
    {SYS_T::print_fatal("Error: print_part_ele is not implemented. \n");}

    virtual void print_part_node() const
    {SYS_T::print_fatal("Error: print_part_node is not implemented. \n");}

    virtual void print_part_ghost_node() const
    {SYS_T::print_fatal("Error: print_part_ghost_node is not implemented. \n");}

    virtual void print_part_local_to_global() const
    {SYS_T::print_fatal("Error: print_part_local_to_global is not implemented. \n");}

    virtual void print_part_LIEN() const
    {SYS_T::print_fatal("Error: print_part_LIEN is not implemented. \n");}

    virtual void print_part_loadbalance_edgecut() const
    {SYS_T::print_fatal("Error: print_part_loadbalance_edgecut is not implemented. \n");}
};

#endif
