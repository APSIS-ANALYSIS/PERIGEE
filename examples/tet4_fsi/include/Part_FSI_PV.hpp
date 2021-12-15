#ifndef PART_FSI_PV_HPP
#define PART_FSI_PV_HPP
// ============================================================================
// Part_FSI_PV.hpp
//
// Date: Dec. 14 2021
// ============================================================================
#include "IMesh.hpp"
#include "IGlobal_Part.hpp"
#include "Map_Node_Index.hpp"
#include "IIEN.hpp"

class Part_FSI_PV
{
  public:
    Part_FSI_PV( const IMesh * const &mesh_p,
    const IMesh * const &mesh_v,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex_p,
    const Map_Node_Index * const &mnindex_v,
    const IIEN * const &IEN_p,
    const IIEN * const &IEN_v,
    const std::vector<double> &ctrlPts,
    const std::vector<int> &phytag,
    const std::vector<int> &p_node_f,
    const std::vector<int> &p_node_s,
    const std::vector<int> &v_node_f,
    const std::vector<int> &v_node_s,
    const int &in_start_idx_p, const int &in_start_idx_v,
    const int &in_cpu_rank, const int &in_cpu_size,
    const int &in_elemType,
    const bool isPrintInfo );

    virtual ~Part_FSI_PV();

  protected:
    // 1. local element (===> ALocal_Elem)
    std::vector<int> elem_loc, elem_phy_tag;
    int nlocalele;

    // 2. local node (===> APart_Node)
    std::vector<int> node_loc_p, node_loc_v;
    std::vector<int> node_loc_original_p, node_loc_original_v;
    std::vector<int> node_ghost_p, node_ghost_v;
    std::vector<int> local_to_global_p, local_to_global_v;

    int nlocalnode_p, nlocalnode_v;
    int nghostnode_p, nghostnode_v;
    int ntotalnode_p, ntotalnode_v;
    int nlocghonode_p, nlocghonode_v;

    // The location in the node_loc for the local fluid/solid nodes
    // e.g. node_loc_fluid.size() gives the number of nodes belonging
    // to a fluid domain in this partition, and node_loc_fluid[ii] gives
    // the location in node_loc for that node's index.
    // Storing the index in node_loc array makes things more flexible,
    // because sometimes, we manage arrays using local partition's numbering
    // sometimes, we directly set entries using the global numbering.
    std::vector<int> node_loc_p_fluid, node_loc_p_solid;
    std::vector<int> node_loc_v_fluid, node_loc_v_solid;

    int nlocalnode_p_fluid, nlocalnode_p_solid;
    int nlocalnode_v_fluid, nlocalnode_v_solid;

    // 3 CPU info and partition parameters (===> APart_Basic_Info)
    const int cpu_rank, cpu_size, dual_edge_ncommon;

    // 4. global mesh info (===> AGlobal_Mesh_Info)
    const int nElem, probDim, elemType;
    const int nFunc_p, sDegree_p, tDegree_p, uDegree_p, nLocBas_p;
    const int nFunc_v, sDegree_v, tDegree_v, uDegree_v, nLocBas_v;

    // 5. LIEN (===> ALocal_IEN)
    std::vector<int> pLIEN, vLIEN; // LIEN for separate fields

    // 6. control points
    std::vector<double> ctrlPts_x_loc, ctrlPts_y_loc, ctrlPts_z_loc;
};

#endif
