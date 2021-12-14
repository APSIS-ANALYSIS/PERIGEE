#ifndef PART_TET_FSI_HPP
#define PART_TET_FSI_HPP
// ============================================================================
// Part_Tet_FSI.hpp
//
// Date: Dec. 14 2021
// ============================================================================
#include "IMesh.hpp"
#include "IGlobal_Part.hpp"
#include "Map_Node_Index.hpp"
#include "IIEN.hpp"

class Part_Tet_FSI
{
  public:
    Part_Tet_FSI( const IMesh * const &mesh,
        const IGlobal_Part * const &gpart,
        const Map_Node_Index * const &mnindex,
        const IIEN * const &IEN,
        const std::vector<double> &ctrlPts,
        const int &in_cpu_rank, const int &in_cpu_size,
        const int &in_dofNum, const int &in_dofMat,
        const int &in_elemType,
        const bool isPrintInfo );

    virtual ~Part_Tet_FSI();

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
    int nbadnode_p, nbadnode_v;
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
    int cpu_rank, cpu_size, dual_edge_ncommon;
    
    // 4. global mesh info (===> AGlobal_Mesh_Info)
    int nElem, nFunc, sDegree, tDegree, uDegree, nLocBas, probDim, dofNum, dofMat, elemType;

    // 5. LIEN (===> ALocal_IEN)
    std::vector<int> pLIEN, vLIEN; // LIEN for separate fields
    
    // 6. control points
    std::vector<double> ctrlPts_x_loc, ctrlPts_y_loc, ctrlPts_z_loc;
};

#endif
