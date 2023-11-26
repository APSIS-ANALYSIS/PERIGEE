#ifndef PART_FEM_FSI_TISSUE_HPP
#define PART_FEM_FSI_TISSUE_HPP
// ============================================================================
// Part_FEM_FSI_Tissue.hpp
//
// Object: Partition 3D tet mesh into subdomains. 
//         The mesh is partitioned into sub-domains that are tagged with three 
//         local basis vectors denoting the radial, circumferential, and
//         longitudinal directions, as well as a mapping from global solid node
//         IEN to the index of the vector that recording the local basis vectors.
//         In the partitioned file, there will be the loc_basisr/c/l vectors 
//         with the length equaling to the number of local+ghost solid nodes, 
//         and the mapping with the length equaling to the number of local+ghost 
//         nodes.
//
// Date Created: Nov. 23 2023
// ============================================================================
#include "Part_FEM_FSI.hpp"

class Part_FEM_FSI_Tissue : public Part_FEM_FSI
{
  public:
    Part_FEM_FSI_Tissue( const IMesh * const &mesh,
        const IGlobal_Part * const &gpart,
        const Map_Node_Index * const &mnindex,
        const IIEN * const &IEN,
        const std::vector<double> &ctrlPts,
        const std::vector<int> &phytag,
        const std::vector<int> &node_f,
        const std::vector<int> &node_s,
        const std::vector<Vector_3> &basis_r,
        const std::vector<Vector_3> &basis_c,
        const std::vector<Vector_3> &basis_l,
        const int &in_cpu_rank, 
        const int &in_cpu_size,
        const int &in_elemType,
        const int &in_start_idx,
        const Field_Property &in_fp );

    virtual ~Part_FEM_FSI_Tissue() = default;

    virtual void write( const std::string &inputFileName ) const;

  protected:
    // local direction vecters. Letters r, c, and l denotes the radial, circumferential,
    // and longitudinal directions, respectively.
    std::vector<Vector_3> loc_basis_r, loc_basis_c, loc_basis_l;

    // A vector with length of nlocghonode, mapping from the total local+ghost nodes
    // to the solid local+ghost nodes. If node [ii] is a solid node, the mapping
    // records the corresponding index of the loc_basis std::vector, otherwise return -1.
    std::vector<int> node_locgho_solid;
    int nlocghonode_s;
};

#endif
