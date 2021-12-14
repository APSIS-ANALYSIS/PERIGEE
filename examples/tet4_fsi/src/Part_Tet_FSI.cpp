#include "Part_Tet_FSI.hpp"

Part_Tet_FSI::Part_Tet_FSI( const IMesh * const &mesh_p,
    const IMesh * const &mesh_v,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex,
    const Map_Node_Index * const &mnindex_p,
    const Map_Node_Index * const &mnindex_v,
    const IIEN * const &IEN_p,
    const IIEN * const &IEN_v,
    const std::vector<double> &ctrlPts,
    const int &in_start_idx_p, const int &in_start_idx_v,
    const int &in_cpu_rank, const int &in_cpu_size,
    const int &in_dofNum, const int &in_dofMat,
    const int &in_elemType,
    const bool isPrintInfo )
: cpu_rank(in_cpu_rank), cpu_size(in_cpu_size),
  dual_edge_ncommon( gpart->get_dual_edge_ncommon() ),
{}


Part_Tet_FSI::~Part_Tet_FSI()
{}



// EOF
