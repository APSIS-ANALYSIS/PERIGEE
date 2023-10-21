#include "ElemBC_3D_weak.hpp"

ElemBC_3D_weak::ElemBC_3D_weak( const std::vector<std::string> &vtkfileList,
    const int &in_weak_bc_type, const double &in_C_bI, const IIEN * const &VIEN, const int &elemtype )
: ElemBC_3D ( elemtype ), weak_bc_type {in_weak_bc_type}, C_bI {in_C_bI}
{
  for(int ebcid{0}; ebcid<num_ebc; ++ebcid)
  {
    face_id[ebcid].resize(num_cell[ebcid]);

    if(elem_type == 501 || elem_type == 502)
    {
      TET_T::Tet4 * tetcell = new TET_T::Tet4();

      for(int ee{0}; ee < num_cell[ebcid]; ++ee)
      {
        const int node_t[3] { get_ien(ebcid, ee, 0), get_ien(ebcid, ee, 1), get_ien(ebcid, ee, 2) };

        const int node_t_gi[3] { get_global_node(ebcid, node_t[0]),
                                 get_global_node(ebcid, node_t[1]),
                                 get_global_node(ebcid, node_t[2]) };

        const int cell_gi = get_global_cell(ebcid, ee);

        const int tet_n[4] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                             VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3) };
        
        tetcell->reset(tet_n[0], tet_n[1], tet_n[2], tet_n[3]);

        face_id[ebcid][ee] = tetcell->get_face_id(node_t_gi[0], node_t_gi[1], node_t_gi[2]);
      }
    }
    else if(elem_type == 601 || elem_type == 602)
    {
      HEX_T::Hex8 * hexcell = new HEX_T::Hex8();

      for(int ee{0}; ee < num_cell[ebcid]; ++ee)
      {
        const int node_q[4] { get_ien(ebcid, ee, 0), get_ien(ebcid, ee, 1),
                              get_ien(ebcid, ee, 2), get_ien(ebcid, ee, 3) };
        
        const int node_q_gi[4] { get_global_node(ebcid, node_q[0]),
                                 get_global_node(ebcid, node_q[1]),
                                 get_global_node(ebcid, node_q[2]),
                                 get_global_node(ebcid, node_q[3]) };
        
        const int cell_gi = get_global_cell(ebcid, ee);

        const int hex_n[8] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                             VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3),
                             VIEN->get_IEN(cell_gi, 4), VIEN->get_IEN(cell_gi, 5),
                             VIEN->get_IEN(cell_gi, 6), VIEN->get_IEN(cell_gi, 7) };
        
        hexcell->reset(hex_n[0], hex_n[1], hex_n[2], hex_n[3],
                       hex_n[4], hex_n[5], hex_n[6], hex_n[7]);

        face_id[ebcid][ee] = hexcell->get_face_id(node_q_gi[0], node_q_gi[1], node_q_gi[2], node_q_gi[3]);
      }
    }
    else
      SYS_T::print_fatal("Error: ElemBC_3D_weak, unknown element type.\n");
  }
}

ElemBC_3D_weak::~ElemBC_3D_weak()
{
  face_id.clear();
}

void ElemBC_3D_weak::set_Q(const std::vector<int> &wall_node_idx, const std::vector<Tensor2_3D> &Q_at_node)
{
  if(weak_bc_type != 2)
  {
    SYS_T::commPrint("No need to set rotation matrices at nodes.\n");
    return;
  }
  else
  {
    Q.resize(num_ebc);
    for(int ebcid {0}; ebcid < num_ebc; ++ebcid)
    {
      Q[ebcid].resize(9 * num_node[ebcid]);
      for(int node{0}; node < num_node[ebcid]; ++node)
      {
        // The location of node in the wall_node_idx vector
        const int location = VEC_T::get_pos(wall_node_idx, global_node[ebcid][node]);

        for(int comp{0}; comp < 9; ++comp)
            Q[ebcid][9 * node + comp] = Q_at_node[location](comp);
      }
    }
  }
}