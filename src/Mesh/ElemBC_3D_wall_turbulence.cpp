#include "ElemBC_3D_wall_turbulence.hpp"

ElemBC_3D_wall_turbulence::ElemBC_3D_wall_turbulence( const std::vector<std::string> &vtkfileList,
    const int &in_weak_bc_type, const double &in_C_bI, const IIEN * const &VIEN, const int &elemtype )
: ElemBC_3D ( elemtype ), weak_bc_type {in_weak_bc_type}, C_bI {in_C_bI}
{
  SYS_T::print_fatal_if(VEC_T::get_size(vtkfileList) > 1,
    "Error, ElemBC_3D_wall_turbulence: The number of wall file should not be more than 1.\n");

  if(VEC_T::get_size(vtkfileList) == 1)
  {
    face_id.resize(num_cell[0]);

    if(elem_type == 501 || elem_type == 502)
    {
      TET_T::Tet4 * tetcell = new TET_T::Tet4();

      for(int ee{0}; ee < num_cell[0]; ++ee)
      {
        const int node_t[3] { get_ien(0, ee, 0), get_ien(0, ee, 1), get_ien(0, ee, 2) };

        const int node_t_gi[3] { get_global_node(0, node_t[0]),
                                 get_global_node(0, node_t[1]),
                                 get_global_node(0, node_t[2]) };

        const int cell_gi = get_global_cell(0, ee);

        const int tet_n[4] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                             VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3) };
        
        tetcell->reset(tet_n[0], tet_n[1], tet_n[2], tet_n[3]);

        face_id[ee] = tetcell->get_face_id(node_t_gi[0], node_t_gi[1], node_t_gi[2]);
      }
    }
    else if(elem_type == 601 || elem_type == 602)
    {
      HEX_T::Hex8 * hexcell = new HEX_T::Hex8();

      for(int ee{0}; ee < num_cell[0]; ++ee)
      {
        const int node_q[4] { get_ien(0, ee, 0), get_ien(0, ee, 1),
                              get_ien(0, ee, 2), get_ien(0, ee, 3) };
        
        const int node_q_gi[4] { get_global_node(0, node_q[0]),
                                 get_global_node(0, node_q[1]),
                                 get_global_node(0, node_q[2]),
                                 get_global_node(0, node_q[3]) };
        
        const int cell_gi = get_global_cell(0, ee);

        const int hex_n[8] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                             VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3),
                             VIEN->get_IEN(cell_gi, 4), VIEN->get_IEN(cell_gi, 5),
                             VIEN->get_IEN(cell_gi, 6), VIEN->get_IEN(cell_gi, 7) };
        
        hexcell->reset(hex_n[0], hex_n[1], hex_n[2], hex_n[3],
                       hex_n[4], hex_n[5], hex_n[6], hex_n[7]);

        face_id[ee] = hexcell->get_face_id(node_q_gi[0], node_q_gi[1], node_q_gi[2], node_q_gi[3]);
      }
    }
    else
      SYS_T::print_fatal("Error: ElemBC_3D_wall_turbulence, unknown element type.\n");
  }
  else
    ; // do nothing
}

ElemBC_3D_wall_turbulence::~ElemBC_3D_wall_turbulence()
{
  face_id.clear();
}