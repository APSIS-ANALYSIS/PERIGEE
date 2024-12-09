#include "ElemBC_3D_turbulence_wall_model.hpp"

ElemBC_3D_turbulence_wall_model::ElemBC_3D_turbulence_wall_model( 
    const std::vector<std::string> &vtkfileList,
    const int &in_wall_model_type, const IIEN * const &VIEN, 
    const FEType &elemtype )
: ElemBC_3D ( vtkfileList, elemtype ), 
  wall_model_type {in_wall_model_type}
{
  SYS_T::print_fatal_if(VEC_T::get_size(vtkfileList) > 1, "Error, ElemBC_3D_turbulence_wall_model: The number of wall file should not be more than 1.\n");

  if(VEC_T::get_size(vtkfileList) == 1)
  {
    face_id.resize(num_cell[0]);

    if(elem_type == FEType::Tet4 || elem_type == FEType::Tet10)
    {
      TET_T::Tet4 * tetcell = new TET_T::Tet4();

      for(int ee{0}; ee < num_cell[0]; ++ee)
      {
        const int node_t[3] { get_ien(0, ee, 0), get_ien(0, ee, 1), get_ien(0, ee, 2) };

        const std::array<int,3> node_t_gi {{ get_global_node(0, node_t[0]),
                                             get_global_node(0, node_t[1]),
                                             get_global_node(0, node_t[2]) }};

        const int cell_gi = get_global_cell(0, ee);

        tetcell->reset( VIEN->get_IEN_array4( cell_gi ) );

        face_id[ee] = tetcell->get_face_id( node_t_gi );
      }

      delete tetcell;
    }
    else if(elem_type == FEType::Hex8 || elem_type == FEType::Hex27)
    {
      HEX_T::Hex8 * hexcell = new HEX_T::Hex8();

      for(int ee{0}; ee < num_cell[0]; ++ee)
      {
        const int node_q[4] { get_ien(0, ee, 0), get_ien(0, ee, 1),
                              get_ien(0, ee, 2), get_ien(0, ee, 3) };
        
        const std::array<int,4> node_q_gi {{ get_global_node(0, node_q[0]),
                                             get_global_node(0, node_q[1]),
                                             get_global_node(0, node_q[2]),
                                             get_global_node(0, node_q[3]) }};
        
        const int cell_gi = get_global_cell(0, ee);

        hexcell->reset( VIEN->get_IEN_array8( cell_gi ) );

        face_id[ee] = hexcell->get_face_id( node_q_gi );
      }

      delete hexcell;
    }
    else
      SYS_T::print_fatal("Error: ElemBC_3D_turbulence_wall_model, unknown element type.\n");
  }
  else
    SYS_T::commPrint("Warning: there is no geometry file provided for the ElemBC_3D_turbulence_wall_model class. \n" );

}

// EOF
