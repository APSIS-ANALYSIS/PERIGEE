#include "ElemBC_3D_tet4.hpp"

ElemBC_3D_tet4::ElemBC_3D_tet4( const std::vector<std::string> &vtpfileList )
{
  num_ebc = static_cast<int>( vtpfileList.size() );

  num_node     = new int [num_ebc];
  num_cell     = new int [num_ebc];
  cell_nLocBas = new int [num_ebc];

  pt_xyz.resize(num_ebc);
  tri_ien.resize(num_ebc);
  global_node.resize(num_ebc);
  global_cell.resize(num_ebc);

  std::cout<<"===> ElemBC_3D_tet4 specified by \n";

  for(int ii=0; ii<num_ebc; ++ii)
  {
    std::cout<<"     ebc_id = "<<ii<<": "<<vtpfileList[ii]<<'\n';

    cell_nLocBas[ii] = 3; // linear tests corresponds to linear triangle
    TET_T::read_vtp_grid( vtpfileList[ii], num_node[ii], num_cell[ii], 
        pt_xyz[ii], tri_ien[ii], global_node[ii], global_cell[ii] );
  }

  std::cout<<"     is generated. \n";
}


ElemBC_3D_tet4::~ElemBC_3D_tet4()
{
  delete [] num_node; num_node = NULL;
  delete [] num_cell; num_cell = NULL;
  delete [] cell_nLocBas; cell_nLocBas = NULL;
}


void ElemBC_3D_tet4::print_info() const
{
  std::cout<<"========================= \n";
  std::cout<<"ElemBC_3D_tet4 : ";
  std::cout<<" num_ebc = "<<num_ebc<<std::endl;
  for(int ii=0; ii<num_ebc; ++ii)
  {
    std::cout<<"== ebc_id = "<<ii<<'\n';
    std::cout<<" num_node = "<<num_node[ii]<<'\t';
    std::cout<<" num_cell = "<<num_cell[ii]<<'\t';
    std::cout<<" cell_nLocBas = "<<cell_nLocBas[ii]<<'\n';
    std::cout<<" global node : ";
    VEC_T::print(global_node[ii]);
    std::cout<<" global elem : ";
    VEC_T::print(global_cell[ii]);
    std::cout<<" surface IEN : ";
    VEC_T::print(tri_ien[ii]);
  }
  std::cout<<"========================= \n";
}


void ElemBC_3D_tet4::resetTriIEN_outwardnormal( const IIEN * const &VIEN )
{
  for(int etype = 0; etype < num_ebc; ++etype)
  {
    int e_num_cell = num_cell[etype];
    std::vector<int> node_t(3, 0);
    // the cell's first three nodes' volumetric indices
    std::vector<int> node_t_gi(3, 0);
    
    int cell_gi; // this cell's global index

    std::vector<int> tet_n(4,0);

    TET_T::Tet4 * tetcell = new TET_T::Tet4();

    for(int ee=0; ee<e_num_cell; ++ee)
    {
      node_t[0] = get_ien(etype, ee, 0);
      node_t[1] = get_ien(etype, ee, 1);
      node_t[2] = get_ien(etype, ee, 2);

      node_t_gi[0] = get_global_node(etype, node_t[0]);
      node_t_gi[1] = get_global_node(etype, node_t[1]);
      node_t_gi[2] = get_global_node(etype, node_t[2]);
      
      cell_gi = get_global_cell(etype, ee);

      tet_n[0] = VIEN->get_IEN(cell_gi, 0);
      tet_n[1] = VIEN->get_IEN(cell_gi, 1);
      tet_n[2] = VIEN->get_IEN(cell_gi, 2);
      tet_n[3] = VIEN->get_IEN(cell_gi, 3);

      tetcell->reset(tet_n[0], tet_n[1], tet_n[2], tet_n[3]);
      
      int tet_face_id = tetcell->get_face_id(node_t_gi[0], 
          node_t_gi[1], node_t_gi[2]);

      int pos0 = -1, pos1 = -1, pos2 = -1;
      switch( tet_face_id )
      {
        case 0:
          pos0 = VEC_T::get_pos(node_t_gi, tet_n[1]);
          pos1 = VEC_T::get_pos(node_t_gi, tet_n[2]);
          pos2 = VEC_T::get_pos(node_t_gi, tet_n[3]);
          break;
        case 1:
          pos0 = VEC_T::get_pos(node_t_gi, tet_n[0]);
          pos1 = VEC_T::get_pos(node_t_gi, tet_n[3]);
          pos2 = VEC_T::get_pos(node_t_gi, tet_n[2]);
          break;
        case 2:
          pos0 = VEC_T::get_pos(node_t_gi, tet_n[0]);
          pos1 = VEC_T::get_pos(node_t_gi, tet_n[1]);
          pos2 = VEC_T::get_pos(node_t_gi, tet_n[3]);
          break;
        case 3:
          pos0 = VEC_T::get_pos(node_t_gi, tet_n[0]);
          pos1 = VEC_T::get_pos(node_t_gi, tet_n[2]);
          pos2 = VEC_T::get_pos(node_t_gi, tet_n[1]);
          break;
        default:
          SYS_T::print_fatal("Error: resetTriIEN_outwardnormal : tet_face_id is out of range. \n");
          break;
      }
      assert(pos0 >=0 && pos0 <=2);
      assert(pos1 >=0 && pos1 <=2);
      assert(pos2 >=0 && pos2 <=2); 

      // Now we have got the corrected ordering of node_t, put them back into
      // tri_ien.
      tri_ien[etype][3*ee+0] = node_t[pos0];
      tri_ien[etype][3*ee+1] = node_t[pos1];
      tri_ien[etype][3*ee+2] = node_t[pos2];
    }
    delete tetcell; 
  }
}

// EOF
