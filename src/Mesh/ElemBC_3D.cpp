#include "ElemBC_3D.hpp"

ElemBC_3D::ElemBC_3D( const FEType &elemtype ) 
: elem_type( elemtype ), num_ebc( 0 )
{
  num_node.clear();
  num_cell.clear();
  cell_nLocBas.clear();

  pt_xyz.clear();
  sur_ien.clear();
  global_node.clear();
  global_cell.clear();

  std::cout<<"===> ElemBC_3D called by an empty constructor is generated. \n";
}

ElemBC_3D::ElemBC_3D( const std::vector<std::string> &vtkfileList,
    const FEType &elemtype ) : elem_type( elemtype ), 
  num_ebc( static_cast<int>( vtkfileList.size() ) )
{
  num_node.resize(num_ebc);  
  num_cell.resize(num_ebc);
  cell_nLocBas.resize(num_ebc);

  pt_xyz.resize(num_ebc);
  sur_ien.resize(num_ebc);
  global_node.resize(num_ebc);
  global_cell.resize(num_ebc);

  std::cout<<"===> ElemBC_3D specified by \n";

  for(int ii=0; ii<num_ebc; ++ii)
  {
    std::cout<<"     ebc_id = "<<ii<<": "<<vtkfileList[ii]<<'\n';

    VTK_T::read_grid( vtkfileList[ii], num_node[ii], num_cell[ii], pt_xyz[ii], sur_ien[ii] );
    
    if(elem_type == FEType::Tet4)
      cell_nLocBas[ii] = 3; // linear triangle
    else if(elem_type == FEType::Tet10)
      cell_nLocBas[ii] = 6; // quadratic triangle
    else if(elem_type == FEType::Hex8)
      cell_nLocBas[ii] = 4; // bilinear quadrangle
    else if(elem_type == FEType::Hex27) 
      cell_nLocBas[ii] = 9; // biquadratic quadrangle
    else
      SYS_T::print_fatal("Error: ElemBC_3D constructor: unknown element type. \n");
    
    global_node[ii] = VTK_T::read_int_PointData( vtkfileList[ii], "GlobalNodeID");
    global_cell[ii] = VTK_T::read_int_CellData( vtkfileList[ii], "GlobalElementID");
  }

  std::cout<<"     is generated. \n";
}

void ElemBC_3D::print_info() const
{
  std::cout<<"========================= \n";
  std::cout<<"ElemBC_3D : ";
  //std::cout<<" elem_type = "<<elem_type<<'\t';
  if(elem_type == FEType::Tet4)
    std::cout<<" elem_type = Tet4"<<'\t';
  else if(elem_type == FEType::Tet10)
    std::cout<<" elem_type = Tet10"<<'\t';
  else if(elem_type == FEType::Hex8)
    std::cout<<" elem_type = Hex8"<<'\t';
  else if(elem_type == FEType::Hex27)
    std::cout<<" elem_type = Hex27"<<'\t';
  else
    std::cout<< " elem_type = Unknown"<<'\t';
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
    VEC_T::print(sur_ien[ii]);
  }
  std::cout<<"========================= \n";
}

void ElemBC_3D::resetSurIEN_outwardnormal( const IIEN * const &VIEN )
{
  if(elem_type == FEType::Tet4)
    reset501IEN_outwardnormal(VIEN); 
  else if(elem_type == FEType::Tet10)
    reset502IEN_outwardnormal(VIEN); 
  else if(elem_type == FEType::Hex8)
    reset601IEN_outwardnormal(VIEN);     
  else if(elem_type == FEType::Hex27)
    reset602IEN_outwardnormal(VIEN);     
  else SYS_T::print_fatal("Error: ElemBC_3D::resetSurIEN_outwardnormal function: unknown element type.\n");
}

void ElemBC_3D::reset501IEN_outwardnormal( const IIEN * const &VIEN )
{
  for(int ebcid=0; ebcid<num_ebc; ++ebcid)
  {
    TET_T::Tet4 * tetcell = new TET_T::Tet4();

    for(int ee=0; ee<num_cell[ebcid]; ++ee)
    {
      // Triangle mesh node index
      const int node_t[3] { get_ien(ebcid, ee, 0), get_ien(ebcid, ee, 1), get_ien(ebcid, ee, 2) };

      // The triangle mesh node's volumetric index
      const std::vector<int> node_t_gi { get_global_node(ebcid, node_t[0]), 
        get_global_node(ebcid, node_t[1]), get_global_node(ebcid, node_t[2]) };

      // cell ee's global/volumetric index  
      const int cell_gi = get_global_cell(ebcid, ee);

      // tet mesh first four node's volumetric index
      const int tet_n[4] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
          VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3) };

      // build the tet object
      tetcell->reset(tet_n[0], tet_n[1], tet_n[2], tet_n[3]);

      // determine the face id for this triangle in the tet object 
      const int tet_face_id = tetcell->get_face_id(node_t_gi[0], node_t_gi[1], node_t_gi[2]);

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
          SYS_T::print_fatal("Error: ElemBC_3D::reset501IEN_outwardnormal function: tet_face_id is out of range. \n");
          break;
      }
      ASSERT(pos0 >=0 && pos0 <=2, "Error: ElemBC_3D::reset501IEN_outwardnormal function logical error.\n" );
      ASSERT(pos1 >=0 && pos1 <=2, "Error: ElemBC_3D::reset501IEN_outwardnormal function logical error.\n" );
      ASSERT(pos2 >=0 && pos2 <=2, "Error: ElemBC_3D::reset501IEN_outwardnormal function logical error.\n" ); 

      // Now we have got the corrected ordering of node_t, put them back into
      // sur_ien.
      sur_ien[ebcid][3*ee+0] = node_t[pos0];
      sur_ien[ebcid][3*ee+1] = node_t[pos1];
      sur_ien[ebcid][3*ee+2] = node_t[pos2];
    }
    delete tetcell; 
  }
}

void ElemBC_3D::reset502IEN_outwardnormal( const IIEN * const &VIEN )
{
  for(int ebcid=0; ebcid<num_ebc; ++ebcid)
  {
    std::vector<int> node_t(6, 0);    // triange node index in 2D mesh
    std::vector<int> node_t_gi(6, 0); // triange node index in 3D mesh
    std::vector<int> tet_n(10, 0);    // tet node index in 3D mesh

    TET_T::Tet4 * tetcell = new TET_T::Tet4();

    for(int ee=0; ee<num_cell[ebcid]; ++ee)
    {
      for(int ii=0; ii<6; ++ii)
      {
        node_t[ii] = get_ien(ebcid, ee, ii);
        node_t_gi[ii] = get_global_node(ebcid, node_t[ii]);
      }

      const int cell_gi = get_global_cell(ebcid, ee);

      for(int ii=0; ii<10; ++ii) tet_n[ii] = VIEN->get_IEN(cell_gi, ii);

      tetcell->reset(tet_n[0], tet_n[1], tet_n[2], tet_n[3]);

      const int tet_face_id = tetcell->get_face_id(node_t_gi[0],
          node_t_gi[1], node_t_gi[2]);

      int pos0 = -1, pos1 = -1, pos2 = -1, pos3 = -1, pos4 = -1, pos5 = -1;

      switch( tet_face_id )
      {
        case 0:
          pos0 = VEC_T::get_pos(node_t_gi, tet_n[1]);
          pos1 = VEC_T::get_pos(node_t_gi, tet_n[2]);
          pos2 = VEC_T::get_pos(node_t_gi, tet_n[3]);
          pos3 = VEC_T::get_pos(node_t_gi, tet_n[5]);
          pos4 = VEC_T::get_pos(node_t_gi, tet_n[9]);
          pos5 = VEC_T::get_pos(node_t_gi, tet_n[8]);
          break;
        case 1:
          pos0 = VEC_T::get_pos(node_t_gi, tet_n[0]);
          pos1 = VEC_T::get_pos(node_t_gi, tet_n[3]);
          pos2 = VEC_T::get_pos(node_t_gi, tet_n[2]);
          pos3 = VEC_T::get_pos(node_t_gi, tet_n[7]);
          pos4 = VEC_T::get_pos(node_t_gi, tet_n[9]);
          pos5 = VEC_T::get_pos(node_t_gi, tet_n[6]);
          break;
        case 2:
          pos0 = VEC_T::get_pos(node_t_gi, tet_n[0]);
          pos1 = VEC_T::get_pos(node_t_gi, tet_n[1]);
          pos2 = VEC_T::get_pos(node_t_gi, tet_n[3]);
          pos3 = VEC_T::get_pos(node_t_gi, tet_n[4]);
          pos4 = VEC_T::get_pos(node_t_gi, tet_n[8]);
          pos5 = VEC_T::get_pos(node_t_gi, tet_n[7]);
          break;
        case 3:
          pos0 = VEC_T::get_pos(node_t_gi, tet_n[0]);
          pos1 = VEC_T::get_pos(node_t_gi, tet_n[2]);
          pos2 = VEC_T::get_pos(node_t_gi, tet_n[1]);
          pos3 = VEC_T::get_pos(node_t_gi, tet_n[6]);
          pos4 = VEC_T::get_pos(node_t_gi, tet_n[5]);
          pos5 = VEC_T::get_pos(node_t_gi, tet_n[4]);
          break;
        default:
          SYS_T::print_fatal("Error: ElemBC_3D::reset502IEN_outwardnormal function: tet_face_id is out of range. \n");
          break;
      }
      ASSERT(pos0 >=0 && pos0 <=5, "Error: ElemBC_3D::reset502IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos1 >=0 && pos1 <=5, "Error: ElemBC_3D::reset502IEN_outwardnormal function logical error.\n" );
      ASSERT(pos2 >=0 && pos2 <=5, "Error: ElemBC_3D::reset502IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos3 >=0 && pos3 <=5, "Error: ElemBC_3D::reset502IEN_outwardnormal function logical error.\n" );
      ASSERT(pos4 >=0 && pos4 <=5, "Error: ElemBC_3D::reset502IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos5 >=0 && pos5 <=5, "Error: ElemBC_3D::reset502IEN_outwardnormal function logical error.\n" );

      sur_ien[ebcid][6*ee+0] = node_t[pos0];
      sur_ien[ebcid][6*ee+1] = node_t[pos1];
      sur_ien[ebcid][6*ee+2] = node_t[pos2];
      sur_ien[ebcid][6*ee+3] = node_t[pos3];
      sur_ien[ebcid][6*ee+4] = node_t[pos4];
      sur_ien[ebcid][6*ee+5] = node_t[pos5];
    }
    delete tetcell;
  }
}

void ElemBC_3D::reset601IEN_outwardnormal( const IIEN * const &VIEN )
{
  for (int ebcid=0; ebcid<num_ebc; ++ebcid)
  {
    HEX_T::Hex8 * hexcell = new HEX_T::Hex8(); 

    for (int ee=0; ee<num_cell[ebcid]; ++ee)
    {
      // Quad mesh node index
      const int node_q[4] { get_ien(ebcid, ee, 0), get_ien(ebcid, ee, 1),  get_ien(ebcid, ee, 2), get_ien(ebcid, ee, 3) };  

      // The quad mesh node's volumetric index
      const std::vector<int> node_q_gi { get_global_node(ebcid, node_q[0]), get_global_node(ebcid, node_q[1]),  
                                         get_global_node(ebcid, node_q[2]), get_global_node(ebcid, node_q[3]) }; 

      // cell ee's global/volumetric index  
      const int cell_gi = get_global_cell(ebcid, ee);

      // hex mesh first eight node's volumetric index
      const int hex_n[8] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
        VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3),
        VIEN->get_IEN(cell_gi, 4), VIEN->get_IEN(cell_gi, 5),
        VIEN->get_IEN(cell_gi, 6), VIEN->get_IEN(cell_gi, 7) };

      // build the hex object
      hexcell->reset(hex_n[0], hex_n[1], hex_n[2], hex_n[3],
                     hex_n[4], hex_n[5], hex_n[6], hex_n[7]);

      // determind the face id for this quadrangle in the hex object
      const int hex_face_id = hexcell->get_face_id(node_q_gi[0], node_q_gi[1], node_q_gi[2], node_q_gi[3]);  

      int pos0 = -1, pos1 = -1, pos2 = -1, pos3 = -1;

      switch (hex_face_id)
      {
        case 0:
          pos0 = VEC_T::get_pos(node_q_gi, hex_n[0]);  
          pos1 = VEC_T::get_pos(node_q_gi, hex_n[3]); 
          pos2 = VEC_T::get_pos(node_q_gi, hex_n[2]); 
          pos3 = VEC_T::get_pos(node_q_gi, hex_n[1]); 
          break;
        case 1:
          pos0 = VEC_T::get_pos(node_q_gi, hex_n[4]);  
          pos1 = VEC_T::get_pos(node_q_gi, hex_n[5]); 
          pos2 = VEC_T::get_pos(node_q_gi, hex_n[6]); 
          pos3 = VEC_T::get_pos(node_q_gi, hex_n[7]); 
          break;
        case 2:
          pos0 = VEC_T::get_pos(node_q_gi, hex_n[0]);  
          pos1 = VEC_T::get_pos(node_q_gi, hex_n[1]); 
          pos2 = VEC_T::get_pos(node_q_gi, hex_n[5]); 
          pos3 = VEC_T::get_pos(node_q_gi, hex_n[4]); 
          break;
        case 3:
          pos0 = VEC_T::get_pos(node_q_gi, hex_n[1]); 
          pos1 = VEC_T::get_pos(node_q_gi, hex_n[2]); 
          pos2 = VEC_T::get_pos(node_q_gi, hex_n[6]); 
          pos3 = VEC_T::get_pos(node_q_gi, hex_n[5]); 
          break;
        case 4:
          pos0 = VEC_T::get_pos(node_q_gi, hex_n[2]);  
          pos1 = VEC_T::get_pos(node_q_gi, hex_n[3]); 
          pos2 = VEC_T::get_pos(node_q_gi, hex_n[7]); 
          pos3 = VEC_T::get_pos(node_q_gi, hex_n[6]); 
          break;
        case 5:
          pos0 = VEC_T::get_pos(node_q_gi, hex_n[0]);  
          pos1 = VEC_T::get_pos(node_q_gi, hex_n[4]); 
          pos2 = VEC_T::get_pos(node_q_gi, hex_n[7]); 
          pos3 = VEC_T::get_pos(node_q_gi, hex_n[3]); 
          break;          
        default:
          SYS_T::print_fatal("Error: ElemBC_3D::reset601IEN_outwardnormal function: hex_face_id is out of range. \n");
          break;
      }
      ASSERT(pos0 >=0 && pos0 <=3, "Error: ElemBC_3D::reset601IEN_outwardnormal function logical error.\n" );
      ASSERT(pos1 >=0 && pos1 <=3, "Error: ElemBC_3D::reset601IEN_outwardnormal function logical error.\n" );
      ASSERT(pos2 >=0 && pos2 <=3, "Error: ElemBC_3D::reset601IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos3 >=0 && pos3 <=3, "Error: ElemBC_3D::reset601IEN_outwardnormal function logical error.\n" ); 

      // Now we have got the corrected ordering of node_q, put them back into
      // sur_ien.
      sur_ien[ebcid][4*ee+0] = node_q[pos0];
      sur_ien[ebcid][4*ee+1] = node_q[pos1];
      sur_ien[ebcid][4*ee+2] = node_q[pos2];
      sur_ien[ebcid][4*ee+3] = node_q[pos3];
    }
    delete hexcell; 
  }
}

void ElemBC_3D::reset602IEN_outwardnormal( const IIEN * const &VIEN )
{
  for (int ebcid=0; ebcid<num_ebc; ++ebcid)
  {
    std::vector<int> node_q(9, 0);    // biquadratic quadrangle node index in 2D mesh
    std::vector<int> node_q_gi(9, 0); // biquadratic quadrangle node index in 3D mesh
    std::vector<int> hex_n(27, 0);    // triquadratic hex node index in 3D mesh

    HEX_T::Hex8 * hexcell = new HEX_T::Hex8();  

    for (int ee=0; ee<num_cell[ebcid]; ++ee)
    {
      for (int ii=0; ii<9; ++ii)
      {
        node_q[ii]    = get_ien(ebcid, ee, ii);
        node_q_gi[ii] = get_global_node(ebcid, node_q[ii]);
      }

      const int cell_gi = get_global_cell(ebcid, ee);

      for(int ii=0; ii<27; ++ii) hex_n[ii] = VIEN->get_IEN(cell_gi, ii);

      // build the hex object
      hexcell->reset(hex_n[0], hex_n[1], hex_n[2], hex_n[3],
                     hex_n[4], hex_n[5], hex_n[6], hex_n[7]);

      // determind the face id for this quadrangle in the hex object
      const int hex_face_id = hexcell->get_face_id(node_q_gi[0], node_q_gi[1], node_q_gi[2], node_q_gi[3]); 

      int pos0 = -1, pos1 = -1, pos2 = -1, pos3 = -1, pos4 = -1, pos5 = -1, pos6 = -1, pos7 = -1, pos8 = -1;

      switch (hex_face_id)
      {
        case 0:
          pos0 = VEC_T::get_pos(node_q_gi, hex_n[0]); 
          pos1 = VEC_T::get_pos(node_q_gi, hex_n[3]); 
          pos2 = VEC_T::get_pos(node_q_gi, hex_n[2]); 
          pos3 = VEC_T::get_pos(node_q_gi, hex_n[1]); 
          pos4 = VEC_T::get_pos(node_q_gi, hex_n[11]); 
          pos5 = VEC_T::get_pos(node_q_gi, hex_n[10]); 
          pos6 = VEC_T::get_pos(node_q_gi, hex_n[9]); 
          pos7 = VEC_T::get_pos(node_q_gi, hex_n[8]); 
          pos8 = VEC_T::get_pos(node_q_gi, hex_n[24]); 
          break;         
        case 1:
          pos0 = VEC_T::get_pos(node_q_gi, hex_n[4]);  
          pos1 = VEC_T::get_pos(node_q_gi, hex_n[5]); 
          pos2 = VEC_T::get_pos(node_q_gi, hex_n[6]); 
          pos3 = VEC_T::get_pos(node_q_gi, hex_n[7]);
          pos4 = VEC_T::get_pos(node_q_gi, hex_n[12]); 
          pos5 = VEC_T::get_pos(node_q_gi, hex_n[13]); 
          pos6 = VEC_T::get_pos(node_q_gi, hex_n[14]); 
          pos7 = VEC_T::get_pos(node_q_gi, hex_n[15]); 
          pos8 = VEC_T::get_pos(node_q_gi, hex_n[25]); 
          break;
        case 2:
          pos0 = VEC_T::get_pos(node_q_gi, hex_n[0]);  
          pos1 = VEC_T::get_pos(node_q_gi, hex_n[1]); 
          pos2 = VEC_T::get_pos(node_q_gi, hex_n[5]); 
          pos3 = VEC_T::get_pos(node_q_gi, hex_n[4]); 
          pos4 = VEC_T::get_pos(node_q_gi, hex_n[8]); 
          pos5 = VEC_T::get_pos(node_q_gi, hex_n[17]); 
          pos6 = VEC_T::get_pos(node_q_gi, hex_n[12]); 
          pos7 = VEC_T::get_pos(node_q_gi, hex_n[16]); 
          pos8 = VEC_T::get_pos(node_q_gi, hex_n[22]); 
          break;
        case 3:
          pos0 = VEC_T::get_pos(node_q_gi, hex_n[1]);  
          pos1 = VEC_T::get_pos(node_q_gi, hex_n[2]); 
          pos2 = VEC_T::get_pos(node_q_gi, hex_n[6]); 
          pos3 = VEC_T::get_pos(node_q_gi, hex_n[5]);
          pos4 = VEC_T::get_pos(node_q_gi, hex_n[9]); 
          pos5 = VEC_T::get_pos(node_q_gi, hex_n[18]); 
          pos6 = VEC_T::get_pos(node_q_gi, hex_n[13]); 
          pos7 = VEC_T::get_pos(node_q_gi, hex_n[17]); 
          pos8 = VEC_T::get_pos(node_q_gi, hex_n[21]);  
          break;
         case 4:
          pos0 = VEC_T::get_pos(node_q_gi, hex_n[2]); 
          pos1 = VEC_T::get_pos(node_q_gi, hex_n[3]); 
          pos2 = VEC_T::get_pos(node_q_gi, hex_n[7]); 
          pos3 = VEC_T::get_pos(node_q_gi, hex_n[6]); 
          pos4 = VEC_T::get_pos(node_q_gi, hex_n[10]); 
          pos5 = VEC_T::get_pos(node_q_gi, hex_n[19]); 
          pos6 = VEC_T::get_pos(node_q_gi, hex_n[14]); 
          pos7 = VEC_T::get_pos(node_q_gi, hex_n[18]); 
          pos8 = VEC_T::get_pos(node_q_gi, hex_n[23]); 
          break;
         case 5:
          pos0 = VEC_T::get_pos(node_q_gi, hex_n[0]);  
          pos1 = VEC_T::get_pos(node_q_gi, hex_n[4]); 
          pos2 = VEC_T::get_pos(node_q_gi, hex_n[7]); 
          pos3 = VEC_T::get_pos(node_q_gi, hex_n[3]); 
          pos4 = VEC_T::get_pos(node_q_gi, hex_n[16]); 
          pos5 = VEC_T::get_pos(node_q_gi, hex_n[15]); 
          pos6 = VEC_T::get_pos(node_q_gi, hex_n[19]); 
          pos7 = VEC_T::get_pos(node_q_gi, hex_n[11]); 
          pos8 = VEC_T::get_pos(node_q_gi, hex_n[20]); 
          break;         
        default:
          SYS_T::print_fatal("Error: ElemBC_3D::reset602IEN_outwardnormal function: hex_face_id is out of range. \n");
          break;
      }
      ASSERT(pos0 >=0 && pos0 <=8, "Error: ElemBC_3D::reset602IEN_outwardnormal function logical error.\n" );
      ASSERT(pos1 >=0 && pos1 <=8, "Error: ElemBC_3D::reset602IEN_outwardnormal function logical error.\n" );
      ASSERT(pos2 >=0 && pos2 <=8, "Error: ElemBC_3D::reset602IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos3 >=0 && pos3 <=8, "Error: ElemBC_3D::reset602IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos4 >=0 && pos4 <=8, "Error: ElemBC_3D::reset602IEN_outwardnormal function logical error.\n" );
      ASSERT(pos5 >=0 && pos5 <=8, "Error: ElemBC_3D::reset602IEN_outwardnormal function logical error.\n" );
      ASSERT(pos6 >=0 && pos6 <=8, "Error: ElemBC_3D::reset602IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos7 >=0 && pos7 <=8, "Error: ElemBC_3D::reset602IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos8 >=0 && pos8 <=8, "Error: ElemBC_3D::reset602IEN_outwardnormal function logical error.\n" ); 

      sur_ien[ebcid][9*ee+0] = node_q[pos0];
      sur_ien[ebcid][9*ee+1] = node_q[pos1];
      sur_ien[ebcid][9*ee+2] = node_q[pos2];
      sur_ien[ebcid][9*ee+3] = node_q[pos3];
      sur_ien[ebcid][9*ee+4] = node_q[pos4];
      sur_ien[ebcid][9*ee+5] = node_q[pos5];
      sur_ien[ebcid][9*ee+6] = node_q[pos6];
      sur_ien[ebcid][9*ee+7] = node_q[pos7];
      sur_ien[ebcid][9*ee+8] = node_q[pos8];
    }
    delete hexcell; 
  }
}

// EOF
