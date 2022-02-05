#include "ElemBC_3D_tet.hpp"

ElemBC_3D_tet::ElemBC_3D_tet( const int &elemtype ) 
: elem_type( elemtype ), num_ebc( 0 )
{
  num_node     = nullptr;
  num_cell     = nullptr;
  cell_nLocBas = nullptr;

  pt_xyz.clear();
  tri_ien.clear();
  global_node.clear();
  global_cell.clear();

  std::cout<<"===> ElemBC_3D_tet called by an empty constructor is generated. \n";
}


ElemBC_3D_tet::ElemBC_3D_tet( const std::string &vtkfile,
   const int &elemtype ) : elem_type( elemtype ), num_ebc( 1 )
{
  num_node     = new int [num_ebc];
  num_cell     = new int [num_ebc];
  cell_nLocBas = new int [num_ebc];

  pt_xyz.resize(num_ebc);
  tri_ien.resize(num_ebc);
  global_node.resize(num_ebc);
  global_cell.resize(num_ebc);

  std::cout<<"===> ElemBC_3D_tet specified by "<<vtkfile<<'\n';

  if(elemtype == 501)
  {
    cell_nLocBas[0] = 3; // linear triangle
    TET_T::read_vtp_grid( vtkfile, num_node[0], num_cell[0],
        pt_xyz[0], tri_ien[0], global_node[0], global_cell[0] );
  }
  else if(elemtype == 502)
  {
    cell_nLocBas[0] = 6; // quadratic triangle
    TET_T::read_vtu_grid( vtkfile, num_node[0], num_cell[0],
        pt_xyz[0], tri_ien[0], global_node[0], global_cell[0] );
  }

  std::cout<<"     is generated. \n";
}


ElemBC_3D_tet::ElemBC_3D_tet( const std::vector<std::string> &vtkfileList,
    const int &elemtype ) : elem_type( elemtype ), 
  num_ebc( static_cast<int>( vtkfileList.size() ) )
{
  num_node     = new int [num_ebc];
  num_cell     = new int [num_ebc];
  cell_nLocBas = new int [num_ebc];

  pt_xyz.resize(num_ebc);
  tri_ien.resize(num_ebc);
  global_node.resize(num_ebc);
  global_cell.resize(num_ebc);

  std::cout<<"===> ElemBC_3D_tet specified by \n";

  for(int ii=0; ii<num_ebc; ++ii)
  {
    std::cout<<"     ebc_id = "<<ii<<": "<<vtkfileList[ii]<<'\n';

    if(elemtype == 501)
    {
      cell_nLocBas[ii] = 3; // linear triangle
      TET_T::read_vtp_grid( vtkfileList[ii], num_node[ii], num_cell[ii],
          pt_xyz[ii], tri_ien[ii], global_node[ii], global_cell[ii] );
    }
    else if(elemtype == 502)
    {
      cell_nLocBas[ii] = 6; // quadratic triangle
      TET_T::read_vtu_grid( vtkfileList[ii], num_node[ii], num_cell[ii],
          pt_xyz[ii], tri_ien[ii], global_node[ii], global_cell[ii] );
    }
  }

  std::cout<<"     is generated. \n";
}


ElemBC_3D_tet::~ElemBC_3D_tet()
{
  delete [] num_node; num_node = nullptr;
  delete [] num_cell; num_cell = nullptr;
  delete [] cell_nLocBas; cell_nLocBas = nullptr;
}


void ElemBC_3D_tet::print_info() const
{
  std::cout<<"========================= \n";
  std::cout<<"ElemBC_3D_tet : ";
  std::cout<<" elem_type = "<<elem_type<<'\t';
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


void ElemBC_3D_tet::resetTriIEN_outwardnormal( const IIEN * const &VIEN )
{
  if(elem_type == 501)
  {
    for(int ebcid = 0; ebcid < num_ebc; ++ebcid)
    {
      TET_T::Tet4 * tetcell = new TET_T::Tet4();

      for(int ee=0; ee<num_cell[ebcid]; ++ee)
      {
        // Triangle mesh node index
        const int node_t[3] { get_ien(ebcid, ee, 0), get_ien(ebcid, ee, 1), get_ien(ebcid, ee, 2) };

        // The triangle mesh node's volumetric index
        const int node_t_gi[3] { get_global_node(ebcid, node_t[0]), 
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
            SYS_T::print_fatal("Error: resetTriIEN_outwardnormal : tet_face_id is out of range. \n");
            break;
        }
        assert(pos0 >=0 && pos0 <=2);
        assert(pos1 >=0 && pos1 <=2);
        assert(pos2 >=0 && pos2 <=2); 

        // Now we have got the corrected ordering of node_t, put them back into
        // tri_ien.
        tri_ien[ebcid][3*ee+0] = node_t[pos0];
        tri_ien[ebcid][3*ee+1] = node_t[pos1];
        tri_ien[ebcid][3*ee+2] = node_t[pos2];
      }
      delete tetcell; 
    }
  }
  else if(elem_type == 502)
  {
    for(int ebcid = 0; ebcid < num_ebc; ++ebcid)
    {
      std::vector<int> node_t(6, 0); // triange node index in 2D mesh
      std::vector<int> node_t_gi(6, 0); // triange node index in 3D mesh
      std::vector<int> tet_n(10,0); // tet node index in 3D mesh

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
            SYS_T::print_fatal("Error: resetTriIEN_outwardnormal : tet_face_id is out of range. \n");
            break;
        }
        assert(pos0 >=0 && pos0 <=5); assert(pos1 >=0 && pos1 <=5);
        assert(pos2 >=0 && pos2 <=5); assert(pos3 >=0 && pos3 <=5);
        assert(pos4 >=0 && pos4 <=5); assert(pos5 >=0 && pos5 <=5);

        tri_ien[ebcid][6*ee+0] = node_t[pos0];
        tri_ien[ebcid][6*ee+1] = node_t[pos1];
        tri_ien[ebcid][6*ee+2] = node_t[pos2];
        tri_ien[ebcid][6*ee+3] = node_t[pos3];
        tri_ien[ebcid][6*ee+4] = node_t[pos4];
        tri_ien[ebcid][6*ee+5] = node_t[pos5];
      }
      delete tetcell;
    }
  }
  else SYS_T::print_fatal("Error: unknown element type.\n");
}

// EOF
