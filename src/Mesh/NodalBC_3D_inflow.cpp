#include "NodalBC_3D_inflow.hpp"

NodalBC_3D_inflow::NodalBC_3D_inflow( const std::vector<std::string> &inffileList,
    const std::string &wallfile,
    const int &nFunc,
    const std::vector<Vector_3> &in_outnormal,
    const int &elemtype ) 
: num_nbc( static_cast<int>( inffileList.size() ) ), elem_type( elemtype )
{
  init( inffileList, wallfile, nFunc, in_outnormal, elemtype );
}

void NodalBC_3D_inflow::init( const std::vector<std::string> &inffileList,
    const std::string &wallfile,
    const int &nFunc,
    const std::vector<Vector_3> &in_outnormal,
    const int &elemtype )
{ 
  // 1. Clear the container for Dirichlet nodes
  dir_nodes.clear();
  num_dir_nodes = 0;

  dir_nodes_on_inlet.resize( num_nbc );
  for(int ii=0; ii<num_nbc; ++ii) dir_nodes_on_inlet[ii].clear();

  num_dir_nodes_on_inlet.resize( num_nbc );

  // 2. Analyze the file type and read in the data
  num_node.resize( num_nbc );
  num_cell.resize( num_nbc );
  nLocBas.resize(  num_nbc );

  sur_ien.resize(num_nbc );
  pt_xyz.resize( num_nbc );

  global_node.resize(  num_nbc );
  global_cell.resize(  num_nbc );

  centroid.resize(  num_nbc );
  outnormal.resize( num_nbc );

  outline_pts.resize(    num_nbc );
  num_out_bc_pts.resize( num_nbc );
  outline_pts_loc_id.resize( num_nbc );

  inf_active_area.resize( num_nbc );
  face_area.resize(       num_nbc );

  intNA.resize( num_nbc ); 

  // Read the wall file
  SYS_T::file_check(wallfile);

  const std::vector<int> wall_gnode = VTK_T::read_int_PointData(wallfile, "GlobalNodeID");

  // Loop over each surface with id ii
  for( int ii=0; ii<num_nbc; ++ii )
  {
    SYS_T::file_check( inffileList[ii] );

    VTK_T::read_grid( inffileList[ii], num_node[ii], num_cell[ii], pt_xyz[ii], sur_ien[ii] );

    if( elemtype == 501 )
      nLocBas[ii] = 3;
    else if( elemtype == 502 )
      nLocBas[ii] = 6;
    else if( elemtype == 601 )
      nLocBas[ii] = 4;
    else if( elemtype == 602 )
      nLocBas[ii] = 9;
    else 
      SYS_T::print_fatal("Error: NodalBC_3D_inflow::init function: unknown element type.\n");

    global_node[ii] = VTK_T::read_int_PointData(inffileList[ii], "GlobalNodeID");
    global_cell[ii] = VTK_T::read_int_CellData(inffileList[ii], "GlobalElementID");

    // Generate the dir-node list. Nodes belonging to the wall are excluded.
    for(unsigned int jj=0; jj<global_node[ii].size(); ++jj)
    {
      SYS_T::print_fatal_if( global_node[ii][jj]<0, "Error: NodalBC_3D_inflow::init function: negative nodal index! \n");

      if( !VEC_T::is_invec( wall_gnode, global_node[ii][jj]) )
      {
        dir_nodes.push_back( global_node[ii][jj] );
        dir_nodes_on_inlet[ii].push_back( global_node[ii][jj] );
      }
    }

    num_dir_nodes_on_inlet[ii] = dir_nodes_on_inlet[ii].size();

    // Calculate the centroid of the surface
    centroid[ii].gen_zero();
    for(int jj=0; jj<num_node[ii]; ++jj)
    {
      centroid[ii](0) += pt_xyz[ii][3*jj+0];
      centroid[ii](1) += pt_xyz[ii][3*jj+1];
      centroid[ii](2) += pt_xyz[ii][3*jj+2];
    }
    centroid[ii] *= (1.0 / (double) num_node[ii]);

    // assign outward normal vector from the input
    outnormal[ii] = in_outnormal[ii];

    // Collect nodes that belong to the wall, and set up a vector that
    // is 1 on the interior nodes and 0 on the wall bc nodes.
    outline_pts[ii].clear();

    num_out_bc_pts[ii] = 0;

    double * temp_sol = new double [num_node[ii]];  

    for(int jj=0; jj<num_node[ii]; ++jj)
    {
      // If the node is not on the wall, it is an interior node, so set
      // the element to 1.
      if( !VEC_T::is_invec(wall_gnode, global_node[ii][jj]) ) 
        temp_sol[jj] = 1.0;

      // otherwise, the node is on the wall surface, so set element to 0.
      else 
      {
        temp_sol[jj] = 0.0;

        // Also store the point's coordinates in outline points
        num_out_bc_pts[ii] += 1;
        outline_pts[ii].push_back( pt_xyz[ii][3*jj+0] );
        outline_pts[ii].push_back( pt_xyz[ii][3*jj+1] );
        outline_pts[ii].push_back( pt_xyz[ii][3*jj+2] );

        outline_pts_loc_id[ii].push_back( jj );
      }
    }

    // If the number of surface nodes matches num_out_bc_pts, the wall
    // mesh contains the inlet surface. This is a common error when
    // adopting sv files, where the user uses the combined exterior surface
    // as the wall mesh. We will throw an error message if detected.
    SYS_T::print_fatal_if( num_out_bc_pts[ii] == num_node[ii], "Error: NodalBC_3D_inflow::init function:"
                          "the number of outline points is %d and the number of total points on the surface is %d."
                          "This is likely due to an improper wall mesh. \n", num_out_bc_pts[ii], num_node[ii] );

    inf_active_area[ii] = 0.0;
    face_area[ii] = 0.0;

    intNA[ii].resize( num_node[ii] );

    // zero the container
    for(int jj=0; jj<num_node[ii]; ++jj) intNA[ii][jj] = 0.0;

    if( elemtype == 501 )
    {
      const int nqp_tri = 3;                       // num qua points
      QuadPts_Gauss_Triangle quad( nqp_tri );      // quadrature rule
      FEAElement_Triangle3_3D_der0 ele( nqp_tri ); // element

      for(int ee=0; ee<num_cell[ii]; ++ee)
      {
        double eptx[3]{}; double epty[3]{}; double eptz[3]{};
        int node_idx[3]{}; 

        for(int jj=0; jj<3; ++jj)
        {
          node_idx[jj] = sur_ien[ii][3*ee+jj];
          eptx[jj] = pt_xyz[ii][ 3*node_idx[jj]+0 ];
          epty[jj] = pt_xyz[ii][ 3*node_idx[jj]+1 ];
          eptz[jj] = pt_xyz[ii][ 3*node_idx[jj]+2 ];
        }

        ele.buildBasis(&quad, eptx, epty, eptz);

        for(int qua=0; qua<nqp_tri; ++qua)
        {
          double R[3]{};
          ele.get_R( qua, R );

          const double gwts = ele.get_detJac(qua) * quad.get_qw(qua);

          for(int jj=0; jj<3; ++jj)
          {
            inf_active_area[ii] += gwts * R[jj] * temp_sol[ sur_ien[ii][3*ee+jj] ];
            face_area[ii] += gwts * R[jj];

            intNA[ii][node_idx[jj]] += gwts * R[jj];
          }
        } // end qua-loop
      } // end ee-loop
    }
    else if( elemtype == 502 )
    {
      const int nqp_tri = 6;                       // num qua points
      QuadPts_Gauss_Triangle quad( nqp_tri );      // quadrature rule
      FEAElement_Triangle6_3D_der0 ele( nqp_tri ); // element

      for(int ee=0; ee<num_cell[ii]; ++ee)
      {
        double eptx[6]{}; double epty[6]{}; double eptz[6]{}; 
        int node_idx[6]{};

        for(int jj=0; jj<6; ++jj)
        {
          node_idx[jj] = sur_ien[ii][6*ee+jj];
          eptx[jj] = pt_xyz[ii][ 3*node_idx[jj]+0 ];
          epty[jj] = pt_xyz[ii][ 3*node_idx[jj]+1 ];
          eptz[jj] = pt_xyz[ii][ 3*node_idx[jj]+2 ];
        }

        ele.buildBasis(&quad, eptx, epty, eptz);

        for(int qua=0; qua<nqp_tri; ++qua)
        {
          double R[6]{};
          ele.get_R( qua, R );

          const double gwts = ele.get_detJac(qua) * quad.get_qw(qua);

          for(int jj=0; jj<6; ++jj)
          {
            inf_active_area[ii] += gwts * R[jj] * temp_sol[ sur_ien[ii][6*ee+jj] ];
            face_area[ii] += gwts * R[jj];

            intNA[ii][node_idx[jj]] += gwts * R[jj];
          }
        } // end qua-loop
      } // end ee-loop
    }
    else if( elemtype == 601 )
    {
      const int nqp_quad = 2;                              // num qua points
      QuadPts_Gauss_Quad quad( nqp_quad );                 // quadrature rule
      FEAElement_Quad4_3D_der0 ele( nqp_quad * nqp_quad ); // element

      for(int ee=0; ee<num_cell[ii]; ++ee)
      {
        double eptx[4]{}; double epty[4]{}; double eptz[4]{}; 
        int node_idx[4]{};

        for(int jj=0; jj<4; ++jj)
        {
          node_idx[jj] = sur_ien[ii][4*ee+jj];
          eptx[jj] = pt_xyz[ii][ 3*node_idx[jj]+0 ];
          epty[jj] = pt_xyz[ii][ 3*node_idx[jj]+1 ];
          eptz[jj] = pt_xyz[ii][ 3*node_idx[jj]+2 ];
        }

        ele.buildBasis(&quad, eptx, epty, eptz);

        for(int qua=0; qua<nqp_quad * nqp_quad; ++qua)
        {
          double R[4]{};
          ele.get_R( qua, R );

          const double gwts = ele.get_detJac(qua) * quad.get_qw(qua);

          for(int jj=0; jj<4; ++jj)
          {
            inf_active_area[ii] += gwts * R[jj] * temp_sol[ sur_ien[ii][4*ee+jj] ];
            face_area[ii] += gwts * R[jj];

            intNA[ii][node_idx[jj]] += gwts * R[jj];
          }
        } // end qua-loop
      } // end ee-loop
    }
    else if( elemtype == 602 )
    {
      const int nqp_quad = 4;                              // num qua points
      QuadPts_Gauss_Quad quad( nqp_quad );                 // quadrature rule
      FEAElement_Quad9_3D_der0 ele( nqp_quad * nqp_quad ); // element

      for(int ee=0; ee<num_cell[ii]; ++ee)
      {
        double eptx[9]{}; double epty[9]{}; double eptz[9]{}; 
        int node_idx[9]{};

        for(int jj=0; jj<9; ++jj)
        {
          node_idx[jj] = sur_ien[ii][9*ee+jj];
          eptx[jj] = pt_xyz[ii][ 3*node_idx[jj]+0 ];
          epty[jj] = pt_xyz[ii][ 3*node_idx[jj]+1 ];
          eptz[jj] = pt_xyz[ii][ 3*node_idx[jj]+2 ];
        }

        ele.buildBasis(&quad, eptx, epty, eptz);

        for(int qua=0; qua<nqp_quad * nqp_quad; ++qua)
        {
          double R[9]{};
          ele.get_R( qua, R );

          const double gwts = ele.get_detJac(qua) * quad.get_qw(qua);

          for(int jj=0; jj<9; ++jj)
          {
            inf_active_area[ii] += gwts * R[jj] * temp_sol[ sur_ien[ii][9*ee+jj] ];
            face_area[ii] += gwts * R[jj];

            intNA[ii][node_idx[jj]] += gwts * R[jj];
          }
        } // end qua-loop
      } // end ee-loop
    }
    else SYS_T::print_fatal("Error: NodalBC_3D_inflow::init function: unknown element type.\n");

    delete [] temp_sol; temp_sol = nullptr;
  } // end ii-loop

  num_dir_nodes = dir_nodes.size();

  VEC_T::sort_unique_resize(dir_nodes);

  SYS_T::print_fatal_if( num_dir_nodes != dir_nodes.size(), "Error: NodalBC_3D_inflow::init function: there are repeated nodes in the inflow file list.\n" );

  // Generate ID array
  Create_ID( nFunc );

  // Finish and print info on screen
  std::cout<<"===> NodalBC_3D_inflow specified by\n";
  for(int ii=0; ii<num_nbc; ++ii)
  {
    std::cout<<"     nbc_id = "<<ii<<": "<<inffileList[ii]<<" with nodes on ";
    std::cout<<"     "<<wallfile<<" excluded.\n";
    std::cout<<"          num_node: "<<num_node[ii]<<", num_cell: "<<num_cell[ii]<<'\n';
    std::cout<<"          centroid: "<<centroid[ii](0)<<'\t'<<centroid[ii](1)<<'\t'<<centroid[ii](2)<<'\n';
    std::cout<<"          number of outline points is "<<num_out_bc_pts[ii]<<'\n';
    std::cout<<"          outward normal is ["<<outnormal[ii](0)<<'\t'<<outnormal[ii](1)<<'\t'<<outnormal[ii](2)<<"]. \n";
    std::cout<<"          area is "<<face_area[ii]<<", and active area is "<<inf_active_area[ii]<<'\n';
  }
}

void NodalBC_3D_inflow::resetSurIEN_outwardnormal( const IIEN * const &VIEN )
{
  if(elem_type == 501)
    reset501IEN_outwardnormal(VIEN); 
  else if(elem_type == 502)
    reset502IEN_outwardnormal(VIEN); 
  else if(elem_type == 601)
    reset601IEN_outwardnormal(VIEN);     
  else if(elem_type == 602)
    reset602IEN_outwardnormal(VIEN);  
  else SYS_T::print_fatal("Error: NodalBC_3D_inflow::resetSurIEN_outwardnormal function: unknown element type.\n");
}

void NodalBC_3D_inflow::reset501IEN_outwardnormal( const IIEN * const &VIEN )
{
  for(int nbcid=0; nbcid<num_nbc; ++nbcid)
  {
    TET_T::Tet4 * tetcell = new TET_T::Tet4();

    for(int ee=0; ee<num_cell[nbcid]; ++ee)
    {
      // Triangle mesh node index
      const int node_t[3] { get_ien(nbcid, ee, 0), get_ien(nbcid, ee, 1), get_ien(nbcid, ee, 2) };

      // The triangle mesh node's volumetric index
      const std::vector<int> node_t_gi { get_global_node(nbcid, node_t[0]), 
        get_global_node(nbcid, node_t[1]), get_global_node(nbcid, node_t[2]) };

      // cell ee's global/volumetric index  
      const int cell_gi = get_global_cell(nbcid, ee);

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
          SYS_T::print_fatal("Error: NodalBC_3D_inflow::reset501IEN_outwardnormal function: tet_face_id is out of range. \n");
          break;
      }
      ASSERT(pos0 >=0 && pos0 <=2, "Error: NodalBC_3D_inflow::reset501IEN_outwardnormal function logical error.\n" );
      ASSERT(pos1 >=0 && pos1 <=2, "Error: NodalBC_3D_inflow::reset501IEN_outwardnormal function logical error.\n" );
      ASSERT(pos2 >=0 && pos2 <=2, "Error: NodalBC_3D_inflow::reset501IEN_outwardnormal function logical error.\n" ); 

      // Now we have got the corrected ordering of node_t, put them back into
      // sur_ien.
      sur_ien[nbcid][3*ee+0] = node_t[pos0];
      sur_ien[nbcid][3*ee+1] = node_t[pos1];
      sur_ien[nbcid][3*ee+2] = node_t[pos2];
    }
    delete tetcell; 
  }
}

void NodalBC_3D_inflow::reset502IEN_outwardnormal( const IIEN * const &VIEN )
{
  for(int nbcid=0; nbcid<num_nbc; ++nbcid)
  {
    std::vector<int> node_t(6, 0);    // triange node index in 2D mesh
    std::vector<int> node_t_gi(6, 0); // triange node index in 3D mesh
    std::vector<int> tet_n(10, 0);    // tet node index in 3D mesh

    TET_T::Tet4 * tetcell = new TET_T::Tet4();

    for(int ee=0; ee<num_cell[nbcid]; ++ee)
    {
      for(int ii=0; ii<6; ++ii)
      {
        node_t[ii] = get_ien(nbcid, ee, ii);
        node_t_gi[ii] = get_global_node(nbcid, node_t[ii]);
      }

      const int cell_gi = get_global_cell(nbcid, ee);

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
          SYS_T::print_fatal("Error: NodalBC_3D_inflow::reset502IEN_outwardnormal function: tet_face_id is out of range. \n");
          break;
      }
      ASSERT(pos0 >=0 && pos0 <=5, "Error: NodalBC_3D_inflow::reset502IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos1 >=0 && pos1 <=5, "Error: NodalBC_3D_inflow::reset502IEN_outwardnormal function logical error.\n" );
      ASSERT(pos2 >=0 && pos2 <=5, "Error: NodalBC_3D_inflow::reset502IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos3 >=0 && pos3 <=5, "Error: NodalBC_3D_inflow::reset502IEN_outwardnormal function logical error.\n" );
      ASSERT(pos4 >=0 && pos4 <=5, "Error: NodalBC_3D_inflow::reset502IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos5 >=0 && pos5 <=5, "Error: NodalBC_3D_inflow::reset502IEN_outwardnormal function logical error.\n" );

      sur_ien[nbcid][6*ee+0] = node_t[pos0];
      sur_ien[nbcid][6*ee+1] = node_t[pos1];
      sur_ien[nbcid][6*ee+2] = node_t[pos2];
      sur_ien[nbcid][6*ee+3] = node_t[pos3];
      sur_ien[nbcid][6*ee+4] = node_t[pos4];
      sur_ien[nbcid][6*ee+5] = node_t[pos5];
    }
    delete tetcell;
  }
}

void NodalBC_3D_inflow::reset601IEN_outwardnormal( const IIEN * const &VIEN )
{
  for (int nbcid=0; nbcid<num_nbc; ++nbcid)
  {
    HEX_T::Hex8 * hexcell = new HEX_T::Hex8(); 

    for (int ee=0; ee<num_cell[nbcid]; ++ee)
    {
      // Quad mesh node index
      const int node_q[4] { get_ien(nbcid, ee, 0), get_ien(nbcid, ee, 1),  get_ien(nbcid, ee, 2), get_ien(nbcid, ee, 3) };  

      // The quad mesh node's volumetric index
      const std::vector<int> node_q_gi { get_global_node(nbcid, node_q[0]), get_global_node(nbcid, node_q[1]),  
                                         get_global_node(nbcid, node_q[2]), get_global_node(nbcid, node_q[3]) }; 

      // cell ee's global/volumetric index  
      const int cell_gi = get_global_cell(nbcid, ee);

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
          SYS_T::print_fatal("Error: NodalBC_3D_inflow::reset601IEN_outwardnormal function: hex_face_id is out of range. \n");
          break;
      }
      ASSERT(pos0 >=0 && pos0 <=3, "Error: NodalBC_3D_inflow::reset601IEN_outwardnormal function logical error.\n" );
      ASSERT(pos1 >=0 && pos1 <=3, "Error: NodalBC_3D_inflow::reset601IEN_outwardnormal function logical error.\n" );
      ASSERT(pos2 >=0 && pos2 <=3, "Error: NodalBC_3D_inflow::reset601IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos3 >=0 && pos3 <=3, "Error: NodalBC_3D_inflow::reset601IEN_outwardnormal function logical error.\n" ); 

      // Now we have got the corrected ordering of node_q, put them back into
      // sur_ien.
      sur_ien[nbcid][4*ee+0] = node_q[pos0];
      sur_ien[nbcid][4*ee+1] = node_q[pos1];
      sur_ien[nbcid][4*ee+2] = node_q[pos2];
      sur_ien[nbcid][4*ee+3] = node_q[pos3];
    }
    delete hexcell; 
  }
}

void NodalBC_3D_inflow::reset602IEN_outwardnormal( const IIEN * const &VIEN )
{
  for (int nbcid=0; nbcid<num_nbc; ++nbcid)
  {
    std::vector<int> node_q(9, 0);    // biquadratic quadrangle node index in 2D mesh
    std::vector<int> node_q_gi(9, 0); // biquadratic quadrangle node index in 3D mesh
    std::vector<int> hex_n(27, 0);    // triquadratic hex node index in 3D mesh

    HEX_T::Hex8 * hexcell = new HEX_T::Hex8();  

    for (int ee=0; ee<num_cell[nbcid]; ++ee)
    {
      for (int ii=0; ii<9; ++ii)
      {
        node_q[ii]    = get_ien(nbcid, ee, ii);
        node_q_gi[ii] = get_global_node(nbcid, node_q[ii]);
      }

      const int cell_gi = get_global_cell(nbcid, ee);

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
          SYS_T::print_fatal("Error: NodalBC_3D_inflow::reset602IEN_outwardnormal function: hex_face_id is out of range. \n");
          break;
      }
      ASSERT(pos0 >=0 && pos0 <=8, "Error: NodalBC_3D_inflow::reset602IEN_outwardnormal function logical error.\n" );
      ASSERT(pos1 >=0 && pos1 <=8, "Error: NodalBC_3D_inflow::reset602IEN_outwardnormal function logical error.\n" );
      ASSERT(pos2 >=0 && pos2 <=8, "Error: NodalBC_3D_inflow::reset602IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos3 >=0 && pos3 <=8, "Error: NodalBC_3D_inflow::reset602IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos4 >=0 && pos4 <=8, "Error: NodalBC_3D_inflow::reset602IEN_outwardnormal function logical error.\n" );
      ASSERT(pos5 >=0 && pos5 <=8, "Error: NodalBC_3D_inflow::reset602IEN_outwardnormal function logical error.\n" );
      ASSERT(pos6 >=0 && pos6 <=8, "Error: NodalBC_3D_inflow::reset602IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos7 >=0 && pos7 <=8, "Error: NodalBC_3D_inflow::reset602IEN_outwardnormal function logical error.\n" ); 
      ASSERT(pos8 >=0 && pos8 <=8, "Error: NodalBC_3D_inflow::reset602IEN_outwardnormal function logical error.\n" ); 

      sur_ien[nbcid][9*ee+0] = node_q[pos0];
      sur_ien[nbcid][9*ee+1] = node_q[pos1];
      sur_ien[nbcid][9*ee+2] = node_q[pos2];
      sur_ien[nbcid][9*ee+3] = node_q[pos3];
      sur_ien[nbcid][9*ee+4] = node_q[pos4];
      sur_ien[nbcid][9*ee+5] = node_q[pos5];
      sur_ien[nbcid][9*ee+6] = node_q[pos6];
      sur_ien[nbcid][9*ee+7] = node_q[pos7];
      sur_ien[nbcid][9*ee+8] = node_q[pos8];
    }
    delete hexcell; 
  }
}

// EOF
