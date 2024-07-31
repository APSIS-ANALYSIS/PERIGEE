#include "Gmsh_FileIO.hpp"

Gmsh_FileIO::Gmsh_FileIO( const std::string &in_file_name )
: filename( in_file_name ), elem_nlocbas{{ 0, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9,
    10, 27, 18, 14, 1, 8, 20, 15, 13, 9, 10, 12, 15, 15, 21, 
    4, 5, 6, 20, 35, 56 }}
{
  // Setup the file instream
  std::ifstream infile( filename.c_str(), std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;

  // First four lines are Gmsh default file format 
  getline(infile, sline); 
  SYS_T::print_fatal_if(sline.compare("$MeshFormat") != 0, 
      "Error: .msh format first line should be $MeshFormat. \n");

  // read the node data and element data in .msh file
  getline(infile, sline);
  if (sline.compare("2.2 0 8") == 0)
    read_msh2(infile);
  else if (sline.compare("4.1 0 8") == 0)
    read_msh4(infile);
  else
    SYS_T::print_fatal("Error: .msh format second line should be '2.2 0 8' or '4.1 0 8'. \n");

  // Finish the file reading, the last line should be $EndElements  
  getline(infile, sline);
  SYS_T::print_fatal_if(sline.compare("$EndElements") != 0, 
      "Error: .msh format line should be $EndElements. \n");

  // When periodic boundary condition is applied in .geo file
  getline(infile, sline);
  if (sline.compare("$Periodic") == 0)
    read_periodic(infile);
  else
    ; // Do nothing, just finish the reading

  // Generate the details of the 1d, 2d and 3d info
  num_phy_domain_3d = 0; num_phy_domain_2d = 0; num_phy_domain_1d = 0;
  phy_3d_index.clear(); phy_2d_index.clear(); phy_1d_index.clear();
  phy_3d_name.clear(); phy_2d_name.clear(); phy_1d_name.clear();
  for(int ii=0; ii<num_phy_domain; ++ii)
  {
    if(phy_dim[ii] == 1)
    {
      num_phy_domain_1d += 1;
      phy_1d_index.push_back( phy_index[ii] );
      phy_1d_name.push_back( phy_name[ii] );
    }
    else if(phy_dim[ii] == 2)
    {
      num_phy_domain_2d += 1;
      phy_2d_index.push_back( phy_index[ii] );
      phy_2d_name.push_back( phy_name[ii] );
    }
    else
    {
      num_phy_domain_3d += 1;
      phy_3d_index.push_back( phy_index[ii] );
      phy_3d_name.push_back( phy_name[ii] );
    }
  }

  phy_3d_nElem.clear(); phy_2d_nElem.clear(); phy_1d_nElem.clear();
  phy_3d_nElem.resize(num_phy_domain_3d);
  phy_2d_nElem.resize(num_phy_domain_2d);
  phy_1d_nElem.resize(num_phy_domain_1d);

  for(int ii=0; ii<num_phy_domain_1d; ++ii)
    phy_1d_nElem[ ii ] = phy_domain_nElem[ phy_1d_index[ii] ];

  for(int ii=0; ii<num_phy_domain_2d; ++ii)
    phy_2d_nElem[ ii ] = phy_domain_nElem[ phy_2d_index[ii] ];

  for(int ii=0; ii<num_phy_domain_3d; ++ii)
    phy_3d_nElem[ ii ] = phy_domain_nElem[ phy_3d_index[ii] ];

  phy_3d_start_index.clear();
  if(num_phy_domain_3d > 0)
  {
    phy_3d_start_index.resize(num_phy_domain_3d);
    phy_3d_start_index[0] = 0;
    for(int ii=1; ii<num_phy_domain_3d; ++ii)
      phy_3d_start_index[ii] = phy_3d_start_index[ii-1] + phy_domain_nElem[ phy_3d_index[ii-1] ];
  }

  phy_2d_start_index.clear();
  if(num_phy_domain_2d > 0)
  {
    phy_2d_start_index.resize(num_phy_domain_2d);
    phy_2d_start_index[0] = 0;
    for(int ii=1; ii<num_phy_domain_2d; ++ii)
      phy_2d_start_index[ii] = phy_2d_start_index[ii-1] + phy_domain_nElem[ phy_2d_index[ii-1] ];
  }
}

Gmsh_FileIO::~Gmsh_FileIO()
{}

void Gmsh_FileIO::print_info() const
{
  std::cout<<"File name : "<<filename<<std::endl;
  std::cout<<"=== Physical domain number: "<<num_phy_domain<<std::endl;
  std::cout<<"    Physical 1d domain number: "<<num_phy_domain_1d<<std::endl;
  std::cout<<"    physical tag: ";
  VEC_T::print(phy_1d_index);
  std::cout<<"    names: ";
  VEC_T::print(phy_1d_name);
  std::cout<<"    nElem: ";
  VEC_T::print(phy_1d_nElem);
  std::cout<<"    etype: ";
  for(int ii=0; ii<num_phy_domain_1d; ++ii)
    std::cout<<ele_type[ phy_1d_index[ii] ]<<'\t';
  std::cout<<std::endl;
  std::cout<<"    nLocBas: ";
  for(int ii=0; ii<num_phy_domain_1d; ++ii)
    std::cout<<ele_nlocbas[ phy_1d_index[ii] ]<<'\t';
  std::cout<<std::endl<<std::endl;

  std::cout<<"    Physical 2d domain number: "<<num_phy_domain_2d<<std::endl;
  std::cout<<"    physical tag: ";
  VEC_T::print(phy_2d_index);
  std::cout<<"    names: ";
  VEC_T::print(phy_2d_name);
  std::cout<<"    nElem: ";
  VEC_T::print(phy_2d_nElem);
  std::cout<<"    etype: ";
  for(int ii=0; ii<num_phy_domain_2d; ++ii)
    std::cout<<ele_type[ phy_2d_index[ii] ]<<'\t';
  std::cout<<std::endl;
  std::cout<<"    nLocBas: ";
  for(int ii=0; ii<num_phy_domain_2d; ++ii)
    std::cout<<ele_nlocbas[ phy_2d_index[ii] ]<<'\t';
  std::cout<<std::endl<<std::endl;
  
  std::cout<<"    Physical 3d domain number: "<<num_phy_domain_3d<<std::endl;
  std::cout<<"    physical tag: ";
  VEC_T::print(phy_3d_index);
  std::cout<<"    names: ";
  VEC_T::print(phy_3d_name);
  std::cout<<"    nElem: ";
  VEC_T::print(phy_3d_nElem);
  std::cout<<"    etype: ";
  for(int ii=0; ii<num_phy_domain_3d; ++ii)
    std::cout<<ele_type[ phy_3d_index[ii] ]<<'\t';
  std::cout<<std::endl;
  std::cout<<"    nLocBas: ";
  for(int ii=0; ii<num_phy_domain_3d; ++ii)
    std::cout<<ele_nlocbas[ phy_3d_index[ii] ]<<'\t';
  std::cout<<std::endl;
  
  std::cout<<"=== Total node number : "<<num_node<<std::endl;
  std::cout<<"=== Total element number : "<<num_elem<<std::endl;
}

void Gmsh_FileIO::write_interior_vtp( const std::string &vtp_filename, 
  const int &index_sur,const int &index_vol1, const int &index_vol2 ) const
{
  SYS_T::print_fatal_if( index_sur >= num_phy_domain_2d || index_sur < 0,
      "Error: Gmsh_FileIO::write_vtp, surface index is wrong. \n");

  SYS_T::print_fatal_if( index_vol1 >= num_phy_domain_3d || index_vol1 < 0,
      "Error: Gmsh_FileIO::write_vtp, volume 1 index is wrong. \n");

  SYS_T::print_fatal_if( index_vol2 >= num_phy_domain_3d || index_vol2 < 0,
      "Error: Gmsh_FileIO::write_vtp, volume 2 index is wrong. \n");

  std::cout<<"=== Gmsh_FileIO::write_interior_vtp for "
    <<phy_2d_name[index_sur]
    <<" associated with "<<phy_3d_name[index_vol1]
    <<" and "<<phy_3d_name[index_vol2]<<std::endl;

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  std::cout<<"-----> write "<<vtp_filename<<".vtp \n";
  mytimer->Reset();
  mytimer->Start();

  // Obtain the surface's IEN array associated with global nodal indices
  const int phy_index_sur = phy_2d_index[index_sur];
  std::vector<int> sur_ien_global( eIEN[phy_index_sur] );

  // obtain the number of local basis function of the surface domain
  const int nlocbas_2d {ele_nlocbas[phy_index_sur]};

  const int bcnumcl = phy_2d_nElem[index_sur];
  SYS_T::print_fatal_if( int(sur_ien_global.size() ) != nlocbas_2d * bcnumcl,
      "Error: Gmsh_FileIO::write_interior_vtp, sur IEN size wrong. \n" );

  // Obtain the volumetric IEN array
  const int phy_index_vol_1 = phy_3d_index[index_vol1];
  std::vector<int> vol_IEN_1( eIEN[phy_index_vol_1] );
  const int nlocbas_3d_1 {ele_nlocbas[phy_index_vol_1]};
  const int numcel_1 = phy_3d_nElem[index_vol1];
  SYS_T::print_fatal_if( VEC_T::get_size(vol_IEN_1) != nlocbas_3d_1 * numcel_1,
      "Error: Gmsh_FileIO::write_interior_vtp, vol_1 IEN size wrong. \n");

  const int phy_index_vol_2 = phy_3d_index[index_vol2];
  std::vector<int> vol_IEN_2( eIEN[phy_index_vol_2] );
  const int nlocbas_3d_2 {ele_nlocbas[phy_index_vol_2]};
  const int numcel_2 = phy_3d_nElem[index_vol2];
  SYS_T::print_fatal_if( VEC_T::get_size(vol_IEN_2) != nlocbas_3d_2 * numcel_2,
      "Error: Gmsh_FileIO::write_interior_vtp, vol_2 IEN size wrong. \n");

  std::string ele_2d {};
  std::string ele_3d {};
  if (nlocbas_2d == 3 && nlocbas_3d_1 == 4 && nlocbas_3d_2 == 4)
  {
    ele_2d = static_cast<std::string>("triangle");
    ele_3d = static_cast<std::string>("tetrahedron");
  }
  else if (nlocbas_2d == 4 && nlocbas_3d_1 == 8 && nlocbas_3d_2 == 8)
  {
    ele_2d = static_cast<std::string>("quadrilateral");
    ele_3d = static_cast<std::string>("hexahedron");
  }
  else
    SYS_T::print_fatal("Error: Gmsh_FileIO::write_interior_vtp, element types of surface and volume donnot match. \n");

  // bcpt stores the global nodal indices in ascending order
  std::vector<int> bcpt( sur_ien_global );
  VEC_T::sort_unique_resize( bcpt ); // unique ascending order nodes
  const int bcnumpt = VEC_T::get_size(bcpt);
  std::cout<<"      num of bc pt = "<<bcnumpt<<'\n';

  // sur_pt stores the xyz coordinates of the points for the surface
  std::vector<double> sur_pt (3*bcnumpt, 0.0);
  for( int ii=0; ii<bcnumpt; ++ii )
  {
    sur_pt[ii*3]   = node[bcpt[ii]*3] ;
    sur_pt[ii*3+1] = node[bcpt[ii]*3+1] ;
    sur_pt[ii*3+2] = node[bcpt[ii]*3+2] ;
  }

  // generate the local surface IEN array
  std::vector<int> sur_ien {};
  for(int ee=0; ee<bcnumcl; ++ee)
  {
    for (int ii{0}; ii < nlocbas_2d; ++ii)
      sur_ien.push_back( VEC_T::get_pos(bcpt, sur_ien_global[nlocbas_2d * ee + ii]) );
  }
  std::cout<<"      "<<ele_2d<<" IEN generated. \n";

  // generate a mapper for the global nodes that returns 1 if the nodes 
  // belong to the surface, 0 otherwise.
  bool * bcmap = new bool [num_node];
  for(int ii=0; ii<num_node; ++ii) bcmap[ii] = 0;
  for(int ii=0; ii<bcnumpt; ++ii) bcmap[bcpt[ii]] = 1;

  // use gelem_1 or 2 to store the vol elements that have face over
  // the surface mesh 
  std::vector<int> gelem_1 {};
  for( int ee=0; ee<numcel_1; ++ee )
  {
    int total = 0;
    for (int jj{0}; jj < nlocbas_3d_1; ++jj)
      total += bcmap[ vol_IEN_1[nlocbas_3d_1 * ee + jj] ];
    if(total >= nlocbas_2d)
      gelem_1.push_back(ee);
  }
  std::cout<<"      vol 1 domain: "<<gelem_1.size()<<" "<<ele_3d<<"s have more than "<<nlocbas_2d<<" points on the surface. \n";

  std::vector<int> gelem_2 {};
  for( int ee=0; ee<numcel_2; ++ee )
  {
    int total = 0;
    for (int jj{0}; jj < nlocbas_3d_2; ++jj)
      total += bcmap[ vol_IEN_2[nlocbas_3d_2 * ee + jj] ];
    if(total >= nlocbas_2d)
      gelem_2.push_back(ee);
  }
  std::cout<<"      vol 2 domain: "<<gelem_2.size()<<" "<<ele_3d<<"s have more than "<<nlocbas_2d<<" points on the surface. \n";

  delete [] bcmap; bcmap = nullptr;

  // determine the face-2-element mapping
  std::vector<int> face2elem_1( bcnumcl, -1 );
  std::vector<int> face2elem_2( bcnumcl, -1) ;
  for(int ff=0; ff<bcnumcl; ++ff)
  {
    std::vector<int> snode( nlocbas_2d, -1 );
    for (int ii{0}; ii < nlocbas_2d; ++ii)
      snode[ii] = sur_ien_global[nlocbas_2d * ff + ii];

    bool vol1_got_sur_elem = false;
    int ee = -1;
    while( !vol1_got_sur_elem && ee < VEC_T::get_size(gelem_1) - 1 )
    {
      ee += 1;
      const int vol_elem = gelem_1[ee];

      std::vector<int> vnode( nlocbas_3d_1, -1 );
      for (int jj{0}; jj < nlocbas_3d_1; ++jj)
        vnode[jj] = vol_IEN_1[nlocbas_3d_1 * vol_elem + jj];

      bool got_all_nodes = true;

        for (int ii{0}; ii < nlocbas_2d; ++ii)
        {
          const bool got_each_node = VEC_T::is_invec(vnode, snode[ii]);
          got_all_nodes = got_all_nodes && got_each_node;
        }
        vol1_got_sur_elem = got_all_nodes;
    }

    if(vol1_got_sur_elem)
      face2elem_1[ff] = gelem_1[ee] + phy_3d_start_index[index_vol1];
    else
    {
      face2elem_1[ff] = -1;
      std::cout<<"Warning: there are surface element not found in the volumetric mesh.\n";
    }

    bool vol2_got_sur_elem = false;
    ee = -1;
    while( !vol2_got_sur_elem && ee < VEC_T::get_size(gelem_2) - 1 )
    {
      ee += 1;
      const int vol_elem = gelem_2[ee];

      std::vector<int> vnode( nlocbas_3d_2, -1 );
      for (int jj{0}; jj < nlocbas_3d_2; ++jj)
        vnode[jj] = vol_IEN_2[nlocbas_3d_2 * vol_elem + jj];

      bool got_all_nodes = true;

        for (int ii{0}; ii < nlocbas_2d; ++ii)
        {
          const bool got_each_node = VEC_T::is_invec(vnode, snode[ii]);
          got_all_nodes = got_all_nodes && got_each_node;
        }
        vol2_got_sur_elem = got_all_nodes;
    }

    if(vol2_got_sur_elem)
      face2elem_2[ff] = gelem_2[ee] + phy_3d_start_index[index_vol2];
    else
    {
      face2elem_2[ff] = -1;
      std::cout<<"Warning: there are surface element not found in the volumetric mesh.\n";
    }
  }

  // Write the mesh file in vtp format
  std::vector<DataVecStr<int>> input_vtk_data {};
  input_vtk_data.push_back({bcpt, "GlobalNodeID", AssociateObject::Node});
  input_vtk_data.push_back({face2elem_1, "GlobalElementID_1", AssociateObject::Cell});
  input_vtk_data.push_back({face2elem_2, "GlobalElementID_2", AssociateObject::Cell});

  if (nlocbas_2d == 3)
  {
    TET_T::write_triangle_grid( vtp_filename, bcnumpt, bcnumcl,
      sur_pt, sur_ien, input_vtk_data );
  }
  else if (nlocbas_2d == 4)
  {
    HEX_T::write_quad_grid( vtp_filename, bcnumpt, bcnumcl,
      sur_pt, sur_ien, input_vtk_data );
  }
  else
    SYS_T::print_fatal("Error: Gmsh_FileIO::write_interior_vtp, undefined element type. \n");
  
  mytimer->Stop();
  std::cout<<"      Time taken "<<mytimer->get_sec()<<" sec. \n";
  delete mytimer;
}

void Gmsh_FileIO::write_interior_vtp( const int &index_sur, 
    const int &index_vol1, const int &index_vol2 ) const
{
  std::string vtp_file_name(phy_2d_name[index_sur]);
  vtp_file_name += "_";
  vtp_file_name += phy_3d_name[index_vol1];
  vtp_file_name += "_";
  vtp_file_name += phy_3d_name[index_vol2];
  
  write_interior_vtp(vtp_file_name, index_sur, index_vol1, index_vol2);
}

void Gmsh_FileIO::write_vtp( const std::string &vtp_filename,
  const int &index_sur, const int &index_vol, const bool &isf2e, const bool &is_slave ) const
{
  SYS_T::print_fatal_if( index_sur >= num_phy_domain_2d || index_sur < 0,
      "Error: Gmsh_FileIO::write_vtp, surface index is wrong. \n");

  SYS_T::print_fatal_if( index_vol >= num_phy_domain_3d || index_vol < 0,
      "Error: Gmsh_FileIO::write_vtp, volume index is wrong. \n");

  std::cout<<"=== Gmsh_FileIO::write_vtp for "
    <<phy_2d_name[index_sur]
    <<" associated with "<<phy_3d_name[index_vol];

  if( isf2e )
    std::cout<<" with face-to-volume element index. \n";
  else
    std::cout<<" without face-to-volume element index. \n";

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  std::cout<<"-----> write "<<vtp_filename<<".vtp \n";

  mytimer->Reset();
  mytimer->Start();

  // Copy the IEN from the whole domain, the nodal indices is from the
  // global domain indices.
  const int phy_index_sur = phy_2d_index[index_sur];
  std::vector<int> sur_ien_global( eIEN[phy_index_sur] );

  // obtain the volumetric mesh IEN array
  const int phy_index_vol = phy_3d_index[index_vol];
  std::vector<int> vol_IEN( eIEN[phy_index_vol] );

  // obtain the number of local basis function of the surface and volume domains
  const int nlocbas_2d {ele_nlocbas[phy_index_sur]};
  const int nlocbas_3d {ele_nlocbas[phy_index_vol]};

  const int bcnumcl = phy_2d_nElem[index_sur];
  SYS_T::print_fatal_if( VEC_T::get_size(sur_ien_global) != nlocbas_2d * bcnumcl,
      "Error: Gmsh_FileIO::write_vtp, sur IEN size wrong. \n" );

  const int numcel = phy_3d_nElem[index_vol];
  SYS_T::print_fatal_if( VEC_T::get_size(vol_IEN) != nlocbas_3d * numcel,
      "Error: Gmsh_FileIO::write_vtp, vol IEN size wrong. \n");

  std::string ele_2d {};
  std::string ele_3d {};
  if (nlocbas_2d == 3 && nlocbas_3d == 4)
  {
    ele_2d = static_cast<std::string>("triangle");
    ele_3d = static_cast<std::string>("tetrahedron");
  }
  else if (nlocbas_2d == 4 && nlocbas_3d == 8)
  {
    ele_2d = static_cast<std::string>("quadrilateral");
    ele_3d = static_cast<std::string>("hexahedron");
  }
  else
    SYS_T::print_fatal("Error: Gmsh_FileIO::write_vtp, element types of surface and volume donnot match. \n");

  // bcpt stores the global node index
  std::vector<int> bcpt( sur_ien_global );

  VEC_T::sort_unique_resize( bcpt ); // unique ascending order nodal id

  const int bcnumpt = VEC_T::get_size( bcpt );

  std::cout<<"      num of bc pt = "<<bcnumpt<<'\n';

  // sur_pt stores the coordinates of the boundary points
  std::vector<double> sur_pt(3*bcnumpt, 0.0);

  PERIGEE_OMP_PARALLEL_FOR
  for( int ii=0; ii<bcnumpt; ++ii )
  {
    sur_pt[ii*3]   = node[bcpt[ii]*3] ;
    sur_pt[ii*3+1] = node[bcpt[ii]*3+1] ;
    sur_pt[ii*3+2] = node[bcpt[ii]*3+2] ;
  }

  // generate the local triangle or quadrilateral IEN array
  std::vector<int> sur_ien {};
  for(int ee=0; ee<bcnumcl; ++ee)
  {
    for (int ii{0}; ii < nlocbas_2d; ++ii)
      sur_ien.push_back( VEC_T::get_pos(bcpt, sur_ien_global[nlocbas_2d * ee + ii]) );
  }
  std::cout<<"      " << ele_2d <<" IEN generated. \n";

  // determine the face-to-element mapping, if we demand it (
  // meaning this face needs boundary integral); otherwise, set
  // the face2elem as -1, since we will only need the nodal indices
  // for Dirichlet type face.
  std::vector<int> face2elem( bcnumcl, -1 );
  if( isf2e )
  {
    // generate a mapper that maps the bc node to 1, other node to 0
    bool * bcmap = new bool [num_node];
    for(int ii=0; ii<num_node; ++ii) bcmap[ii] = 0;
    for(int ii=0; ii<bcnumpt; ++ii) bcmap[bcpt[ii]] = 1;

    // use the bcmap to obtain the vol element that has its face on this surface
    std::vector<int> gelem {};

    PERIGEE_OMP_PARALLEL
    {
      std::vector<int> temp_gelem {};
      PERIGEE_OMP_FOR
      for( int ee=0; ee<numcel; ++ee )
      {
        int total = 0;
        for (int jj=0; jj < nlocbas_3d; ++jj)
          total += bcmap[ vol_IEN[nlocbas_3d  * ee + jj] ];
        if(total >= nlocbas_2d)
          temp_gelem.push_back(ee);
      }
      PERIGEE_OMP_CRITICAL
      VEC_T::insert_end(gelem, temp_gelem);
    }

    delete [] bcmap; bcmap = nullptr;
    std::cout<<"      "<<gelem.size()<<" "<<ele_3d<<"s have faces over the surface. \n";

    PERIGEE_OMP_PARALLEL_FOR
    for(int ff=0; ff<bcnumcl; ++ff)
    {
      std::vector<int> snode( nlocbas_2d, -1);
      for (int ii{0}; ii < nlocbas_2d; ++ii)
        snode[ii] = sur_ien_global[nlocbas_2d * ff + ii];

      bool got_sur_elem = false;
      int ee = -1;
      while( !got_sur_elem && ee < VEC_T::get_size(gelem) - 1 )
      {
        ee += 1;
        const int vol_elem = gelem[ee];

        std::vector<int> vnode( nlocbas_3d, -1 );
        for (int jj{0}; jj < nlocbas_3d; ++jj)
          vnode[jj] = vol_IEN[nlocbas_3d * vol_elem + jj];

        bool got_all_nodes = true;

        for (int ii{0}; ii < nlocbas_2d; ++ii)
        {
          const bool got_each_node = VEC_T::is_invec(vnode, snode[ii]);
          got_all_nodes = got_all_nodes && got_each_node;
        }
        got_sur_elem = got_all_nodes;
      }

      // If the boundary surface element is not found, 
      // we write -1 as the mapping value
      if(got_sur_elem)
        face2elem[ff] = gelem[ee] + phy_3d_start_index[index_vol];
      else
        face2elem[ff] = -1;
    }
    std::cout<<"      face2elem mapping generated. \n";
  }

  // Write the mesh file in vtp format
  std::vector<DataVecStr<int>> input_vtk_data {};
  input_vtk_data.push_back({bcpt, "GlobalNodeID", AssociateObject::Node});
  input_vtk_data.push_back({face2elem, "GlobalElementID", AssociateObject::Cell});

  // Write the master nodes' global id
  if (is_slave)
  {
    std::vector<int> master_id ( bcnumpt, -1 );
    for(int ii{0}; ii < bcnumpt; ++ii)
    {
      // The position of slave node in per_slave vector
      const int pos_slave = VEC_T::get_pos(per_slave, bcpt[ii]);
      SYS_T::print_fatal_if( pos_slave == -1,
        "Error: Gmsh_FileIO::write_vtp, node %d of boundary %s is not a slave node.\n", bcpt[ii], phy_2d_name[index_sur].c_str());
      
      master_id[ii] = per_master[pos_slave];
    }
    std::cout<<"      master-slave mapping generated. \n";
    input_vtk_data.push_back({master_id, "MasterNodeID", AssociateObject::Node});
  }

  if (nlocbas_2d == 3)
  {
    TET_T::write_triangle_grid( vtp_filename, bcnumpt, bcnumcl,
      sur_pt, sur_ien, input_vtk_data );
  }
  else if (nlocbas_2d == 4)
  {
    HEX_T::write_quad_grid( vtp_filename, bcnumpt, bcnumcl,
      sur_pt, sur_ien, input_vtk_data );
  }
  else
    SYS_T::print_fatal("Error: Gmsh_FileIO::write_vtp, undefined element type. \n");

  mytimer->Stop();
  std::cout<<"      Time taken "<<mytimer->get_sec()<<" sec. \n";
  delete mytimer;
}

void Gmsh_FileIO::write_vtp(const std::string &vtp_filename,
  const std::string &phy_name_sur, const std::string &phy_name_vol,
  const bool &isf2e, const bool &is_slave) const
{
  const int index_sur = VEC_T::get_pos(phy_2d_name, phy_name_sur);
  SYS_T::print_fatal_if(index_sur == -1,
    "Error: Gmsh_FileIO::write_vtp, wrong physical name of surface.\n");

  const int index_vol = VEC_T::get_pos(phy_3d_name, phy_name_vol);
  SYS_T::print_fatal_if(index_vol == -1,
    "Error: Gmsh_FileIO::write_vtp, wrong physical name of volume.\n");

  write_vtp(vtp_filename, index_sur, index_vol, isf2e, is_slave);
}

void Gmsh_FileIO::write_each_vtu( const std::vector<std::string> name_list) const
{
  std::cout<<"=== Gmsh_FileIO::write_each_vtu.\n";
  std::cout<<"--- There are "<<num_phy_domain_3d<<" 3D physical domains.\n";

  SYS_T::print_fatal_if(VEC_T::get_size(name_list) != num_phy_domain_3d,
    "Error: Gmsh_FileIO::write_each_vtu, the number of files should match the number of 3d domains");

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  for(int ii=0; ii<num_phy_domain_3d; ++ii)
  {
    const std::string vtu_file_name = name_list[ ii ];
    std::cout<<"-----> write "<<vtu_file_name<<".vtu \t";

    mytimer->Reset();
    mytimer->Start();

    std::cout<<" nElem = "<<phy_3d_nElem[ii]<<'\t';

    const int start_eindex = phy_3d_start_index[ii];

    std::cout<<" starting e index = "<<start_eindex<<'\t';

    // generate physics tag
    std::vector<int> ptag(phy_3d_nElem[ii], ii);

    // collect the FSI indices of the sub-domain nodes
    const int domain_index = phy_3d_index[ii];
    std::vector<int> local_node_idx = eIEN[ domain_index ];
    VEC_T::sort_unique_resize( local_node_idx );

    const int num_local_node = VEC_T::get_size(local_node_idx);

    // collect those points' coordinates
    std::vector<double> local_coor( 3 * num_local_node , 0.0 );
    PERIGEE_OMP_PARALLEL_FOR
    for(int jj=0; jj<num_local_node; ++jj )
    {
      local_coor[ 3*jj+0 ] = node[ 3*local_node_idx[jj] + 0 ];
      local_coor[ 3*jj+1 ] = node[ 3*local_node_idx[jj] + 1 ];
      local_coor[ 3*jj+2 ] = node[ 3*local_node_idx[jj] + 2 ];
    }

    // generate a local IEN
    const int nloc = ele_nlocbas[domain_index];
    std::vector<int> domain_IEN( phy_3d_nElem[ii] * nloc, 0 );
    PERIGEE_OMP_PARALLEL_FOR
    for(int ee=0; ee<phy_3d_nElem[ii]; ++ee)
    {
      for(int jj=0; jj<nloc; ++jj)
      {
        const int target = eIEN[ domain_index ][ ee*nloc + jj ];
        domain_IEN[ ee * nloc + jj ] = VEC_T::get_pos( local_node_idx, target );
      }
    } 

    // Element index (using the start_eindex)
    std::vector<int> local_cell_idx( phy_3d_nElem[ii], 0 );
    
    PERIGEE_OMP_PARALLEL_FOR
    for(int ee=0; ee<phy_3d_nElem[ii]; ++ee)
      local_cell_idx[ee] = start_eindex + ee;

    // write the sub-volumetric domain's vtk/vtu file
    // the subdomain element index start with the start_eindex
    std::vector<DataVecStr<int>> input_vtk_data {};
    input_vtk_data.push_back({local_node_idx, "GlobalNodeID", AssociateObject::Node});
    input_vtk_data.push_back({local_cell_idx, "GlobalElementID", AssociateObject::Cell});
    input_vtk_data.push_back({ptag, "Physics_tag", AssociateObject::Cell});

    // Element type of this domain
    // Element type is defined by Gmsh, different from vtk->GetCellType()!
    const int Etype = ele_type[phy_3d_index[ii]];
    if ( Etype == 4 || Etype == 11 )
    {
      TET_T::write_tet_grid( vtu_file_name, num_local_node, 
        phy_3d_nElem[ ii ], local_coor, domain_IEN, input_vtk_data, true);
    }
    else if ( Etype == 5 || Etype == 12 )
    {
      HEX_T::write_hex_grid( vtu_file_name, num_local_node, 
        phy_3d_nElem[ ii ], local_coor, domain_IEN, input_vtk_data, true);
    }
    else
      SYS_T::print_fatal("Error: Gmsh_FileIO::write_each_vtu, undefined element type of domain %d. \n", phy_3d_index[ii] + 1);

    mytimer->Stop();

    std::cout<<mytimer->get_sec()<<" sec. \n";
  }

  delete mytimer;
}

void Gmsh_FileIO::write_each_vtu() const
{
  write_each_vtu(phy_3d_name);
}

void Gmsh_FileIO::write_vtu( const std::string &in_fname, 
    const bool &isXML ) const
{
  std::cout<<"=== Gmsh_FileIO::write_vtu. \n";
  std::cout<<"    There are "<<num_phy_domain_3d
    <<" 3D physical domains, with indices: \n";

  SYS_T::Timer * mytimer = new SYS_T::Timer();
  mytimer->Reset();
  mytimer->Start();

  // Prepare for the whole mesh's IEN, and physical tag
  std::vector<int> wIEN {}; // whole mesh IEN
  std::vector<int> wtag {}; // whole vol mehs phys tag
  int wnElem = 0; // whole mesh number of elements

  // whole mesh num of node is assumed to be num_node 
  const int wnNode = num_node;

  // Element type of whole mesh.
  // Now we need all 3d domain use the same element.
  const int wElemtype = ele_type[phy_3d_index[0]];

  // The 3d volumetric elements are list in the order of phy_3d_index.
  // That means, phy_3d_index[0]'s elements come first, then phy_3d_index[1],
  // etc. 
  for(int ii=0; ii<num_phy_domain_3d; ++ii)
  {
    const int domain_index = phy_3d_index[ii];
    std::cout<<"    "<<domain_index<<'\t';
    wnElem += phy_3d_nElem[ ii ];
    VEC_T::insert_end( wIEN, eIEN[domain_index] );

    for(int jj=0; jj<phy_3d_nElem[ii]; ++jj)
      wtag.push_back(ii);

    SYS_T::print_fatal_if(ele_type[phy_3d_index[ii]] != wElemtype, 
      "Error: Gmsh_FileIO::write_vtu, 3d domain use different elements. \n" );
  }

  std::cout<<"\n    "<<wnElem<<" total elems and "<<wnNode<<" total nodes. \n";

  // write whole domain
  std::vector<DataVecStr<int>> input_vtk_data {};
  input_vtk_data.push_back({wtag, "Physics_tag", AssociateObject::Cell});
  
  std::vector<int> temp_nid(wnNode, 0);
  for(int ii=0; ii<wnNode; ++ii) temp_nid[ii] = ii;
  input_vtk_data.push_back({temp_nid, "GlobalNodeID", AssociateObject::Node});

  std::vector<int> temp_eid(wnElem, 0);
  for(int ii=0; ii<wnElem; ++ii) temp_eid[ii] = ii;
  input_vtk_data.push_back({temp_eid, "GlobalElementID", AssociateObject::Cell});
  
  if ( wElemtype == 4 || wElemtype == 11 )
  {
    TET_T::write_tet_grid( in_fname, wnNode, wnElem, node,
      wIEN, input_vtk_data, isXML ); 
  }
  else if ( wElemtype == 5 || wElemtype == 12 )
  {
    HEX_T::write_hex_grid( in_fname, wnNode, wnElem, node,
      wIEN, input_vtk_data, isXML ); 
  }
  else
    SYS_T::print_fatal("Error: Gmsh_FileIO::write_vtu, undefined element type. \n" );

  mytimer->Stop();
  std::cout<<"    Time taken "<<mytimer->get_sec()<<" sec. \n";
  delete mytimer;
}

void Gmsh_FileIO::write_solid_vtu( const std::string &in_fname, 
    const bool &isXML ) const
{
  std::cout<<"=== Gmsh_FileIO::write_solid_vtu.\n";
  std::cout<<"--- There are "<<num_phy_domain_3d-1<<" 3D solid physical domains.\n";

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  mytimer->Reset();
  mytimer->Start();

  std::cout<<"-----> write solid.vtu \n";
  
  std::vector<int> stag {};
  std::vector<int> solid_node_idx {};
  int snElem = 0;
  std::vector<int> sIEN {};
  std::vector<int> solid_cell_idx {};

  for(int ii=1; ii<num_phy_domain_3d; ++ii)
  {
    const int start_eindex = phy_3d_start_index[ii];
    std::vector<int> local_cell_idx( phy_3d_nElem[ii], 0 );
    PERIGEE_OMP_PARALLEL_FOR
    for(int ee=0; ee<phy_3d_nElem[ii]; ++ee)
      local_cell_idx[ee] = start_eindex + ee;
    VEC_T::insert_end( solid_cell_idx, local_cell_idx );
    snElem += phy_3d_nElem[ii];

    // generate physics tag
    std::vector<int> ptag(phy_3d_nElem[ii], ii);
    VEC_T::insert_end( stag, ptag );

    const int domain_index = phy_3d_index[ii];
    std::vector<int> local_node_idx = eIEN[ domain_index ];
    VEC_T::sort_unique_resize( local_node_idx );
    VEC_T::insert_end( solid_node_idx, local_node_idx );

    VEC_T::insert_end(sIEN, eIEN[domain_index]);
  }

  VEC_T::sort_unique_resize( solid_node_idx );
  const int num_solid_node = VEC_T::get_size(solid_node_idx);

  std::vector<double> solid_coor( 3 * num_solid_node , 0.0 );
  PERIGEE_OMP_PARALLEL_FOR
  for(int jj=0; jj<num_solid_node; ++jj )
  {
    solid_coor[ 3*jj+0 ] = node[ 3*solid_node_idx[jj] + 0 ];
    solid_coor[ 3*jj+1 ] = node[ 3*solid_node_idx[jj] + 1 ];
    solid_coor[ 3*jj+2 ] = node[ 3*solid_node_idx[jj] + 2 ];
  }

  const int nloc = ele_nlocbas[phy_3d_index[0]];
  std::vector<int> domain_IEN( snElem * nloc, 0 );
  PERIGEE_OMP_PARALLEL_FOR
  for(int ee=0; ee<snElem; ++ee)
  {
    for(int jj=0; jj<nloc; ++jj)
    {
      const int target = sIEN[ ee*nloc + jj ];
      domain_IEN[ ee * nloc + jj ] = VEC_T::get_pos( solid_node_idx, target );
    }
  }

  // write the solid domain's vtk/vtu file
  // the subdomain element index start with the start_eindex
  std::vector<DataVecStr<int>> input_vtk_data {};
  input_vtk_data.push_back({solid_node_idx, "GlobalNodeID", AssociateObject::Node});
  input_vtk_data.push_back({solid_cell_idx, "GlobalElementID", AssociateObject::Cell});
  input_vtk_data.push_back({stag, "Physics_tag", AssociateObject::Cell});

  // Element type of this domain
  // Element type is defined by Gmsh, different from vtk->GetCellType()!
  const int Etype = ele_type[phy_3d_index[0]];
  if ( Etype == 4 || Etype == 11 )
  {
    TET_T::write_tet_grid( in_fname, num_solid_node, 
      snElem, solid_coor, domain_IEN, input_vtk_data, isXML);
  }
  else if ( Etype == 5 || Etype == 12 )
  {
    HEX_T::write_hex_grid( in_fname, num_solid_node, 
      snElem, solid_coor, domain_IEN, input_vtk_data, isXML);
  }
  else
    SYS_T::print_fatal("Error: Gmsh_FileIO::write_solid_vtu, undefined element type of domain %d. \n", phy_3d_index[0] + 1);

  mytimer->Stop();

  std::cout<<mytimer->get_sec()<<" sec. \n";
  
  delete mytimer;
}

// void Gmsh_FileIO::check_FSI_ordering( const std::string &phy1,
//     const std::string &phy2 ) const
// {
//   SYS_T::print_fatal_if( num_phy_domain_3d != 2, "Error: Gmsh_FileIO FSI mesh should contain only 2 physical domains.\n");

//   const std::string name0 = phy_3d_name[ 0 ];
//   const std::string name1 = phy_3d_name[ 1 ];

//   SYS_T::print_fatal_if( name0.compare(phy1), "Error: Gmsh_FileIO FSI mesh 3d subdomain index 0 should be fluid domain.\n" );
//   SYS_T::print_fatal_if( name1.compare(phy2), "Error: Gmsh_FileIO FSI mesh 3d subdomain index 1 should be solid domain.\n" );
// }

void Gmsh_FileIO::check_FSI_ordering( const std::string &phy1,
    const std::string &phy2, const std::string &phy3 ) const
{
  SYS_T::print_fatal_if( num_phy_domain_3d != 3, "Error: Gmsh_FileIO FSI mesh should contain 3 physical domains.\n");

  const std::string name0 = phy_3d_name[ 0 ];
  const std::string name1 = phy_3d_name[ 1 ];
  const std::string name2 = phy_3d_name[ 2 ];

  SYS_T::print_fatal_if( name0.compare(phy1), "Error: Gmsh_FileIO FSI mesh 3d subdomain index 0 should be fluid domain.\n" );
  SYS_T::print_fatal_if( name1.compare(phy2), "Error: Gmsh_FileIO FSI mesh 3d subdomain index 1 should be solid_1 domain.\n" );
  SYS_T::print_fatal_if( name2.compare(phy3), "Error: Gmsh_FileIO FSI mesh 3d subdomain index 2 should be solid_2 domain.\n" );
}

void Gmsh_FileIO::write_sur_h5( const int &index_2d, 
    const std::vector<int> &index_1d ) const
{
  // Perform basic logical checks
  SYS_T::print_fatal_if( index_2d >= num_phy_domain_2d || index_2d < 0,
      "Error: Gmsh_FileIO::write_sur_h5, surface index is wrong. \n");

  const unsigned int num_1d_edge = index_1d.size();

  for(unsigned int ii=0; ii<num_1d_edge; ++ii)
    SYS_T::print_fatal_if( index_1d[ii] >= num_phy_domain_1d || index_1d[ii] < 0,
        "Error: Gmsh_FileIO::write_sur_h5, edge index is wrong. \n");

  // Open an HDF5 file
  std::string h5_file_name( "Gmsh_" );
  h5_file_name.append( phy_2d_name[index_2d] );
  h5_file_name.append(".h5");

  std::cout<<"=== Gmsh_FileIO::write_sur_h5 for "
    <<phy_2d_name[index_2d]<<" associated with ";

  for(unsigned int ii=0; ii<num_1d_edge; ++ii)
    std::cout<<phy_1d_name[ index_1d[ii] ]<<'\t';
  std::cout<<std::endl;

  hid_t file_id = H5Fcreate(h5_file_name.c_str(), 
      H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  // Write 2D domain first
  const std::string slash("/");
  const int domain_2d_idx = phy_2d_index[index_2d];
  std::cout<<"-----> write "<<phy_2d_name[index_2d]<<'\n';

  const int num_2d_cell = phy_2d_nElem[index_2d];
  const int nLocBas_2d  = ele_nlocbas[domain_2d_idx];

  h5w->write_intScalar(file_id, "num_edge", num_1d_edge);
  h5w->write_intScalar(file_id, "num_node", num_node);
  h5w->write_intScalar(file_id, "num_cell", num_2d_cell); 
  h5w->write_intScalar(file_id, "ele_type", ele_type[domain_2d_idx]);
  h5w->write_intScalar(file_id, "nLocBas",  nLocBas_2d);
  h5w->write_intScalar(file_id, "index_2d", index_2d);
  h5w->write_intScalar(file_id, "phy_tag",  phy_2d_index[index_2d]);
  h5w->write_intScalar(file_id, "start_index_2d", phy_2d_start_index[index_2d]);
  h5w->write_doubleVector(file_id, "node", node);
  h5w->write_intVector(file_id, "IEN", eIEN[domain_2d_idx]);

  // Write 1D domain next
  for(unsigned int ii=0; ii<num_1d_edge; ++ii)
  {
    const int domain_1d_idx = phy_1d_index[ index_1d[ii] ];

    std::cout<<"-----> write /"<<phy_1d_name[ index_1d[ii] ]<<'\n';

    const int nLocBas_1d = ele_nlocbas[domain_1d_idx];

    const int num_1d_cell = phy_1d_nElem[ index_1d[ii] ];

    std::vector<int> edge_ien_global( eIEN[domain_1d_idx] );

    std::vector<int> bcpt( edge_ien_global ); 

    VEC_T::sort_unique_resize( bcpt );

    const int bcnumpt = VEC_T::get_size( bcpt );

    std::cout<<"      num of bc pt = "<<bcnumpt<<'\n';    

    // surpt stores the coordinates of the boundary points 
    std::vector<double> surpt( 3*bcnumpt, 0.0 );
    for( int jj=0; jj<bcnumpt; ++jj )
    {
      surpt[jj*3]   = node[bcpt[jj]*3] ;
      surpt[jj*3+1] = node[bcpt[jj]*3+1] ;
      surpt[jj*3+2] = node[bcpt[jj]*3+2] ;
    } 

    // generate a mapper that maps the bc node to 1; other node to 0
    bool * bcmap = new bool [num_node]; 
    for(int jj=0; jj<num_node; ++jj) bcmap[jj] = 0;
    for(int jj=0; jj<bcnumpt; ++jj) bcmap[bcpt[jj]] = 1;

    // generate a list of surface elements that has boundary on this edge 
    std::vector<int> gelem {};
    for( int ee=0; ee<num_2d_cell; ++ee )
    {
      int total = 0;
      for(int jj=0; jj<nLocBas_2d; ++jj)
        total += bcmap[ eIEN[domain_2d_idx][nLocBas_2d*ee+jj] ];

      if(total >= 2) gelem.push_back(ee); 
    } 
    delete [] bcmap; bcmap = nullptr;
    std::cout<<"      "<<gelem.size()<<" elems have edge over the boundary.\n";

    // generate the local edge element IEN array
    std::vector<int> edge_ien_local {};
    for(int ee=0; ee<num_1d_cell; ++ee)
    {
      for(int jj=0; jj<nLocBas_1d; ++jj)
        edge_ien_local.push_back( VEC_T::get_pos(bcpt, edge_ien_global[nLocBas_1d*ee+jj]) );
    }
    std::cout<<"      edge IEN generated. \n";

    // Locate the surface element that the edge belongs to. In Gmsh,
    // all elements are defined in this way. The edge two end points are
    // the first two points in the IEN array, no mater what the order is.
    // The surface element's vertices lay in the front of the surface IEN,
    // no mather how many interior points there are.
    // Hence, we only need to treat all elements as if they are linear
    // line/triangle/quadrilateral elements. Once the end points match with
    // the vertices, the edge is on the surface boundary.
    std::vector<int> face2elem(num_1d_cell, -1);
    for(int ff=0; ff<num_1d_cell; ++ff)
    {
      const int node0 = edge_ien_global[ nLocBas_1d * ff + 0 ];
      const int node1 = edge_ien_global[ nLocBas_1d * ff + 1 ];
      bool gotit = false;
      int ee = -1;
      while( !gotit && ee < int(gelem.size()) - 1 )
      {
        ee += 1;
        const int sur_elem = gelem[ee];
        int snode[3] { eIEN[domain_2d_idx][nLocBas_2d * sur_elem],
                       eIEN[domain_2d_idx][nLocBas_2d * sur_elem + 1],
                       eIEN[domain_2d_idx][nLocBas_2d * sur_elem + 2] };
        std::sort(snode, snode+3);

        const bool got0 = ( std::find(snode, snode+3, node0) != snode+3 );
        const bool got1 = ( std::find(snode, snode+3, node1) != snode+3 );
        gotit = got0 && got1;
      }
      if(gotit)
        face2elem[ff] = gelem[ee] + phy_2d_start_index[index_2d];
      else
        face2elem[ff] = -1;
    }
    std::cout<<"      face2elem mapping generated. \n"; 

    // Record info
    std::string name_1d_domain(slash);
    name_1d_domain.append( std::to_string( static_cast<int>(ii) ) );
    hid_t g_id = H5Gcreate( file_id, name_1d_domain.c_str(), 
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    h5w->write_string(g_id, "name", phy_1d_name[ index_1d[ii] ] );
    h5w->write_intScalar(g_id, "phy_tag", index_1d[ii]);
    h5w->write_intScalar(g_id, "num_cell", num_1d_cell);
    h5w->write_intScalar(g_id, "num_node", bcnumpt );
    h5w->write_intScalar(g_id, "ele_type", ele_type[domain_1d_idx]);
    h5w->write_intScalar(g_id, "nLocBas", nLocBas_1d);
    h5w->write_intVector(g_id, "IEN_glo", edge_ien_global);
    h5w->write_intVector(g_id, "IEN_loc", edge_ien_local);
    h5w->write_intVector(g_id, "pt_idx", bcpt);
    h5w->write_intVector(g_id, "edge2elem", face2elem);
    h5w->write_doubleVector(g_id, "pt_coor", surpt);

    H5Gclose(g_id);
  }

  // Close the HDF5 file
  delete h5w;
  H5Fclose( file_id );
}

void Gmsh_FileIO::write_vol_h5( const int &index_3d,
    const std::vector<int> &index_2d ) const
{
  // Perform basic index boundary check
  SYS_T::print_fatal_if( index_3d >= num_phy_domain_3d || index_3d < 0,
      "Error: Gmsh_FileIO::write_vol_h5, surface index is wrong.\n");

  const unsigned int num_2d_face = index_2d.size();

  for(unsigned int ii=0; ii<num_2d_face; ++ii)
    SYS_T::print_fatal_if( index_2d[ii] >= num_phy_domain_2d || index_2d[ii] < 0,
        "Error: Gmsh_FileIO::write_vol_h5, face index is wrong. \n");

  // Open an HDF5 file
  std::string h5_file_name( "Gmsh_" );
  h5_file_name.append( phy_3d_name[index_3d] );
  h5_file_name.append(".h5");

  std::cout<<"=== Gmsh_FileIO::write_vol_h5 for "
    <<phy_3d_name[index_3d]<<" associated with ";

  for(unsigned int ii=0; ii<num_2d_face; ++ii)
    std::cout<<phy_2d_name[ index_2d[ii] ]<<'\t';
  std::cout<<std::endl;

  hid_t file_id = H5Fcreate(h5_file_name.c_str(),
      H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  // Write the 3D domain at the root
  const std::string slash("/");
  const int domain_3d_idx = phy_3d_index[index_3d];
  std::cout<<"-----> write "<<phy_3d_name[ index_3d ]<<'\n';

  const int num_3d_cell = phy_3d_nElem[index_3d];
  const int nLocBas_3d  = ele_nlocbas[ domain_3d_idx ];

  h5w->write_intScalar(file_id, "num_face", num_2d_face);
  h5w->write_intScalar(file_id, "num_node", num_node);
  h5w->write_intScalar(file_id, "num_cell", num_3d_cell);
  h5w->write_intScalar(file_id, "ele_type", ele_type[domain_3d_idx]);
  h5w->write_intScalar(file_id, "nLocBas", nLocBas_3d);
  h5w->write_intScalar(file_id, "index_3d", index_3d);
  h5w->write_intScalar(file_id, "phy_tag", phy_3d_index[index_3d]);
  h5w->write_intScalar(file_id, "start_index_3d", phy_3d_start_index[index_3d]);
  h5w->write_doubleVector(file_id, "node", node);
  h5w->write_intVector(file_id, "IEN", eIEN[domain_3d_idx]);

  // Write the 2D domain in groups
  for(unsigned int ii=0; ii<num_2d_face; ++ii)
  {
    std::cout<<"-----> write /"<<phy_2d_name[ index_2d[ii] ]<<'\n';

    const int domain_2d_idx = phy_2d_index[ index_2d[ii] ];
    const int nLocBas_2d = ele_nlocbas[ domain_2d_idx ];
    const int num_2d_cell = phy_2d_nElem[ index_2d[ii] ];

    std::vector<int> face_ien_global( eIEN[domain_2d_idx] );
    std::vector<int> bcpt( face_ien_global );
    VEC_T::sort_unique_resize(bcpt);
    const int bcnumpt = VEC_T::get_size( bcpt );
    std::cout<<"      num of bc pt = "<<bcnumpt<<'\n';

    // tript stores the coordinates of the boundary points
    std::vector<double> surpt(3*bcnumpt, 0.0);
    for( int jj=0; jj<bcnumpt; ++jj )
    {
      surpt[jj*3]   = node[bcpt[jj]*3];
      surpt[jj*3+1] = node[bcpt[jj]*3+1];
      surpt[jj*3+2] = node[bcpt[jj]*3+2];
    }

    // Generate a mapper that maps the bc node to 1; other node to 0
    bool * bcmap = new bool [num_node];
    for(int jj=0; jj<num_node; ++jj) bcmap[jj] = 0;
    for(int jj=0; jj<bcnumpt; ++jj) bcmap[bcpt[jj]] = 1;

    int nVertex_2d {0}, nVertex_3d {0};
    if(nLocBas_2d == 3 || nLocBas_2d == 6) // tri-tet
    {
      nVertex_2d = 3;
      nVertex_3d = 4;
    }
    else if(nLocBas_2d == 4 || nLocBas_2d == 9) // quad-hex
    {
      nVertex_2d = 4;
      nVertex_3d = 8;
    }
    else
      SYS_T::print_fatal("Error: Gmsh_FileIO::write_vol_h5, unknown element type.\n");

    // Generate a list of volume elements that have boundary over the face
    std::vector<int> gelem {};
    for( int ee=0; ee<num_3d_cell; ++ee )
    {
      int total = 0;
      for(int jj=0; jj<nVertex_3d; ++jj)
        total += bcmap[ eIEN[domain_3d_idx][nLocBas_3d*ee+jj] ];

      if(total >= nVertex_2d) gelem.push_back(ee);
    }
    delete [] bcmap; bcmap = nullptr;
    std::cout<<"      "<<gelem.size()<<" elems have face over the boundary.\n";

    // Generate the local face element IEN array
    std::vector<int> face_ien_local; face_ien_local.clear();
    for(int ee=0; ee<num_2d_cell; ++ee)
    {
      for(int jj=0; jj<nLocBas_2d; ++jj)
        face_ien_local.push_back(VEC_T::get_pos(bcpt, face_ien_global[nLocBas_2d*ee+jj]));
    }
    std::cout<<"      edge IEN generated. \n";

    // Loacate the volumetric element that the face element belongs to 
    std::vector<int> face2elem(num_2d_cell, -1);
    for(int ff=0; ff<num_2d_cell; ++ff)
    { 
      std::vector<int> sur_node (nVertex_2d, -1); // the vertex indices of a surface element
      
      for(int kk {0}; kk < nVertex_2d; ++kk)
        sur_node[kk] = face_ien_global[ nLocBas_2d * ff + kk ];
      
      bool gotit = false;
      int ee = -1;
      while( !gotit && ee < int(gelem.size()) - 1 )
      {
        ee += 1;
        const int vol_elem = gelem[ee];
        std::vector<int> vol_node (nVertex_3d, -1); // the vertex indices of a volume element
        for(int jj {0}; jj < nVertex_3d; ++jj)
          vol_node[jj] = eIEN[domain_3d_idx][ nLocBas_3d * vol_elem + jj ];
        
        bool got_all_node = true;
        for(int kk {0}; kk < nVertex_2d; ++kk)
          got_all_node = got_all_node && VEC_T::is_invec(vol_node, sur_node[kk]);

        gotit = got_all_node;
      }
      if(gotit)
        face2elem[ff] = gelem[ee] + phy_3d_start_index[index_3d];
      else
        face2elem[ff] = -1;
    }
    std::cout<<"      face2elem mapping generated.\n";

    // Record info
    std::string name_2d_domain(slash);
    name_2d_domain.append( std::to_string( static_cast<int>(ii) ) );
    hid_t g_id = H5Gcreate( file_id, name_2d_domain.c_str(),
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    h5w->write_string(g_id, "name", phy_2d_name[ index_2d[ii] ] );
    h5w->write_intScalar(g_id, "phy_tag", index_2d[ii]);
    h5w->write_intScalar(g_id, "num_cell", num_2d_cell);
    h5w->write_intScalar(g_id, "num_node", bcnumpt );
    h5w->write_intScalar(g_id, "ele_type", ele_type[domain_2d_idx]);
    h5w->write_intScalar(g_id, "nLocBas", nLocBas_2d);
    h5w->write_intVector(g_id, "IEN_glo", face_ien_global);
    h5w->write_intVector(g_id, "IEN_loc", face_ien_local);
    h5w->write_intVector(g_id, "pt_idx", bcpt);
    h5w->write_intVector(g_id, "face2elem", face2elem);
    h5w->write_doubleVector(g_id, "pt_coor", surpt);

    H5Gclose(g_id);
  } // End-loop-over-2d-face

  // Close the HDF5 file
  delete h5w;
  H5Fclose( file_id ); 
}

void Gmsh_FileIO::write_vol_h5( const int &index_3d,
    const std::vector<int> &index_2d,
    const std::vector<int> &index_2d_need_facemap ) const
{
  // Perform basic index boundary check
  SYS_T::print_fatal_if( index_3d >= num_phy_domain_3d || index_3d < 0,
      "Error: Gmsh_FileIO::write_vol_h5, surface index is wrong.\n");

  const unsigned int num_2d_face = index_2d.size();

  for(unsigned int ii=0; ii<num_2d_face; ++ii)
    SYS_T::print_fatal_if( index_2d[ii] >= num_phy_domain_2d || index_2d[ii] < 0,
        "Error: Gmsh_FileIO::write_vol_h5, face index is wrong. \n");

  // Open an HDF5 file
  std::string h5_file_name( "Gmsh_" );
  h5_file_name.append( phy_3d_name[index_3d] );
  h5_file_name.append(".h5");

  std::cout<<"=== Gmsh_FileIO::write_vol_h5 for "
    <<phy_3d_name[index_3d]<<" associated with ";

  for(unsigned int ii=0; ii<num_2d_face; ++ii)
    std::cout<<phy_2d_name[ index_2d[ii] ]<<'\t';
  std::cout<<std::endl;

  hid_t file_id = H5Fcreate(h5_file_name.c_str(),
      H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  // Write the 3D domain at the root
  const std::string slash("/");
  const int domain_3d_idx = phy_3d_index[index_3d];
  std::cout<<"-----> write "<<phy_3d_name[ index_3d ]<<'\n';

  const int num_3d_cell = phy_3d_nElem[index_3d];
  const int nLocBas_3d  = ele_nlocbas[ domain_3d_idx ];

  h5w->write_intScalar(file_id, "num_face", num_2d_face);
  h5w->write_intScalar(file_id, "num_node", num_node);
  h5w->write_intScalar(file_id, "num_cell", num_3d_cell);
  h5w->write_intScalar(file_id, "ele_type", ele_type[domain_3d_idx]);
  h5w->write_intScalar(file_id, "nLocBas", nLocBas_3d);
  h5w->write_intScalar(file_id, "index_3d", index_3d);
  h5w->write_intScalar(file_id, "phy_tag", phy_3d_index[index_3d]);
  h5w->write_intScalar(file_id, "start_index_3d", phy_3d_start_index[index_3d]);
  h5w->write_doubleVector(file_id, "node", node);
  h5w->write_intVector(file_id, "IEN", eIEN[domain_3d_idx]);

  // Write the 2D domain in groups
  for(unsigned int ii=0; ii<num_2d_face; ++ii)
  {
    std::cout<<"-----> write /"<<phy_2d_name[ index_2d[ii] ]<<'\n';

    const int domain_2d_idx = phy_2d_index[ index_2d[ii] ];
    const int nLocBas_2d = ele_nlocbas[ domain_2d_idx ];
    const int num_2d_cell = phy_2d_nElem[ index_2d[ii] ];

    std::vector<int> face_ien_global( eIEN[domain_2d_idx] );
    std::vector<int> bcpt( face_ien_global );
    VEC_T::sort_unique_resize(bcpt);
    const int bcnumpt = VEC_T::get_size( bcpt );
    std::cout<<"      num of bc pt = "<<bcnumpt<<'\n';

    // tript stores the coordinates of the boundary points
    std::vector<double> surpt(3*bcnumpt, 0.0);
    for( int jj=0; jj<bcnumpt; ++jj )
    {
      surpt[jj*3]   = node[bcpt[jj]*3];
      surpt[jj*3+1] = node[bcpt[jj]*3+1];
      surpt[jj*3+2] = node[bcpt[jj]*3+2];
    }

    // Generate a mapper that maps the bc node to 1; other node to 0
    bool * bcmap = new bool [num_node];
    for(int jj=0; jj<num_node; ++jj) bcmap[jj] = 0;
    for(int jj=0; jj<bcnumpt; ++jj) bcmap[bcpt[jj]] = 1;

    int nVertex_2d {0}, nVertex_3d {0};
    if(nLocBas_2d == 3 || nLocBas_2d == 6) // tri-tet
    {
      nVertex_2d = 3;
      nVertex_3d = 4;
    }
    else if(nLocBas_2d == 4 || nLocBas_2d == 9) // quad-hex
    {
      nVertex_2d = 4;
      nVertex_3d = 8;
    }
    else
      SYS_T::print_fatal("Error: Gmsh_FileIO::write_vol_h5, unknown element type.\n");

    // Generate a list of volume elements that have boundary over the face
    std::vector<int> gelem {};
    for( int ee=0; ee<num_3d_cell; ++ee )
    {
      int total = 0;
      for(int jj=0; jj<nVertex_3d; ++jj)
        total += bcmap[ eIEN[domain_3d_idx][nLocBas_3d*ee+jj] ];

      if(total >= nVertex_2d) gelem.push_back(ee);
    }
    delete [] bcmap; bcmap = nullptr;
    std::cout<<"      "<<gelem.size()<<" elems have face over the boundary.\n";

    // Generate the local face element IEN array
    std::vector<int> face_ien_local {};
    for(int ee=0; ee<num_2d_cell; ++ee)
    {
      for(int jj=0; jj<nLocBas_2d; ++jj)
        face_ien_local.push_back(VEC_T::get_pos(bcpt, face_ien_global[nLocBas_2d*ee+jj]));
    }
    std::cout<<"      edge IEN generated. \n";

    // Loacate the volumetric element that the face element belongs to 
    std::vector<int> face2elem(num_2d_cell, -1);
    if( VEC_T::is_invec( index_2d_need_facemap, index_2d[ii] ) )
    {
      for(int ff=0; ff<num_2d_cell; ++ff)
      {
        std::vector<int> sur_node (nVertex_2d, -1); // the vertex indices of a surface element
        for(int kk {0}; kk < nVertex_2d; ++kk)
          sur_node[kk] = face_ien_global[ nLocBas_2d * ff + kk ];

        bool gotit = false;
        int ee = -1;
        while( !gotit && ee < int(gelem.size()) - 1 )
        {
          ee += 1;
          const int vol_elem = gelem[ee];
          std::vector<int> vol_node (nVertex_3d, -1); // the vertex indices of a volume element
          for(int jj {0}; jj < nVertex_3d; ++jj)
            vol_node[jj] = eIEN[domain_3d_idx][ nLocBas_3d * vol_elem + jj ];
          
          bool got_all_node = true;
          for(int kk {0}; kk < nVertex_2d; ++kk)
            got_all_node = got_all_node && VEC_T::is_invec(vol_node, sur_node[kk]);

          gotit = got_all_node;
        }
        if(gotit)
          face2elem[ff] = gelem[ee] + phy_3d_start_index[index_3d];
        else
          face2elem[ff] = -1;
      }
      std::cout<<"      face2elem mapping generated.\n";
    }

    // Record info
    std::string name_2d_domain(slash);
    name_2d_domain.append( std::to_string( static_cast<int>(ii) ) );
    hid_t g_id = H5Gcreate( file_id, name_2d_domain.c_str(),
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    h5w->write_string(g_id, "name", phy_2d_name[ index_2d[ii] ] );
    h5w->write_intScalar(g_id, "phy_tag", index_2d[ii]);
    h5w->write_intScalar(g_id, "num_cell", num_2d_cell);
    h5w->write_intScalar(g_id, "num_node", bcnumpt );
    h5w->write_intScalar(g_id, "ele_type", ele_type[domain_2d_idx]);
    h5w->write_intScalar(g_id, "nLocBas", nLocBas_2d);
    h5w->write_intVector(g_id, "IEN_glo", face_ien_global);
    h5w->write_intVector(g_id, "IEN_loc", face_ien_local);
    h5w->write_intVector(g_id, "pt_idx", bcpt);
    h5w->write_intVector(g_id, "face2elem", face2elem);
    h5w->write_doubleVector(g_id, "pt_coor", surpt);

    H5Gclose(g_id);
  } // End-loop-over-2d-face

  // Close the HDF5 file
  delete h5w;
  H5Fclose( file_id ); 
}

void Gmsh_FileIO::update_FSI_nodal_ordering()
{
  // First generate the fluid and solid node lists
  std::vector<int> new2old = eIEN[ phy_3d_index[0] ];
  VEC_T::sort_unique_resize( new2old );

  std::vector<int> snode = eIEN[ phy_3d_index[1] ];
  VEC_T::sort_unique_resize(snode);

  // append the node unique in the solid domain after the fluid node to generate
  // a new2old mapping 
  for( const int ii : snode )
  {
    // if solid node is NOT in the fluid node, append it
    if( !VEC_T::is_invec(new2old, ii) )
      new2old.push_back( ii );
  }

  // Now clean the snode vector to save memory
  VEC_T::clean( snode );

  SYS_T::print_fatal_if( static_cast<int>( new2old.size() ) != num_node, "Error: Gmsh_FildIO::update_FSI_nodal_ordering the number of nodes in the first two sub-domain does match the num_node!\n" );

  // Now generate the old2new mapping
  std::vector<int> old2new( num_node, 0 );

  for(int ii=0; ii<num_node; ++ii)
    old2new[ new2old[ii] ] = ii;

  // Now clean the new2old mapper to save memory
  VEC_T::clean( new2old );

  // Now update the nodal x-y-z coordinates
  std::vector<double> temp( 3 * num_node, 0.0 ); // temporary xyz coordinates based on new indices

  for( int ii=0; ii<num_node; ++ii )
  {
    const int new_idx = old2new[ ii ]; // ii is the old index

    temp[3 * new_idx + 0] = node[3 * ii + 0];
    temp[3 * new_idx + 1] = node[3 * ii + 1];
    temp[3 * new_idx + 2] = node[3 * ii + 2];
  }

  node = temp; // now update the node vector

  // Now update all the IEN arrays by the new nodal index set
  for(int ii=0; ii<num_phy_domain; ++ii)
  {
    const int len = VEC_T::get_size(eIEN[ii]);
    for( int jj=0; jj<len; ++jj )
      eIEN[ii][jj] = old2new[ eIEN[ii][jj] ];
  }
}

void Gmsh_FileIO::update_quadratic_tet_IEN( const int &index_3d )
{
  SYS_T::print_fatal_if(index_3d < 0 || index_3d >= num_phy_domain_3d,
      "Error: input index_3d is out of range.\n");

  const int domain_index = phy_3d_index[ index_3d ];

  const int ne = phy_3d_nElem[ index_3d ];

  const int nlocbas = ele_nlocbas[ domain_index ];

  SYS_T::print_fatal_if(nlocbas != 10,
      "Error: Gmsh_FileIO updata_quadratic_tet_IEN only works for 10-node quadratic element. \n");

  std::cout<<"=== Gmsh_FileIO::update_quadratic_tet_IEN for "<<phy_3d_name[index_3d]<<std::endl;;

  // Now updating the eIEN array
  for(int ee=0; ee<ne; ++ee)
  {
    const int temp = eIEN[domain_index][10*ee+8];
    eIEN[domain_index][10*ee+8] = eIEN[domain_index][10*ee+9];
    eIEN[domain_index][10*ee+9] = temp;
  }
}

void Gmsh_FileIO::update_quadratic_hex_IEN( const int &index_3d )
{
  SYS_T::print_fatal_if(index_3d < 0 || index_3d >= num_phy_domain_3d,
      "Error: input index_3d is out of range.\n");

  const int domain_index = phy_3d_index[ index_3d ];

  const int ne = phy_3d_nElem[ index_3d ];

  const int nlocbas = ele_nlocbas[ domain_index ];

  SYS_T::print_fatal_if(nlocbas != 27 && nlocbas != 20, "Error: Gmsh_FileIO updata_quadratic_hex_IEN only works for 27-node or 20-node quadratic element. \n");

  // Now updating the eIEN array
  for(int ee=0; ee<ne; ++ee)
  {
    const int temp9  = eIEN[domain_index][nlocbas * ee + 9];
    const int temp10 = eIEN[domain_index][nlocbas * ee + 10];
    const int temp11 = eIEN[domain_index][nlocbas * ee + 11];
    const int temp12 = eIEN[domain_index][nlocbas * ee + 12];
    const int temp13 = eIEN[domain_index][nlocbas * ee + 13];
    const int temp14 = eIEN[domain_index][nlocbas * ee + 14];
    const int temp15 = eIEN[domain_index][nlocbas * ee + 15];
    const int temp16 = eIEN[domain_index][nlocbas * ee + 16];
    const int temp17 = eIEN[domain_index][nlocbas * ee + 17];
    const int temp18 = eIEN[domain_index][nlocbas * ee + 18];
    const int temp19 = eIEN[domain_index][nlocbas * ee + 19];
    eIEN[domain_index][nlocbas * ee + 9] = temp11;
    eIEN[domain_index][nlocbas * ee + 10] = temp13;
    eIEN[domain_index][nlocbas * ee + 11] = temp9;
    eIEN[domain_index][nlocbas * ee + 12] = temp16;
    eIEN[domain_index][nlocbas * ee + 13] = temp18;
    eIEN[domain_index][nlocbas * ee + 14] = temp19;
    eIEN[domain_index][nlocbas * ee + 15] = temp17;
    eIEN[domain_index][nlocbas * ee + 16] = temp10;
    eIEN[domain_index][nlocbas * ee + 17] = temp12;
    eIEN[domain_index][nlocbas * ee + 18] = temp14;
    eIEN[domain_index][nlocbas * ee + 19] = temp15;
    if (nlocbas == 27)
    {
      const int temp20 = eIEN[domain_index][27 * ee + 20];
      const int temp21 = eIEN[domain_index][27 * ee + 21];
      const int temp22 = eIEN[domain_index][27 * ee + 22];
      const int temp23 = eIEN[domain_index][27 * ee + 23];
      const int temp24 = eIEN[domain_index][27 * ee + 24];
      eIEN[domain_index][27 * ee + 20] = temp22;
      eIEN[domain_index][27 * ee + 21] = temp23;
      eIEN[domain_index][27 * ee + 22] = temp21;
      eIEN[domain_index][27 * ee + 23] = temp24;
      eIEN[domain_index][27 * ee + 24] = temp20;
    }
    else
    {
      // for the case of nlocbas = 20, we do not need to extra work for the 7
      // nodes that are associated with the tri-quadratic element.
    }
  }
}

void Gmsh_FileIO::write_quadratic_sur_vtu( const std::string &vtu_filename,
    const int &index_sur, const int &index_vol, const bool &isf2e, const bool &is_slave ) const
{
  SYS_T::print_fatal_if( index_sur >= num_phy_domain_2d || index_sur < 0,
      "Error: Gmsh_FileIO::write_quadratic_sur_vtu, surface index is wrong. \n");

  SYS_T::print_fatal_if( index_vol >= num_phy_domain_3d || index_vol < 0,
      "Error: Gmsh_FileIO::write_quadratic_sur_vtu, volume index is wrong. \n");

  std::cout<<"=== Gmsh_FileIO::write_quadratic_sur_vtu for "
    <<phy_2d_name[index_sur]
    <<" associated with "<<phy_3d_name[index_vol];

  if( isf2e )
    std::cout<<" with face-to-volume element index. \n";
  else
    std::cout<<" without face-to-volume element index. \n";

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  std::cout<<"-----> write "<<vtu_filename<<".vtu \n";

  mytimer->Reset(); mytimer->Start();

  const int phy_index_sur = phy_2d_index[index_sur];
  const int phy_index_vol = phy_3d_index[index_vol];

  // obtain the number of local basis function of the surface and volume domains
  const int nlocbas_2d {ele_nlocbas[phy_index_sur]};
  const int nlocbas_3d {ele_nlocbas[phy_index_vol]};

  // number of vertices of a element of the surface and volume 
  int nVertex_2d {-1};
  int nVertex_3d {-1};

  std::string ele_2d {};
  std::string ele_3d {};
  if (nlocbas_2d == 6 && nlocbas_3d == 10)
  { 
    nVertex_2d = 3;
    nVertex_3d = 4;
    ele_2d = static_cast<std::string>("triangle");
    ele_3d = static_cast<std::string>("tetrahedron");
  }
  else if (nlocbas_2d == 9 && nlocbas_3d == 27)
  {
    nVertex_2d = 4;
    nVertex_3d = 8;
    ele_2d = static_cast<std::string>("quadrilateral");
    ele_3d = static_cast<std::string>("hexahedron");
  }
  else
    SYS_T::print_fatal("Error: Gmsh_FileIO::write_quadratic_sur_vtu, element types of surface and volume donnot match. \n");

  // surface mesh ien copied from eIEN
  std::vector<int> sur_ien_global( eIEN[phy_index_sur] );

  // global node index
  std::vector<int> bcpt( sur_ien_global );
  const int bcnumcl = phy_2d_nElem[index_sur];

  SYS_T::print_fatal_if( VEC_T::get_size(sur_ien_global) != nlocbas_2d * bcnumcl,
      "Error: Gmsh_FileIO::write_quadratic_sur_vtu, sur IEN size wrong. \n" );

  VEC_T::sort_unique_resize( bcpt ); // unique ascending order nodes

  const int bcnumpt = VEC_T::get_size( bcpt );

  std::cout<<"      num of bc pt = "<<bcnumpt<<'\n';

  // tript stores the coordinates of the boundary points
  std::vector<double> sur_pt( 3*bcnumpt, 0.0 );

  PERIGEE_OMP_PARALLEL_FOR
  for( int ii=0; ii<bcnumpt; ++ii )
  {
    sur_pt[ii*3]   = node[bcpt[ii]*3] ;
    sur_pt[ii*3+1] = node[bcpt[ii]*3+1] ;
    sur_pt[ii*3+2] = node[bcpt[ii]*3+2] ;
  }

  // Volume mesh IEN
  std::vector<int> vol_IEN( eIEN[phy_index_vol] );
  const int numcel = phy_3d_nElem[index_vol];

  SYS_T::print_fatal_if( int( vol_IEN.size() ) != nlocbas_3d * numcel,
      "Error: Gmsh_FileIO::write_quadratic_sur_vtu, vol IEN size wrong. \n");

  // generate local triangle IEN array
  std::vector<int> sur_ien {};
  for(int ee=0; ee<bcnumcl; ++ee)
  { 
    for (int ii{0}; ii < nlocbas_2d; ++ii)
      sur_ien.push_back( VEC_T::get_pos(bcpt, sur_ien_global[nlocbas_2d * ee + ii]) );
  }
  std::cout<<"      " << ele_2d <<" IEN generated. \n";

  std::vector<int> face2elem( bcnumcl, -1 );
  if( isf2e )
  {
    // A mapper that maps bc node to 1 other to 0
    bool * bcmap = new bool [num_node];
    for(int ii=0; ii<num_node; ++ii) bcmap[ii] = 0;
    for(int ii=0; ii<bcnumpt; ++ii) bcmap[bcpt[ii]] = 1;

    std::vector<int> gelem {};

    PERIGEE_OMP_PARALLEL
    {
      std::vector<int> temp_gelem {};
      PERIGEE_OMP_FOR
      for( int ee=0; ee<numcel; ++ee )
      {
        int total = 0;
        for (int jj{0}; jj < nVertex_3d; ++jj)
          total += bcmap[ vol_IEN[nlocbas_3d * ee + jj] ];
        if(total >= nVertex_2d) 
          temp_gelem.push_back(ee);
      }
      PERIGEE_OMP_CRITICAL
      VEC_T::insert_end(gelem, temp_gelem);
    }

    delete [] bcmap; bcmap = nullptr;
    std::cout<<"      "<<gelem.size()<<" "<<ele_3d<<"s have faces over the surface. \n";


    PERIGEE_OMP_PARALLEL_FOR
    for(int ff=0; ff<bcnumcl; ++ff)
    {
      std::vector<int> snode( nVertex_2d, -1 );
      for (int ii{0}; ii < nVertex_2d; ++ii)
        snode[ii] = sur_ien_global[nlocbas_2d * ff + ii];
      
      bool got_sur_elem = false;
      int ee = -1;
      while( !got_sur_elem && ee < VEC_T::get_size(gelem) - 1 )
      {
        ee += 1;
        const int vol_elem = gelem[ee];

        std::vector<int> vnode( nVertex_3d, -1 );
        for (int jj{0}; jj < nVertex_3d; ++jj)
          vnode[jj] = vol_IEN[nlocbas_3d * vol_elem + jj];

        bool got_all_vertices = true;
        for (int ii{0}; ii < nVertex_2d; ++ii)
        {
          const bool got_each_vertex = VEC_T::is_invec(vnode, snode[ii]);
          got_all_vertices = got_all_vertices && got_each_vertex;
        }
        got_sur_elem = got_all_vertices;
      }

      // If the boundary surface element is not found,
      // we write -1 as the mapping value
      if(got_sur_elem)
        face2elem[ff] = gelem[ee] + phy_3d_start_index[index_vol];
      else
        face2elem[ff] = -1;
    }
    std::cout<<"      face2elem mapping generated. \n";
  }
  
  std::vector<DataVecStr<int>> input_vtk_data {};
  input_vtk_data.push_back({bcpt, "GlobalNodeID", AssociateObject::Node});
  input_vtk_data.push_back({face2elem, "GlobalElementID", AssociateObject::Cell});

  if (is_slave)
  {
    std::vector<int> master_id ( bcnumpt, -1 );
    for(int ii{0}; ii < bcnumpt; ++ii)
    {
      // The position of slave node in per_slave vector
      const int pos_slave = VEC_T::get_pos(per_slave, bcpt[ii]);
      SYS_T::print_fatal_if( pos_slave == -1,
        "Error: Gmsh_FileIO::write_write_quadratic_sur_vtu, node %d of boundary %s is not a slave node.\n", bcpt[ii], phy_2d_name[index_sur].c_str());
      
      master_id[ii] = per_master[pos_slave];
    }
    std::cout<<"      master-slave mapping generated. \n";
    input_vtk_data.push_back({master_id, "MasterNodeID", AssociateObject::Node});
  }

  if (nlocbas_2d == 6){
    TET_T::write_quadratic_triangle_grid( vtu_filename, bcnumpt, bcnumcl,
      sur_pt, sur_ien, input_vtk_data );
  }
  else if (nlocbas_2d == 9){
    HEX_T::write_quadratic_quad_grid( vtu_filename, bcnumpt, bcnumcl,
      sur_pt, sur_ien, input_vtk_data );
  }
  else
    SYS_T::print_fatal("Error: Gmsh_FileIO::write_quadratic_sur_vtu, undefined element type. \n");

  mytimer->Stop();
  std::cout<<"      Time taken "<<mytimer->get_sec()<<" sec. \n";
  delete mytimer;
}

void Gmsh_FileIO::write_quadratic_sur_vtu(const std::string &vtu_filename,
  const std::string &phy_name_sur, const std::string &phy_name_vol,
  const bool &isf2e, const bool &is_slave) const
{
  const int index_sur = VEC_T::get_pos(phy_2d_name, phy_name_sur);
  SYS_T::print_fatal_if(index_sur == -1,
    "Error: Gmsh_FileIO::write_quadratic_sur_vtu, wrong physical name of surface.\n");

  const int index_vol = VEC_T::get_pos(phy_3d_name, phy_name_vol);
  SYS_T::print_fatal_if(index_vol == -1,
    "Error: Gmsh_FileIO::write_quadratic_sur_vtu, wrong physical name of volume.\n");

  write_quadratic_sur_vtu(vtu_filename, index_sur, index_vol, isf2e, is_slave);
}

void Gmsh_FileIO::read_msh2(std::ifstream &infile)
{
  std::istringstream sstrm;
  std::string sline;

  getline(infile, sline); 
  SYS_T::print_fatal_if(sline.compare("$EndMeshFormat") != 0, 
      "Error: .msh format third line should be $EndMeshFormat. \n");
  
  getline(infile, sline);
  SYS_T::print_fatal_if(sline.compare("$PhysicalNames") != 0, 
      "Error: .msh format fourth line should be $PhysicalNames. \n");
  
  sstrm.clear();
  getline(infile, sline); sstrm.str(sline); sstrm>>num_phy_domain;

  // For each physical domain, read their index, we assume that
  // in the gmsh file, the index ranges from [1, num_phy_domain].
  // read the domain's dimension, should be 2 or 3;
  // read the domain's name, and remove the " " at the names.
  int pdim, pidx; std::string pname;
  for(int ii=0; ii<num_phy_domain; ++ii)
  {
    sstrm.clear(); getline(infile, sline); sstrm.str(sline);
    sstrm >> pdim; sstrm >> pidx; sstrm >> pname;

    // Gmsh have "name" as the physical name, first removes
    // the two primes in the name.
    pname.erase( pname.begin() ); 
    pname.erase( pname.end()-1 );
    
    // minus 1 because .msh file index starts from 1
    phy_dim.push_back(pdim);
    phy_index.push_back(pidx-1);
    phy_name.push_back(pname);
  }
  
  // Check the phy_index is within the rage
  // Make sure the physical domain index is in the range [1, num_phy_domain].
  std::vector<int> temp_phy_idx( phy_index );
  VEC_T::sort_unique_resize(temp_phy_idx);
  for(int ii=0; ii<num_phy_domain; ++ii)
  {
    SYS_T::print_fatal_if(temp_phy_idx[ii] != static_cast<int>(ii), 
      "Error: in the .msh file, the physical domain index should be in the rage [1, num_phy_domain]. \n");
  }

  // file syntax $EndPhysicalNames $Nodes 
  getline(infile, sline); 
  SYS_T::print_fatal_if(sline.compare("$EndPhysicalNames") != 0, 
      "Error: .msh format line should be $EndPhysicalNames. \n");
  
  getline(infile, sline);
  SYS_T::print_fatal_if(sline.compare("$Nodes") != 0, 
      "Error: .msh format line should be $Nodes. \n");
  
  // get the number of nodes
  sstrm.clear();
  getline(infile, sline); sstrm.str(sline); sstrm>>num_node;

  // Record x-y-z coordinates for the nodes into the vector node
  node.resize(3*num_node);

  for(int ii=0; ii<num_node; ++ii)
  {
    sstrm.clear(); getline(infile, sline); sstrm.str(sline);
    int nidx; // node index
    sstrm >> nidx;
  
    SYS_T::print_fatal_if( nidx != ii+1, "Error: .msh file, the nodal index should be in the range [1, num_node]. \n");

    sstrm >> node[ii*3];
    sstrm >> node[ii*3+1];
    sstrm >> node[ii*3+2];
  }

  // file syntax $EndNodes, $Elements
  getline(infile, sline); 
  SYS_T::print_fatal_if(sline.compare("$EndNodes") != 0, 
      "Error: .msh format line should be $EndNodes. \n");
  
  getline(infile, sline);
  SYS_T::print_fatal_if(sline.compare("$Elements") != 0, 
      "Error: .msh format line should be $Elements. \n");

  // get the number of elements
  sstrm.clear();
  getline(infile, sline); sstrm.str(sline); sstrm>>num_elem;

  // element IEN
  int eidx, etype, num_tag, phy_tag, geo_tag, enum_node = 0;

  phy_domain_nElem.assign(num_phy_domain, 0);

  // We assume that the physical tag ranges from 1 to num_phy_domain
  eIEN.resize(num_phy_domain);
  ele_nlocbas.resize(num_phy_domain);
  ele_type.resize(num_phy_domain);

  for(int ii=0; ii<num_phy_domain; ++ii)
  {
    eIEN[ii].clear();
    ele_nlocbas[ii] = -1;
  }

  for(int ii=0; ii<num_elem; ++ii)
  {
    sstrm.clear(); getline(infile, sline); sstrm.str(sline);
    sstrm >> eidx; sstrm >> etype; sstrm >> num_tag;
    
    // number of tag has to be 2, physical tag and geometrical tag,
    // based on the Gmsh default setting.
    SYS_T::print_fatal_if( num_tag!=2,
        "Error: .msh file number of tag for element is not 2.\n");

    // The pre-defined const array in the beginning of the constructor
    // gives the element number of nodes
    enum_node = elem_nlocbas[ etype ];

    // elem_phy_tag stores phy_tag-1 to be compatible with phy_index
    sstrm >> phy_tag; //elem_phy_tag.push_back(phy_tag -1);
    sstrm >> geo_tag; //elem_geo_tag.push_back(geo_tag);

    // Add the number of element to the physical domain
    phy_domain_nElem[phy_tag - 1] += 1;
    
    // Record the number of basis function (element type) in 
    // the physical domain. If there are different type element in the
    // physical subdomain, throw an error message.
    if( ele_nlocbas[phy_tag-1] == -1 )
    {
      ele_nlocbas[phy_tag-1] = enum_node;
      ele_type[phy_tag-1] = etype;
    }
    else 
    {
      SYS_T::print_fatal_if( ele_type[phy_tag-1] != etype,
        "Error: the physical domain have mixed type of elements. \n" );
    }
    
    // Record the IEN array for the element
    for(int jj=0; jj<enum_node; ++jj)
    {
      int temp_index;
      sstrm >> temp_index;
      // to make the IEN compatible with the c array: we correct
      // the node index by minus 1.
      eIEN[phy_tag-1].push_back( temp_index - 1 );
    }
  }
}

void Gmsh_FileIO::read_msh4(std::ifstream &infile)
{
  std::istringstream sstrm;
  std::string sline;

    getline(infile, sline); 
  SYS_T::print_fatal_if(sline.compare("$EndMeshFormat") != 0, 
      "Error: .msh format third line should be $EndMeshFormat. \n");
  
  getline(infile, sline);
  SYS_T::print_fatal_if(sline.compare("$PhysicalNames") != 0, 
      "Error: .msh format fourth line should be $PhysicalNames. \n");

  getline(infile, sline); sstrm.str(sline); sstrm>>num_phy_domain;

  // For each physical domain, read their index, we assume that
  // in the gmsh file, the index ranges from [1, num_phy_domain].
  // read the domain's dimension, should be 2 or 3;
  // read the domain's name, and remove the " " at the names.
  int pdim, pidx; std::string pname;
  for(int ii=0; ii<num_phy_domain; ++ii)
  {
    sstrm.clear(); getline(infile, sline); sstrm.str(sline);
    sstrm >> pdim; sstrm >> pidx; sstrm >> pname;

    // Gmsh have "name" as the physical name, first removes
    // the two primes in the name.
    pname.erase( pname.begin() ); 
    pname.erase( pname.end()-1 );
    
    // minus 1 because .msh file index starts from 1
    phy_dim.push_back(pdim);
    phy_index.push_back(pidx-1);
    phy_name.push_back(pname);
  }

  // Check the phy_index is within the rage
  // Make sure the physical domain index is in the range [1, num_phy_domain].
  // Note: here our phy_index ranges in [ 0, num_phy_domain-1], as we have made
  // a correction in above to make it start from zero.
  std::vector<int> temp_phy_idx( phy_index );
  VEC_T::sort_unique_resize(temp_phy_idx);
  for(int ii=0; ii<num_phy_domain; ++ii)
  {
    SYS_T::print_fatal_if(temp_phy_idx[ii] != static_cast<int>(ii), 
      "Error: in the .msh file, the physical domain index should be in the rage [1, num_phy_domain]. \n");
  }

  // file syntax $EndPhysicalNames $Nodes 
  getline(infile, sline); 
  SYS_T::print_fatal_if(sline.compare("$EndPhysicalNames") != 0, 
      "Error: .msh format line should be $EndPhysicalNames. \n");
  
  getline(infile, sline);
  SYS_T::print_fatal_if(sline.compare("$Entities") != 0, 
      "Error: .msh format line should be $Entities. \n");

  // get the number of original points, curves, surfaces and volumes
  int numPoints, numCurves, numSurfaces, numVolumes;
  sstrm.clear(); getline(infile, sline); sstrm.str(sline);
  sstrm >> numPoints; sstrm >> numCurves; sstrm >> numSurfaces; sstrm >> numVolumes;

  // skip over the information of original points
  for (int point{0}; point < numPoints; ++point)
    getline(infile, sline);
  // If we assign physical groups to these points, we cannot skip the information above

  // build up the relationship of entities(GeoTag) and physical groups(PhyTag)
  int GeoTag, PhyTag, num_phy_tag; 
  double minX, minY, minZ, maxX, maxY, maxZ;

  // read each line for the information of original curves
  std::vector<int> curveTag (numCurves, -1);
  std::vector<int> curvePhyTag (numCurves, -1);
  for (int curve{0}; curve < numCurves; ++curve)
  {
    sstrm.clear(); getline(infile, sline); sstrm.str(sline);

    sstrm >> GeoTag;
    curveTag[curve] = GeoTag;

    sstrm >> minX; sstrm >> minY; sstrm >> minZ;
    sstrm >> maxX; sstrm >> maxY; sstrm >> maxZ;

    sstrm >> num_phy_tag;
    if (num_phy_tag == 0)
    {
      // some curves are not in any physical group, set phy_tag = -1 by default
      continue;
    }
    else
    {
      // here we suppose one entity belongs to at most one physical group
      SYS_T::print_fatal_if( num_phy_tag != 1,
        "Error: .msh file number of physical tag for a curve is not 1.\n");
      sstrm >> PhyTag;
      curvePhyTag[curve] = PhyTag;
    }
  }

  // read each line for the information of original surfaces
  std::vector<int> surfaceTag (numSurfaces, -1);
  std::vector<int> surfacePhyTag (numSurfaces, -1);
  for (int surface{0}; surface < numSurfaces; ++surface)
  {
    sstrm.clear(); getline(infile, sline); sstrm.str(sline);

    sstrm >> GeoTag;
    surfaceTag[surface] = GeoTag;

    sstrm >> minX; sstrm >> minY; sstrm >> minZ;
    sstrm >> maxX; sstrm >> maxY; sstrm >> maxZ;

    sstrm >> num_phy_tag;
    if (num_phy_tag == 0)
    {
      // some surfaces are not in any physical group, set phy_tag = -1 by default
      continue;
    }
    else
    {
      // here we suppose one entity belongs to at most one physical group
      SYS_T::print_fatal_if( num_phy_tag != 1,
        "Error: .msh file number of physical tag for a surface is not 1.\n");
      sstrm >> PhyTag;
      surfacePhyTag[surface] = PhyTag;
    }
  }

  // read each line for the information of original volumes
  std::vector<int> volumeTag (numVolumes, -1);
  std::vector<int> volumePhyTag (numVolumes, -1);
  for (int volume{0}; volume < numVolumes; ++ volume)
  {
    sstrm.clear(); getline(infile, sline); sstrm.str(sline);

    sstrm >> GeoTag;
    volumeTag[volume] = GeoTag;

    sstrm >> minX; sstrm >> minY; sstrm >> minZ;
    sstrm >> maxX; sstrm >> maxY; sstrm >> maxZ;

    sstrm >> num_phy_tag;
    if (num_phy_tag == 0)
    {
      // some volumes are not in any physical group, set phy_tag = -1 by default
      continue;
    }
    else
    {
      // here we suppose one entity belongs to at most one physical group
      SYS_T::print_fatal_if( num_phy_tag != 1,
        "Error: .msh file number of physical tag for a volume is not 1.\n");
      sstrm >> PhyTag;
      volumePhyTag[volume] = PhyTag;
    }
  }

  getline(infile, sline);
  SYS_T::print_fatal_if(sline.compare("$EndEntities") != 0, 
    "Error: .msh format line should be $EndEntities. \n");

  // here we suppose there is no partitioned entities, the next line should be '$Nodes'
  getline(infile, sline);
  SYS_T::print_fatal_if(sline.compare("$Nodes") != 0, 
     "Error: .msh format line should be $Nodes. \n");
  
  // get the number of nodes
  int num_blocks;    // numEntityBlocks
  sstrm.clear(); getline(infile, sline); sstrm.str(sline);
  sstrm >> num_blocks; sstrm>>num_node;

  // Record x-y-z coordinates for the nodes into the vector node
  node.resize(3*num_node);

  int recorded_node_num {0};
  for (int block{0}; block < num_blocks; ++block)
  {
    sstrm.clear(); getline(infile, sline); sstrm.str(sline);
    int entity_dim, entity_tag, parametric, num_node_in_block;
    sstrm >> entity_dim; sstrm >> entity_tag; sstrm >> parametric; sstrm >> num_node_in_block;

    // here we suppose parametric = 0, then there is no <u>, <v>, <w>
    SYS_T::print_fatal_if( parametric != 0, 
      "Error: .msh file, the third parameter of nodal blocks should be 0 in block %d.\n", block);

    // here we suppose the nodal index should be in [1, num_node]
    for (int ii{0}; ii < num_node_in_block; ++ii)
    {
      sstrm.clear(); getline(infile, sline); sstrm.str(sline);
      int nidx;
      sstrm >> nidx;
      SYS_T::print_fatal_if( nidx != recorded_node_num + ii + 1, 
        "Error: .msh file, the nodal index should be in the range [1, num_node]. \n");
    }

    // record coordinates
    for (int ii{0}; ii < num_node_in_block; ++ii)
    {
      sstrm.clear(); getline(infile, sline); sstrm.str(sline);
      int coor_idx{recorded_node_num + ii};
      sstrm >> node[coor_idx * 3];
      sstrm >> node[coor_idx * 3 + 1];
      sstrm >> node[coor_idx * 3 + 2];
    }
    
    // finish record in this block
    recorded_node_num += num_node_in_block;
  }
  
  // check
  SYS_T::print_fatal_if( recorded_node_num != num_node, 
    "Error: .msh file, the number of recorded nodes does not match with the number of nodes . \n");

  // file syntax $EndNodes, $Elements
  getline(infile, sline); 
  SYS_T::print_fatal_if(sline.compare("$EndNodes") != 0, 
      "Error: .msh format line should be $EndNodes. \n");
  
  getline(infile, sline);
  SYS_T::print_fatal_if(sline.compare("$Elements") != 0, 
      "Error: .msh format line should be $Elements. \n");

  // get the number of elements
  sstrm.clear(); getline(infile, sline); sstrm.str(sline);
  sstrm >> num_blocks; sstrm>>num_elem;

  // We assume that the physical tag ranges from 1 to num_phy_domain
  eIEN.resize(num_phy_domain);
  phy_domain_nElem.resize(num_phy_domain);
  ele_nlocbas.resize(num_phy_domain);
  ele_type.resize(num_phy_domain);

  for(int ii=0; ii<num_phy_domain; ++ii)
  {
    eIEN[ii].clear();
    phy_domain_nElem[ii] = 0;
    ele_nlocbas[ii] = -1;
  }

  int recorded_ele_num {0};
  for (int block{0}; block < num_blocks; ++block)
  {
    sstrm.clear(); getline(infile, sline); sstrm.str(sline);
    int entity_dim, geo_tag, etype, num_ele_in_block;
    sstrm >> entity_dim; sstrm >> geo_tag; sstrm >> etype; sstrm >> num_ele_in_block;

    // The pre-defined const array in the beginning of the constructor
    // gives the element number of nodes
    int enum_node {elem_nlocbas [etype]};

    // find physical tag according to the entity information
    int phy_tag {-1}, geo_index {-1};
    if (entity_dim == 1)  // line element
    {
      geo_index = VEC_T::get_pos(curveTag, geo_tag);
      phy_tag = curvePhyTag[geo_index];
    }
    else if (entity_dim == 2) // surface element
    {
      geo_index = VEC_T::get_pos(surfaceTag, geo_tag);
      phy_tag = surfacePhyTag[geo_index];
    }
    else if (entity_dim == 3) // volume element
    {
      geo_index = VEC_T::get_pos(volumeTag, geo_tag);
      phy_tag = volumePhyTag[geo_index];
    }
    else
      SYS_T::print_fatal( "Error: .msh file, the dimension of element should be 1, 2 or 3.\n" );

    // Record the number of basis function (element type) in 
    // the physical domain. If there are different type element in the
    // physical subdomain, throw an error message.
    if( ele_nlocbas[phy_tag-1] == -1 )
    {
      ele_nlocbas[phy_tag-1] = enum_node;
      ele_type[phy_tag-1] = etype;
    }
    else 
    {
      SYS_T::print_fatal_if( ele_type[phy_tag-1] != etype,
        "Error: the physical domain have mixed type of elements. \n" );
    }

    // read each line of element information in a block
    for (int ii{0}; ii < num_ele_in_block; ++ii)
    {
      sstrm.clear(); getline(infile, sline); sstrm.str(sline);
      int eidx;
      sstrm >> eidx;
      // Record the IEN array for the element
      for(int jj=0; jj<enum_node; ++jj)
      {
        int temp_index;
        sstrm >> temp_index;
        // to make the IEN compatible with the c array: we correct
        // the node index by minus 1.
        eIEN[phy_tag-1].push_back( temp_index - 1 );
      } 
    }
    // Add the number of element to the physical domain
    phy_domain_nElem[phy_tag - 1] += num_ele_in_block;
    recorded_ele_num += num_ele_in_block;
  }

  // check
  SYS_T::print_fatal_if( recorded_ele_num != num_elem, 
    "Error: .msh file, the number of recorded elements does not match with the number of elements. \n");
}

void Gmsh_FileIO::read_periodic(std::ifstream &infile)
{
  std::istringstream sstrm;
  std::string sline;

  getline(infile, sline); sstrm.str(sline);
  int numPeriodicLinks;
  sstrm >> numPeriodicLinks;

  for(int link{0}; link < numPeriodicLinks; ++link)
  {
    getline(infile, sline); // skip: entityDim, entityTag, entityTagMaster
    getline(infile, sline); // skip: numAffine and affine

    sstrm.clear(); getline(infile, sline); sstrm.str(sline);
    int numCorrespondingNodes;
    sstrm >> numCorrespondingNodes;

    for(int node_pair{0}; node_pair < numCorrespondingNodes; ++node_pair)
    {
      sstrm.clear(); getline(infile, sline); sstrm.str(sline);
      int nodeTag, nodeTagMaster;
      sstrm >> nodeTag; sstrm >> nodeTagMaster;

      if(VEC_T::is_invec(per_slave, nodeTag - 1))
        ; // When we use periodic BC, if a slave node have more than one master, they will follow 
        // a common primary master. Hence there is no need to assign more than one master to a slave.
      else
      {
        per_slave.push_back(nodeTag - 1);
        per_master.push_back(nodeTagMaster - 1);
        // Node index of .msh file starts from 1,
        // but in ien array our node index starts from 0.
      }
    }
  }

  // file syntax $EndPeriodic
  getline(infile, sline);
  SYS_T::print_fatal_if( sline.compare("$EndPeriodic") != 0, 
     "Error: .msh format line should be $EndPeriodic. \n");

  // find primary masters
  for(int &master : per_master)
  {
    int pos_as_slave = VEC_T::get_pos(per_slave, master);
    while(pos_as_slave != -1) // the master is others' slave
    {
      master = per_master[pos_as_slave];
      pos_as_slave = VEC_T::get_pos(per_slave, master);
    }
  }
}

// EOF
