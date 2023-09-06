#include "Gmsh_FileIO.hpp"

Gmsh_FileIO::Gmsh_FileIO( const std::string &in_file_name )
: filename( in_file_name )
{
  // This is the element-type-to-num-of-local-basis mapping
  // based on the Gmsh format. The first is zero because Gmsh
  // start the element type number with 1. Detailed definition
  // of the element type is in elm-type of the MSH ASCII file
  // format.
  // The elm type 1 is a two-node line
  // The elm type 2 is a three-node triangle
  // The elm type 3 is a 4-node quadrangle
  // The elm type 4 is a 4-node tetrahedron
  // ...
  // The elm type 31 is a 56-node fifth-order tetrahedron
  const int elem_nlocbas[] = {0, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9,
    10, 27, 18, 14, 1, 8, 20, 15, 13, 9, 10, 12, 15, 15, 21, 
    4, 5, 6, 20, 35, 56};

  // Setup the file instream
  std::ifstream infile( filename.c_str(), std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;

  // First four lines are Gmsh default file format 
  getline(infile, sline); 
  SYS_T::print_fatal_if(sline.compare("$MeshFormat") != 0, 
      "Error: .msh format first line should be $MeshFormat. \n");
  
  getline(infile, sline);
  SYS_T::print_fatal_if(sline.compare("2.2 0 8") != 0, 
      "Error: .msh format second line should be 2.2 0 8. \n");
  
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
  std::vector<std::string> phy_name;
  std::vector<int> phy_index, phy_dim;
  phy_index.clear(); phy_dim.clear(); phy_name.clear();
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
    if( temp_phy_idx[ii] != static_cast<int>(ii) )
      SYS_T::print_fatal("Error: in the .msh file, the physical domain index should be in the rage [1, num_phy_domain]. \n");
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

  // stores the ii-th domain's number of elements. ii is the physical tag
  std::vector<int> phy_domain_nElem;

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
    if( ele_nlocbas[phy_tag-1] == -1 ){
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

  // Finish the file reading, the last line should be $EndElements  
  getline(infile, sline);
  SYS_T::print_fatal_if(sline.compare("$EndElements") != 0, 
      "Error: .msh format line should be $EndElements. \n");

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

void Gmsh_FileIO::write_interior_vtp( const int &index_sur,
    const int &index_vol1, const int &index_vol2 ) const
{
  SYS_T::print_fatal_if( index_sur >= num_phy_domain_2d || index_sur < 0,
      "Error: Gmsh_FileIO::write_vtp, surface index is wrong. \n");

  SYS_T::print_fatal_if( index_vol1 >= num_phy_domain_3d || index_vol1 < 0,
      "Error: Gmsh_FileIO::write_vtp, volume 1 index is wrong. \n");

  SYS_T::print_fatal_if( index_vol2 >= num_phy_domain_3d || index_vol2 < 0,
      "Error: Gmsh_FileIO::write_vtp, volume 2 index is wrong. \n");

  const int phy_index_sur = phy_2d_index[index_sur];
  const int phy_index_vol_1 = phy_3d_index[index_vol1];
  const int phy_index_vol_2 = phy_3d_index[index_vol2];
  const int bcnumcl = phy_2d_nElem[index_sur];
  const int numcel_1 = phy_3d_nElem[index_vol1];
  const int numcel_2 = phy_3d_nElem[index_vol2];

  std::cout<<"=== Gmsh_FileIO::write_interior_vtp for "
    <<phy_2d_name[index_sur]
    <<" associated with "<<phy_3d_name[index_vol1]
    <<" and "<<phy_3d_name[index_vol2]<<std::endl;

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  std::string vtp_file_name(phy_2d_name[index_sur]);
  vtp_file_name += "_";
  vtp_file_name += phy_3d_name[index_vol1];
  vtp_file_name += "_";
  vtp_file_name += phy_3d_name[index_vol2];
  std::cout<<"-----> write "<<vtp_file_name<<".vtp \n";
  mytimer->Reset();
  mytimer->Start();

  // Obtain the surface's IEN array associated with global nodal indices
  std::vector<int> trien_global( eIEN[phy_index_sur] );

  SYS_T::print_fatal_if( int(trien_global.size() ) != 3 * bcnumcl,
      "Error: Gmsh_FileIO::write_vtp, sur IEN size wrong. \n" );

  // bcpt stores the global nodal indices in ascending order
  std::vector<int> bcpt( trien_global );
  VEC_T::sort_unique_resize( bcpt ); // unique ascending order nodes
  const int bcnumpt = static_cast<int>( bcpt.size() );
  std::cout<<"      num of bc pt = "<<bcnumpt<<'\n';

  // Obtain the volumetric IEN array
  std::vector<int> vol_IEN_1( eIEN[phy_index_vol_1] );
  std::vector<int> vol_IEN_2( eIEN[phy_index_vol_2] );

  // tript stores the xyz coordinates of the points for the surface
  std::vector<double> tript; tript.clear(); tript.resize(3*bcnumpt);
  for( int ii=0; ii<bcnumpt; ++ii )
  {
    tript[ii*3]   = node[bcpt[ii]*3] ;
    tript[ii*3+1] = node[bcpt[ii]*3+1] ;
    tript[ii*3+2] = node[bcpt[ii]*3+2] ;
  }

  // generate a mapper for the global nodes that returns 1 if the nodes 
  // belong to the surface, 0 otherwise.
  bool * bcmap = new bool [num_node];
  for(int ii=0; ii<num_node; ++ii) bcmap[ii] = 0;
  for(int ii=0; ii<bcnumpt; ++ii) bcmap[bcpt[ii]] = 1;

  // use gelem_1[2] to store the vol elements that have face over the surface
  // mesh 
  std::vector<int> gelem_1; gelem_1.clear();
  for( int ee=0; ee<numcel_1; ++ee )
  {
    int total = 0;
    total += bcmap[ vol_IEN_1[4*ee] ];
    total += bcmap[ vol_IEN_1[4*ee+1] ];
    total += bcmap[ vol_IEN_1[4*ee+2] ];
    total += bcmap[ vol_IEN_1[4*ee+3] ];
    if(total >= 3) gelem_1.push_back(ee);
  }
  std::cout<<"      vol 1 domain: "<<gelem_1.size()<<" tets have more than 3 points on the surface. \n";

  std::vector<int> gelem_2; gelem_2.clear();
  for( int ee=0; ee<numcel_2; ++ee )
  {
    int total = 0;
    total += bcmap[ vol_IEN_2[4*ee] ];
    total += bcmap[ vol_IEN_2[4*ee+1] ];
    total += bcmap[ vol_IEN_2[4*ee+2] ];
    total += bcmap[ vol_IEN_2[4*ee+3] ];
    if(total >= 3) gelem_2.push_back(ee);
  }
  std::cout<<"      vol 2 domain: "<<gelem_2.size()<<" tets have more than 3 points on the surface. \n";

  // generate the local triangle IEN array
  std::vector<int> trien; trien.clear();
  int node0, node1, node2;
  std::vector<int>::iterator it;
  for(int ee=0; ee<bcnumcl; ++ee)
  {
    node0 = trien_global[3*ee];
    node1 = trien_global[3*ee+1];
    node2 = trien_global[3*ee+2];

    it = find(bcpt.begin(), bcpt.end(), node0);
    trien.push_back( it - bcpt.begin() );

    it = find(bcpt.begin(), bcpt.end(), node1);
    trien.push_back( it - bcpt.begin() );

    it = find(bcpt.begin(), bcpt.end(), node2);
    trien.push_back( it - bcpt.begin() );
  }
  std::cout<<"      triangle IEN generated. \n";

  // determine the face-2-element mapping
  std::vector<int> face2elem_1; face2elem_1.resize( bcnumcl );
  int vol_elem;
  int vnode[4];
  bool got0, got1, got2, gotit;
  for(int ff=0; ff<bcnumcl; ++ff)
  {
    node0 = trien_global[3*ff];
    node1 = trien_global[3*ff+1];
    node2 = trien_global[3*ff+2];
    gotit = false;
    int ee = -1;
    while( !gotit && ee < int(gelem_1.size()) - 1 )
    {
      ee += 1;
      vol_elem = gelem_1[ee];
      vnode[0] = vol_IEN_1[4*vol_elem];
      vnode[1] = vol_IEN_1[4*vol_elem+1];
      vnode[2] = vol_IEN_1[4*vol_elem+2];
      vnode[3] = vol_IEN_1[4*vol_elem+3];
      std::sort(vnode, vnode+4);

      got0 = ( std::find(vnode, vnode+4, node0) != vnode+4 );
      got1 = ( std::find(vnode, vnode+4, node1) != vnode+4 );
      got2 = ( std::find(vnode, vnode+4, node2) != vnode+4 );
      gotit = got0 && got1 && got2;
    }

    if(gotit)
      face2elem_1[ff] = gelem_1[ee] + phy_3d_start_index[index_vol1];
    else
    {
      face2elem_1[ff] = -1;
      std::cout<<"Warning: there are surface element not found in the volumetric mesh.\n";
    }
  }

  std::vector<int> face2elem_2; face2elem_2.resize( bcnumcl );
  for(int ff=0; ff<bcnumcl; ++ff)
  {
    node0 = trien_global[3*ff];
    node1 = trien_global[3*ff+1];
    node2 = trien_global[3*ff+2];
    gotit = false;
    int ee = -1;
    while( !gotit && ee < int(gelem_2.size()) - 1 )
    {
      ee += 1;
      vol_elem = gelem_2[ee];
      vnode[0] = vol_IEN_2[4*vol_elem];
      vnode[1] = vol_IEN_2[4*vol_elem+1];
      vnode[2] = vol_IEN_2[4*vol_elem+2];
      vnode[3] = vol_IEN_2[4*vol_elem+3];
      std::sort(vnode, vnode+4);

      got0 = ( std::find(vnode, vnode+4, node0) != vnode+4 );
      got1 = ( std::find(vnode, vnode+4, node1) != vnode+4 );
      got2 = ( std::find(vnode, vnode+4, node2) != vnode+4 );
      gotit = got0 && got1 && got2;
    }

    if(gotit)
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
  TET_T::write_triangle_grid( vtp_file_name, bcnumpt, bcnumcl,
      tript, trien, input_vtk_data );

  delete [] bcmap; bcmap = nullptr;
  mytimer->Stop();
  std::cout<<"      Time taken "<<mytimer->get_sec()<<" sec. \n";
  delete mytimer;
}

void Gmsh_FileIO::write_vtp(const int &index_sur, 
    const int &index_vol, const bool &isf2e ) const
{
  SYS_T::print_fatal_if( index_sur >= num_phy_domain_2d || index_sur < 0,
      "Error: Gmsh_FileIO::write_vtp, surface index is wrong. \n");

  SYS_T::print_fatal_if( index_vol >= num_phy_domain_3d || index_vol < 0,
      "Error: Gmsh_FileIO::write_vtp, volume index is wrong. \n");

  const int phy_index_sur = phy_2d_index[index_sur];
  const int phy_index_vol = phy_3d_index[index_vol];
  const int bcnumcl = phy_2d_nElem[index_sur];
  const int numcel = phy_3d_nElem[index_vol];

  std::cout<<"=== Gmsh_FileIO::write_vtp for "
    <<phy_2d_name[index_sur]
    <<" associated with "<<phy_3d_name[index_vol];

  if( isf2e )
    std::cout<<" with face-to-volume element index. \n";
  else
    std::cout<<" without face-to-volume element index. \n";

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  std::string vtp_file_name(phy_2d_name[index_sur]);
  vtp_file_name += "_";
  vtp_file_name += phy_3d_name[index_vol];
  std::cout<<"-----> write "<<vtp_file_name<<".vtp \n";
  mytimer->Reset();
  mytimer->Start();

  // Copy the IEN from the whole domain, the nodal indices is from the
  // global domain indices.
  std::vector<int> trien_global( eIEN[phy_index_sur] );

  // bcpt stores the global node index
  std::vector<int> bcpt( trien_global );

  SYS_T::print_fatal_if( int(trien_global.size() ) != 3 * bcnumcl,
      "Error: Gmsh_FileIO::write_vtp, sur IEN size wrong. \n" );

  // obtain the volumetric mesh IEN array
  std::vector<int> vol_IEN( eIEN[phy_index_vol] );

  SYS_T::print_fatal_if( int( vol_IEN.size() ) != 4 * numcel,
      "Error: Gmsh_FileIO::write_vtp, vol IEN size wrong. \n");

  VEC_T::sort_unique_resize( bcpt ); // unique ascending order nodes

  const int bcnumpt = static_cast<int>( bcpt.size() );

  std::cout<<"      num of bc pt = "<<bcnumpt<<'\n';

  // tript stores the coordinates of the boundary points
  std::vector<double> tript; tript.clear(); tript.resize(3*bcnumpt);
  for( int ii=0; ii<bcnumpt; ++ii )
  {
    tript[ii*3]   = node[bcpt[ii]*3] ;
    tript[ii*3+1] = node[bcpt[ii]*3+1] ;
    tript[ii*3+2] = node[bcpt[ii]*3+2] ;
  }

  // generate a mapper that maps the bc node to 1, other node to 0
  bool * bcmap = new bool [num_node];
  for(int ii=0; ii<num_node; ++ii) bcmap[ii] = 0;
  for(int ii=0; ii<bcnumpt; ++ii) bcmap[bcpt[ii]] = 1;

  std::vector<int> gelem; gelem.clear();
  for( int ee=0; ee<numcel; ++ee )
  {
    int total = 0;
    total += bcmap[ vol_IEN[4*ee] ];
    total += bcmap[ vol_IEN[4*ee+1] ];
    total += bcmap[ vol_IEN[4*ee+2] ];
    total += bcmap[ vol_IEN[4*ee+3] ];
    if(total >= 3) gelem.push_back(ee);
  }
  delete [] bcmap; bcmap = nullptr;
  std::cout<<"      "<<gelem.size()<<" tets have faces over the surface. \n";

  // generate the local triangle IEN array
  std::vector<int> trien; trien.clear();
  for(int ee=0; ee<bcnumcl; ++ee)
  {
    trien.push_back( VEC_T::get_pos(bcpt, trien_global[3*ee  ]) );
    trien.push_back( VEC_T::get_pos(bcpt, trien_global[3*ee+1]) );
    trien.push_back( VEC_T::get_pos(bcpt, trien_global[3*ee+2]) );
  }
  std::cout<<"      triangle IEN generated. \n";

  // determine the face-to-element mapping, if we demand it (
  // meaning this face needs boundary integral); otherwise, set
  // the face2elem as -1, since we will only need the nodal indices
  // for Dirichlet type face.
  std::vector<int> face2elem; face2elem.resize( bcnumcl, -1 );
  if( isf2e )
  {
    for(int ff=0; ff<bcnumcl; ++ff)
    {
      const int node0 = trien_global[3*ff];
      const int node1 = trien_global[3*ff+1];
      const int node2 = trien_global[3*ff+2];
      bool gotit = false;
      int ee = -1;
      while( !gotit && ee < int(gelem.size()) - 1 )
      {
        ee += 1;
        const int vol_elem = gelem[ee];
        int vnode[4] { vol_IEN[4*vol_elem], vol_IEN[4*vol_elem+1],
         vol_IEN[4*vol_elem+2], vol_IEN[4*vol_elem+3] };
        std::sort(vnode, vnode+4);

        const bool got0 = ( std::find(vnode, vnode+4, node0) != vnode+4 );
        const bool got1 = ( std::find(vnode, vnode+4, node1) != vnode+4 );
        const bool got2 = ( std::find(vnode, vnode+4, node2) != vnode+4 );
        gotit = got0 && got1 && got2;
      }

      // If the boundary surface element is not found, 
      // we write -1 as the mapping value
      if(gotit)
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
  TET_T::write_triangle_grid( vtp_file_name, bcnumpt, bcnumcl,
      tript, trien, input_vtk_data );

  mytimer->Stop();
  std::cout<<"      Time taken "<<mytimer->get_sec()<<" sec. \n";
  delete mytimer;
}

void Gmsh_FileIO::write_each_vtu() const
{
  std::cout<<"=== Gmsh_FileIO::wirte_each_vtu.\n";
  std::cout<<"--- There are "<<num_phy_domain_3d<<" 3D physical domains.\n";

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  for(int ii=0; ii<num_phy_domain_3d; ++ii)
  {
    const int domain_index = phy_3d_index[ii];
    const std::string vtu_file_name = phy_3d_name[ ii ];
    std::cout<<"-----> write "<<vtu_file_name<<".vtu \t";

    mytimer->Reset();
    mytimer->Start();

    std::cout<<" nElem = "<<phy_3d_nElem[ii]<<'\t';

    const int start_eindex = phy_3d_start_index[ii];

    std::cout<<" starting e index = "<<start_eindex<<'\t';

    // generate physics tag
    std::vector<int> ptag; ptag.clear();
    ptag.assign(phy_3d_nElem[ii], ii);
    
    // collect the FSI indices of the sub-domain nodes
    std::vector<int> local_node_idx = eIEN[ domain_index ];
    VEC_T::sort_unique_resize( local_node_idx );

    const int num_local_node = static_cast<int>( local_node_idx.size() );

    // collect those points' coordinates
    std::vector<double> local_coor; 
    local_coor.resize( 3 * num_local_node );
    for(int jj=0; jj<num_local_node; ++jj )
    {
      local_coor[ 3*jj+0 ] = node[ 3*local_node_idx[jj] + 0 ];
      local_coor[ 3*jj+1 ] = node[ 3*local_node_idx[jj] + 1 ];
      local_coor[ 3*jj+2 ] = node[ 3*local_node_idx[jj] + 2 ];
    }
    
    // generate a local IEN
    std::vector<int> domain_IEN;
    const int nloc = ele_nlocbas[domain_index];
    domain_IEN.resize( phy_3d_nElem[ii] * nloc );
    for(int ee=0; ee<phy_3d_nElem[ii]; ++ee)
    {
      for(int jj=0; jj<nloc; ++jj)
      {
        const int target = eIEN[ domain_index ][ ee*nloc + jj ];
        domain_IEN[ ee * nloc + jj ] = VEC_T::get_pos( local_node_idx, target );
      }
    } 

    // Element index (using the start_eindex)
    std::vector<int> local_cell_idx; local_cell_idx.resize( phy_3d_nElem[ii] );
    for(int ee=0; ee<phy_3d_nElem[ii]; ++ee)
      local_cell_idx[ee] = start_eindex + ee;

    // write the sub-volumetric domain's vtk/vtu file
    // the subdomain element index start with the start_eindex
    std::vector<DataVecStr<int>> input_vtk_data {};
    input_vtk_data.push_back({local_node_idx, "GlobalNodeID", AssociateObject::Node});
    input_vtk_data.push_back({local_cell_idx, "GlobalElementID", AssociateObject::Cell});
    input_vtk_data.push_back({ptag, "Physics_tag", AssociateObject::Cell});
    TET_T::write_tet_grid( vtu_file_name, num_local_node, 
        phy_3d_nElem[ ii ], local_coor, domain_IEN, input_vtk_data, true);

    mytimer->Stop();

    std::cout<<mytimer->get_sec()<<" sec. \n";
  }

  delete mytimer;
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
  std::vector<int> wIEN; // whole mesh IEN
  wIEN.clear();
  std::vector<int> wtag; // whole vol mehs phys tag
  wtag.clear();
  int wnElem = 0; // whole mesh number of elements

  // whole mesh num of node is assumed to be num_node 
  const int wnNode = num_node;

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
  
  
  TET_T::write_tet_grid( in_fname, wnNode, wnElem, node,
      wIEN, input_vtk_data, isXML ); 

  mytimer->Stop();
  std::cout<<"    Time taken "<<mytimer->get_sec()<<" sec. \n";
  delete mytimer;
}

void Gmsh_FileIO::check_FSI_ordering( const std::string &phy1,
   const std::string &phy2 ) const
{
  SYS_T::print_fatal_if( num_phy_domain_3d != 2, "Error: Gmsh_FileIO FSI mesh should contain only 2 physical domains.\n");

  const std::string name0 = phy_3d_name[ 0 ];
  const std::string name1 = phy_3d_name[ 1 ];

  SYS_T::print_fatal_if( name0.compare(phy1), "Error: Gmsh_FileIO FSI mesh 3d subdomain index 0 should be fluid domain.\n" );
  SYS_T::print_fatal_if( name1.compare(phy2), "Error: Gmsh_FileIO FSI mesh 3d subdomain index 1 should be solid domain.\n" );
}

void Gmsh_FileIO::write_tri_h5( const int &index_2d, 
    const std::vector<int> &index_1d ) const
{
  // Perform basic logical checks
  SYS_T::print_fatal_if( index_2d >= num_phy_domain_2d || index_2d < 0,
      "Error: Gmsh_FileIO::write_2d_h5, surface index is wrong. \n");

  const unsigned int num_1d_edge = index_1d.size();

  for(unsigned int ii=0; ii<num_1d_edge; ++ii)
    SYS_T::print_fatal_if( index_1d[ii] >= num_phy_domain_1d || index_1d[ii] < 0,
              "Error: Gmsh_FileIO::write_2d_h5, edge index is wrong. \n");

  // Open an HDF5 file
  std::string h5_file_name( "Gmsh_" );
  h5_file_name.append( phy_2d_name[index_2d] );
  h5_file_name.append(".h5");

  std::cout<<"=== Gmsh_FileIO::write_tri_h5 for "
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
    
    const int bcnumpt = static_cast<int>( bcpt.size() );
    
    std::cout<<"      num of bc pt = "<<bcnumpt<<'\n';    
   
    // tript stores the coordinates of the boundary points 
    std::vector<double> tript; tript.clear(); tript.resize(3*bcnumpt);
    for( int jj=0; jj<bcnumpt; ++jj )
    {
      tript[jj*3]   = node[bcpt[jj]*3] ;
      tript[jj*3+1] = node[bcpt[jj]*3+1] ;
      tript[jj*3+2] = node[bcpt[jj]*3+2] ;
    } 

    // generate a mapper that maps the bc node to 1; other node to 0
    bool * bcmap = new bool [num_node]; 
    for(int jj=0; jj<num_node; ++jj) bcmap[jj] = 0;
    for(int jj=0; jj<bcnumpt; ++jj) bcmap[bcpt[jj]] = 1;

    // generate a list of surface elements that has boundary on this edge 
    std::vector<int> gelem; gelem.clear();
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
    std::vector<int> edge_ien_local; edge_ien_local.clear();
    for(int ee=0; ee<num_1d_cell; ++ee)
    {
      for(int jj=0; jj<nLocBas_1d; ++jj)
      {
        std::vector<int>::iterator it = find( bcpt.begin(),
            bcpt.end(), edge_ien_global[nLocBas_1d*ee+jj] );
        edge_ien_local.push_back( it - bcpt.begin() );
      }
    }
    std::cout<<"      edge IEN generated. \n";

    // Locate the surface element that the edge belongs to. In Gmsh,
    // all elements are defined in this way. The edge two end points are
    // the first two points in the IEN array, no mater what the order is.
    // The triangle's three corner points are the first three points in the
    // triangle's IEN, no mather how many interior points there are.
    // Hence, we only need to treat all elements as if they are linear
    // line/triangle elements. Once the end points match with the corner
    // points, the edge is on the triangle boundary.
    std::vector<int> face2elem; face2elem.resize(num_1d_cell, -1);
    int sur_elem, node0, node1;
    int snode[3];
    bool got0, got1, gotit;
    for(int ff=0; ff<num_1d_cell; ++ff)
    {
      node0 = edge_ien_global[ nLocBas_1d * ff + 0 ];
      node1 = edge_ien_global[ nLocBas_1d * ff + 1 ];
      gotit = false;
      int ee = -1;
      while( !gotit && ee < int(gelem.size()) - 1 )
      {
        ee += 1;
        sur_elem = gelem[ee];
        snode[0] = eIEN[domain_2d_idx][nLocBas_2d * sur_elem + 0];
        snode[1] = eIEN[domain_2d_idx][nLocBas_2d * sur_elem + 1];
        snode[2] = eIEN[domain_2d_idx][nLocBas_2d * sur_elem + 2];
        std::sort(snode, snode+3);

        got0 = ( std::find(snode, snode+3, node0) != snode+3 );
        got1 = ( std::find(snode, snode+3, node1) != snode+3 );
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
    h5w->write_doubleVector(g_id, "pt_coor", tript);

    H5Gclose(g_id);
  }

  // Close the HDF5 file
  delete h5w;
  H5Fclose( file_id );
}

void Gmsh_FileIO::write_tet_h5( const int &index_3d,
    const std::vector<int> &index_2d ) const
{
  // Perform basic index boundary check
  SYS_T::print_fatal_if( index_3d >= num_phy_domain_3d || index_3d < 0,
      "Error: Gmsh_FileIO::write_tet_h5, surface index is wrong.\n");

  const unsigned int num_2d_face = index_2d.size();

  for(unsigned int ii=0; ii<num_2d_face; ++ii)
    SYS_T::print_fatal_if( index_2d[ii] >= num_phy_domain_2d || index_2d[ii] < 0,
        "Error: Gmsh_FileIO::write_tet_h5, face index is wrong. \n");

  // Open an HDF5 file
  std::string h5_file_name( "Gmsh_" );
  h5_file_name.append( phy_3d_name[index_3d] );
  h5_file_name.append(".h5");

  std::cout<<"=== Gmsh_FileIO::write_tet_h5 for "
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
    const int bcnumpt = static_cast<int>( bcpt.size() );
    std::cout<<"      num of bc pt = "<<bcnumpt<<'\n';
    
    // tript stores the coordinates of the boundary points
    std::vector<double> tript; tript.clear(); tript.resize(3*bcnumpt);
    for( int jj=0; jj<bcnumpt; ++jj )
    {
      tript[jj*3]   = node[bcpt[jj]*3];
      tript[jj*3+1] = node[bcpt[jj]*3+1];
      tript[jj*3+2] = node[bcpt[jj]*3+2];
    }
    
    // Generate a mapper that maps the bc node to 1; other node to 0
    bool * bcmap = new bool [num_node];
    for(int jj=0; jj<num_node; ++jj) bcmap[jj] = 0;
    for(int jj=0; jj<bcnumpt; ++jj) bcmap[bcpt[jj]] = 1;

    // Generate a list of volume elements that have boundary over the face
    std::vector<int> gelem; gelem.clear();
    for( int ee=0; ee<num_3d_cell; ++ee )
    {
      int total = 0;
      for(int jj=0; jj<nLocBas_3d; ++jj)
        total += bcmap[ eIEN[domain_3d_idx][nLocBas_3d*ee+jj] ];

      if(total >= 3) gelem.push_back(ee);
    }
    delete [] bcmap; bcmap = nullptr;
    std::cout<<"      "<<gelem.size()<<" elems have face over the boundary.\n";
    
    // Generate the local face element IEN array
    std::vector<int> face_ien_local; face_ien_local.clear();
    for(int ee=0; ee<num_2d_cell; ++ee)
    {
      for(int jj=0; jj<nLocBas_2d; ++jj)
      {
        std::vector<int>::iterator it = find( bcpt.begin(),
            bcpt.end(), face_ien_global[nLocBas_2d*ee+jj] );
        face_ien_local.push_back( it - bcpt.begin() );
      }
    }
    std::cout<<"      edge IEN generated. \n";
     
    // Loacate the volumetric element that the face element belongs to 
    std::vector<int> face2elem; face2elem.resize(num_2d_cell, -1);
    int vol_elem, node0, node1, node2;
    int vnode[4];
    bool got0, got1, got2, gotit;
    for(int ff=0; ff<num_2d_cell; ++ff)
    {
      node0 = face_ien_global[ nLocBas_2d * ff + 0 ];
      node1 = face_ien_global[ nLocBas_2d * ff + 1 ];
      node2 = face_ien_global[ nLocBas_2d * ff + 2 ];
      gotit = false;
      int ee = -1;
      while( !gotit && ee < int(gelem.size()) - 1 )
      {
        ee += 1;
        vol_elem = gelem[ee];
        vnode[0] = eIEN[domain_3d_idx][nLocBas_3d * vol_elem + 0];
        vnode[1] = eIEN[domain_3d_idx][nLocBas_3d * vol_elem + 1];
        vnode[2] = eIEN[domain_3d_idx][nLocBas_3d * vol_elem + 2];
        vnode[3] = eIEN[domain_3d_idx][nLocBas_3d * vol_elem + 3];
        std::sort(vnode, vnode+4);
        
        got0 = ( std::find(vnode, vnode+4, node0) != vnode+4 );
        got1 = ( std::find(vnode, vnode+4, node1) != vnode+4 );
        got2 = ( std::find(vnode, vnode+4, node2) != vnode+4 );
        gotit = got0 && got1 && got2;
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
    h5w->write_doubleVector(g_id, "pt_coor", tript);

    H5Gclose(g_id);
  } // End-loop-over-2d-face

  // Close the HDF5 file
  delete h5w;
  H5Fclose( file_id ); 
}

void Gmsh_FileIO::write_tet_h5( const int &index_3d,
    const std::vector<int> &index_2d,
    const std::vector<int> &index_2d_need_facemap ) const
{
  // Perform basic index boundary check
  SYS_T::print_fatal_if( index_3d >= num_phy_domain_3d || index_3d < 0,
      "Error: Gmsh_FileIO::write_tet_h5, surface index is wrong.\n");

  const unsigned int num_2d_face = index_2d.size();

  for(unsigned int ii=0; ii<num_2d_face; ++ii)
    SYS_T::print_fatal_if( index_2d[ii] >= num_phy_domain_2d || index_2d[ii] < 0,
        "Error: Gmsh_FileIO::write_tet_h5, face index is wrong. \n");

  // Open an HDF5 file
  std::string h5_file_name( "Gmsh_" );
  h5_file_name.append( phy_3d_name[index_3d] );
  h5_file_name.append(".h5");

  std::cout<<"=== Gmsh_FileIO::write_tet_h5 for "
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
    const int bcnumpt = static_cast<int>( bcpt.size() );
    std::cout<<"      num of bc pt = "<<bcnumpt<<'\n';
    
    // tript stores the coordinates of the boundary points
    std::vector<double> tript; tript.clear(); tript.resize(3*bcnumpt);
    for( int jj=0; jj<bcnumpt; ++jj )
    {
      tript[jj*3]   = node[bcpt[jj]*3];
      tript[jj*3+1] = node[bcpt[jj]*3+1];
      tript[jj*3+2] = node[bcpt[jj]*3+2];
    }
    
    // Generate a mapper that maps the bc node to 1; other node to 0
    bool * bcmap = new bool [num_node];
    for(int jj=0; jj<num_node; ++jj) bcmap[jj] = 0;
    for(int jj=0; jj<bcnumpt; ++jj) bcmap[bcpt[jj]] = 1;

    // Generate a list of volume elements that have boundary over the face
    std::vector<int> gelem; gelem.clear();
    for( int ee=0; ee<num_3d_cell; ++ee )
    {
      int total = 0;
      for(int jj=0; jj<nLocBas_3d; ++jj)
        total += bcmap[ eIEN[domain_3d_idx][nLocBas_3d*ee+jj] ];

      if(total >= 3) gelem.push_back(ee);
    }
    delete [] bcmap; bcmap = nullptr;
    std::cout<<"      "<<gelem.size()<<" elems have face over the boundary.\n";
    
    // Generate the local face element IEN array
    std::vector<int> face_ien_local; face_ien_local.clear();
    for(int ee=0; ee<num_2d_cell; ++ee)
    {
      for(int jj=0; jj<nLocBas_2d; ++jj)
      {
        std::vector<int>::iterator it = find( bcpt.begin(),
            bcpt.end(), face_ien_global[nLocBas_2d*ee+jj] );
        face_ien_local.push_back( it - bcpt.begin() );
      }
    }
    std::cout<<"      edge IEN generated. \n";
     
    // Loacate the volumetric element that the face element belongs to 
    std::vector<int> face2elem; face2elem.resize(num_2d_cell, -1);
    if( VEC_T::is_invec( index_2d_need_facemap, index_2d[ii] ) )
    {
      int vol_elem, node0, node1, node2;
      int vnode[4];
      bool got0, got1, got2, gotit;
      for(int ff=0; ff<num_2d_cell; ++ff)
      {
        node0 = face_ien_global[ nLocBas_2d * ff + 0 ];
        node1 = face_ien_global[ nLocBas_2d * ff + 1 ];
        node2 = face_ien_global[ nLocBas_2d * ff + 2 ];
        gotit = false;
        int ee = -1;
        while( !gotit && ee < int(gelem.size()) - 1 )
        {
          ee += 1;
          vol_elem = gelem[ee];
          vnode[0] = eIEN[domain_3d_idx][nLocBas_3d * vol_elem + 0];
          vnode[1] = eIEN[domain_3d_idx][nLocBas_3d * vol_elem + 1];
          vnode[2] = eIEN[domain_3d_idx][nLocBas_3d * vol_elem + 2];
          vnode[3] = eIEN[domain_3d_idx][nLocBas_3d * vol_elem + 3];
          std::sort(vnode, vnode+4);

          got0 = ( std::find(vnode, vnode+4, node0) != vnode+4 );
          got1 = ( std::find(vnode, vnode+4, node1) != vnode+4 );
          got2 = ( std::find(vnode, vnode+4, node2) != vnode+4 );
          gotit = got0 && got1 && got2;
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
    h5w->write_doubleVector(g_id, "pt_coor", tript);

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
    if( !VEC_T::is_invec(new2old, ii) ) new2old.push_back( ii );
  }

  // Now clean the snode vector to save memory
  VEC_T::clean( snode );

  SYS_T::print_fatal_if( static_cast<int>( new2old.size() ) != num_node, "Error: Gmsh_FildIO::update_FSI_nodal_ordering the number of nodes in the first two sub-domain does match the num_node!\n" );

  // Now generate the old2new mapping
  std::vector<int> old2new;
  old2new.resize( num_node );

  for(int ii=0; ii<num_node; ++ii) old2new[ new2old[ii] ] = ii;

  // Now clean the new2old mapper to save memory
  VEC_T::clean( new2old );

  // Now update the nodal x-y-z coordinates
  std::vector<double> temp; // temporary xyz coordinates based on new indices
  
  temp.resize( 3 * num_node );

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
    const int len = static_cast<int>( eIEN[ii].size() );
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

  SYS_T::print_fatal_if(nlocbas != 10, "Error: Gmsh_FileIO updata_quadratic_tet_IEN only works for 10-node quadratic element. \n");
  
  std::cout<<"=== Gmsh_FileIO::update_quadratic_tet_IEN for "<<phy_3d_name[index_3d]<<std::endl;;

  // Now upateing the eIEN array
  for(int ee=0; ee<ne; ++ee)
  {
    const double temp = eIEN[domain_index][10*ee+8];
    eIEN[domain_index][10*ee+8] = eIEN[domain_index][10*ee+9];
    eIEN[domain_index][10*ee+9] = temp;
  }
}

void Gmsh_FileIO::write_quadratic_sur_vtu( const int &index_sur,
    const int &index_vol, const bool &isf2e ) const
{
  SYS_T::print_fatal_if( index_sur >= num_phy_domain_2d || index_sur < 0,
      "Error: Gmsh_FileIO::write_vtp, surface index is wrong. \n");

  SYS_T::print_fatal_if( index_vol >= num_phy_domain_3d || index_vol < 0,
      "Error: Gmsh_FileIO::write_vtp, volume index is wrong. \n");

  const int phy_index_sur = phy_2d_index[index_sur];
  const int phy_index_vol = phy_3d_index[index_vol];
  const int bcnumcl = phy_2d_nElem[index_sur];
  const int numcel = phy_3d_nElem[index_vol];

  SYS_T::print_fatal_if( ele_nlocbas[phy_index_sur] != 6, "Error: Gmsh_FileIO write_quadratic_sur_vtu only works for 6-node triangle surface mesh.\n");
  
  SYS_T::print_fatal_if( ele_nlocbas[phy_index_vol] != 10, "Error: Gmsh_FileIO write_quadratic_sur_vtu only works for 10-node tetrahedral volumetric mesh.\n");

  std::cout<<"=== Gmsh_FileIO::write_quadratuc_sur_vtu for "
    <<phy_2d_name[index_sur]
    <<" associated with "<<phy_3d_name[index_vol];

  if( isf2e )
    std::cout<<" with face-to-volume element index. \n";
  else
    std::cout<<" without face-to-volume element index. \n";

  SYS_T::Timer * mytimer = new SYS_T::Timer();

  std::string vtu_file_name(phy_2d_name[index_sur]);
  vtu_file_name += "_";
  vtu_file_name += phy_3d_name[index_vol];
  std::cout<<"-----> write "<<vtu_file_name<<".vtu \n";
  
  mytimer->Reset(); mytimer->Start();

  // triangle mesh ien copied from eIEN
  std::vector<int> trien_global( eIEN[phy_index_sur] );
  
  // global node index
  std::vector<int> bcpt( trien_global );

  SYS_T::print_fatal_if( int(trien_global.size() ) != 6 * bcnumcl,
      "Error: Gmsh_FileIO::write_quadratic_sur_vtu, sur IEN size wrong. \n" );

  VEC_T::sort_unique_resize( bcpt ); // unique ascending order nodes
 
  const int bcnumpt = static_cast<int>( bcpt.size() );

  std::cout<<"      num of bc pt = "<<bcnumpt<<'\n';

  // tript stores the coordinates of the boundary points
  std::vector<double> tript; tript.clear(); tript.resize(3*bcnumpt);
  for( int ii=0; ii<bcnumpt; ++ii )
  {
    tript[ii*3]   = node[bcpt[ii]*3] ;
    tript[ii*3+1] = node[bcpt[ii]*3+1] ;
    tript[ii*3+2] = node[bcpt[ii]*3+2] ;
  }

  // A mapper that maps bc node to 1 other to 0
  bool * bcmap = new bool [num_node];
  for(int ii=0; ii<num_node; ++ii) bcmap[ii] = 0;
  for(int ii=0; ii<bcnumpt; ++ii) bcmap[bcpt[ii]] = 1;

  // Volume mesh IEN
  std::vector<int> vol_IEN( eIEN[phy_index_vol] );

  SYS_T::print_fatal_if( int( vol_IEN.size() ) != 10 * numcel,
      "Error: Gmsh_FileIO::write_quadratic_sur_vtu, vol IEN size wrong. \n");

  std::vector<int> gelem; gelem.clear();
  for( int ee=0; ee<numcel; ++ee )
  {
    int total = 0;
    total += bcmap[ vol_IEN[10*ee] ];
    total += bcmap[ vol_IEN[10*ee+1] ];
    total += bcmap[ vol_IEN[10*ee+2] ];
    total += bcmap[ vol_IEN[10*ee+3] ];
    if(total >= 3) gelem.push_back(ee);
  }
  delete [] bcmap; bcmap = nullptr;
  std::cout<<"      "<<gelem.size()<<" tets have faces over the surface. \n";

  // generate local triangle IEN array
  std::vector<int> trien; trien.clear();
  for(int ee=0; ee<bcnumcl; ++ee)
  {
    auto it = find(bcpt.begin(), bcpt.end(), trien_global[6*ee]);
    trien.push_back( it - bcpt.begin() );

    it = find(bcpt.begin(), bcpt.end(), trien_global[6*ee+1]);
    trien.push_back( it - bcpt.begin() );

    it = find(bcpt.begin(), bcpt.end(), trien_global[6*ee+2]);
    trien.push_back( it - bcpt.begin() );
    
    it = find(bcpt.begin(), bcpt.end(), trien_global[6*ee+3]);
    trien.push_back( it - bcpt.begin() );

    it = find(bcpt.begin(), bcpt.end(), trien_global[6*ee+4]);
    trien.push_back( it - bcpt.begin() );

    it = find(bcpt.begin(), bcpt.end(), trien_global[6*ee+5]);
    trien.push_back( it - bcpt.begin() );
  }
  std::cout<<"      triangle IEN generated. \n"; 

  std::vector<int> face2elem; face2elem.resize( bcnumcl, -1 );
  if( isf2e )
  {
    int vnode[4];
    bool got0, got1, got2, gotit;
    for(int ff=0; ff<bcnumcl; ++ff)
    {
      const int node0 = trien_global[6*ff];
      const int node1 = trien_global[6*ff+1];
      const int node2 = trien_global[6*ff+2];
      gotit = false;
      int ee = -1;
      while( !gotit && ee < int(gelem.size()) - 1 )
      {
        ee += 1;
        const int vol_elem = gelem[ee];
        vnode[0] = vol_IEN[10*vol_elem];
        vnode[1] = vol_IEN[10*vol_elem+1];
        vnode[2] = vol_IEN[10*vol_elem+2];
        vnode[3] = vol_IEN[10*vol_elem+3];
        std::sort(vnode, vnode+4);

        got0 = ( std::find(vnode, vnode+4, node0) != vnode+4 );
        got1 = ( std::find(vnode, vnode+4, node1) != vnode+4 );
        got2 = ( std::find(vnode, vnode+4, node2) != vnode+4 );
        gotit = got0 && got1 && got2;
      }

      // If the boundary surface element is not found,
      // we write -1 as the mapping value
      if(gotit)
        face2elem[ff] = gelem[ee] + phy_3d_start_index[index_vol];
      else
        face2elem[ff] = -1;
    }
    std::cout<<"      face2elem mapping generated. \n";
  }
  
  std::vector<DataVecStr<int>> input_vtk_data {};
  input_vtk_data.push_back({bcpt, "GlobalNodeID", AssociateObject::Node});
  input_vtk_data.push_back({face2elem, "GlobalElementID", AssociateObject::Cell});
  TET_T::write_quadratic_triangle_grid( vtu_file_name, bcnumpt, bcnumcl,
      tript, trien, input_vtk_data );

  mytimer->Stop();
  std::cout<<"      Time taken "<<mytimer->get_sec()<<" sec. \n";
  delete mytimer;
}

// EOF
