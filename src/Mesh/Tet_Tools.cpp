#include "Tet_Tools.hpp"

void TET_T::gen_tet_grid( vtkUnstructuredGrid * const &grid_w,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array )
{
  // Check the input data compatibility
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: TET_T::gen_tet_grid point vector size does not match the number of points. \n");

  // detect the element type
  int nlocbas = -1;
  if( int(ien_array.size()) == 4*numcels ) nlocbas = 4;
  else if( int(ien_array.size()) == 10*numcels ) nlocbas = 10;
  else SYS_T::print_fatal("Error: TET_T::gen_tet_grid ien array size does not match the number of cells. \n");

  // Check the connectivity array
  std::vector<int> temp = ien_array;
  VEC_T::sort_unique_resize(temp);
  if( int(temp.size()) != numpts ) SYS_T::print_fatal("Error: TET_T::gen_tet_grid numpts does not match the number of unique points in the ien array. Please re-organize the input. \n");
  VEC_T::clean(temp);

  // 1. nodal points coordinates
  vtkPoints * ppt = vtkPoints::New(); 
  ppt->SetDataTypeToDouble();
  for(int ii=0; ii<numpts; ++ii)
  {
    const double coor[3] = { pt[3*ii], pt[3*ii+1], pt[3*ii+2] };
    ppt -> InsertPoint(ii, coor);
  }

  grid_w -> SetPoints(ppt);
  ppt -> Delete();

  // 2. Cell
  // Assign an additional data type for element aspect ratio
  vtkDoubleArray * edge_aspect_ratio = vtkDoubleArray::New();
  edge_aspect_ratio -> SetName("Aspect_ratio");
  edge_aspect_ratio -> SetNumberOfComponents(1);

  if(nlocbas == 4) 
  {
    vtkCell * cl = vtkTetra::New();

    for(int ii=0; ii<numcels; ++ii)
    {
      for(int lnode=0; lnode<4; ++lnode)
        cl->GetPointIds()->SetId( lnode, ien_array[4*ii + lnode] );

      grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );

      std::vector<double> cell_node; cell_node.clear();
      for(int lnode=0; lnode<4; ++lnode)
      {
        int node_offset = 3 * ien_array[4*ii + lnode];
        cell_node.push_back(pt[node_offset]);
        cell_node.push_back(pt[node_offset+1]);
        cell_node.push_back(pt[node_offset+2]);
      }
      edge_aspect_ratio -> InsertNextValue( TET_T::get_aspect_ratio(cell_node) );
    }

    cl -> Delete();
  }
  else if(nlocbas == 10)
  {
    vtkCell * cl = vtkQuadraticTetra::New();

    for(int ii=0; ii<numcels; ++ii)
    {
      for(int lnode=0; lnode<10; ++lnode)
        cl->GetPointIds()->SetId( lnode, ien_array[10*ii + lnode] );

      grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );

      std::vector<double> cell_node; cell_node.clear();
      for(int lnode=0; lnode<4; ++lnode)
      {
        int node_offset = 3 * ien_array[10*ii + lnode];
        cell_node.push_back(pt[node_offset]);
        cell_node.push_back(pt[node_offset+1]);
        cell_node.push_back(pt[node_offset+2]);
      }
      edge_aspect_ratio -> InsertNextValue( TET_T::get_aspect_ratio(cell_node) );
    }

    cl -> Delete();
  }
  else SYS_T::print_fatal("Error: TET_T::gen_tet_grid unknown local basis number.\n");

  // Add the asepct-ratio to grid_w
  grid_w -> GetCellData() -> AddArray( edge_aspect_ratio );
  edge_aspect_ratio -> Delete();
}

void TET_T::write_tet_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, const std::vector<int> &ien_array,
      const std::vector< DataVecStr<int> > &IOdata, const bool &isXML )
{
  // Setup the VTK objects
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  // Generate the mesh and compute aspect ratios
  gen_tet_grid( grid_w, numpts, numcels, pt, ien_array );

  // We need to make sure there are no data in IOdata that have the same name
  std::vector<std::string> name_list {};
  for( auto data : IOdata ) name_list.push_back( data.get_name() );

  VEC_T::sort_unique_resize( name_list );

  SYS_T::print_exit_if( name_list.size() != IOdata.size(), "Error: In TET_T::write_tet_grid, there are %d data in the IOdata that have the same name.\n", IOdata.size() - name_list.size() + 1 );
  
  // We add the IOdata for VTK
  for( auto data : IOdata )
  {
    if( data.get_object() == AssociateObject::Node )
      VTK_T::add_int_PointData( grid_w, data.get_data(), data.get_name() );
    else if( data.get_object() == AssociateObject::Cell )
      VTK_T::add_int_CellData( grid_w, data.get_data(), data.get_name() );
    else
      SYS_T::print_exit( "Error: In TET_T::write_tet_grid, there is an unknown object type in DataVecStr %s", data.get_name().c_str() );
  }

  // Write the prepared grid_w to a vtu or vtk file on disk
  VTK_T::write_vtkPointSet(filename, grid_w, isXML);

  grid_w->Delete();
}

void TET_T::write_triangle_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, 
      const std::vector<int> &ien_array,
      const std::vector<DataVecStr<int>> &IOdata )
{
  // Setup the VTK objects
  vtkPolyData * grid_w = vtkPolyData::New();

  // Generate the mesh
  gen_triangle_grid( grid_w, numpts, numcels, pt, ien_array );

  // We need to make sure there are no data in IOdata that have the same name
  std::vector<std::string> name_list {};
  for( auto data : IOdata ) name_list.push_back( data.get_name() );

  VEC_T::sort_unique_resize( name_list );

  SYS_T::print_exit_if( name_list.size() != IOdata.size(), "Error: In TET_T::write_triangle_grid, there are %d data in the IOdata that have the same name.\n", IOdata.size() - name_list.size() + 1 );

// We add the IOdata for VTK
  for( auto data : IOdata )
  {
    if( data.get_object() == AssociateObject::Node )
      VTK_T::add_int_PointData( grid_w, data.get_data(), data.get_name() );
    else if( data.get_object() == AssociateObject::Cell )
      VTK_T::add_int_CellData( grid_w, data.get_data(), data.get_name() );
    else
      SYS_T::print_exit( "Error: In TET_T::write_triangle_grid, there is an unknown object type in DataVecStr %s", data.get_name().c_str() );
  }

  // Write the prepared grid_w to a vtu or vtk file on disk
  VTK_T::write_vtkPointSet(filename, grid_w);

  grid_w->Delete();
}

void TET_T::gen_triangle_grid( vtkPolyData * const &grid_w,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array )
{
  // check the input data compatibility
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: TET_T::gen_triangle_grid point vector size does not match the number of points. \n");

  if(int(ien_array.size()) != 3*numcels) SYS_T::print_fatal("Error: TET_T::gen_triangle_grid ien array size does not match the number of cells. \n");

  // 1. nodal points
  vtkPoints * ppt = vtkPoints::New();
  ppt->SetDataTypeToDouble();
  for(int ii=0; ii<numpts; ++ii)
  {
    const double coor[3] = { pt[3*ii], pt[3*ii+1], pt[3*ii+2] };
    ppt -> InsertPoint(ii, coor);
  }

  grid_w -> SetPoints(ppt);
  ppt -> Delete();

  // 2. Cell
  vtkCellArray * cl = vtkCellArray::New();
  for(int ii=0; ii<numcels; ++ii)
  {
    vtkTriangle * tr = vtkTriangle::New();

    tr->GetPointIds()->SetId( 0, ien_array[3*ii] );
    tr->GetPointIds()->SetId( 1, ien_array[3*ii+1] );
    tr->GetPointIds()->SetId( 2, ien_array[3*ii+2] );
    cl -> InsertNextCell(tr);
    tr -> Delete();
  }
  grid_w->SetPolys(cl);
  cl->Delete();
}

void TET_T::write_quadratic_triangle_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, 
      const std::vector<int> &ien_array,
      const std::vector<DataVecStr<int>> &IOdata )
{
  // Setup the VTK objects
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  // Generate the mesh
  gen_quadratic_triangle_grid( grid_w, numpts, numcels, pt, ien_array );

  // We need to make sure there are no data in IOdata that have the same name
  std::vector<std::string> name_list {};
  for( auto data : IOdata ) name_list.push_back( data.get_name() );

  VEC_T::sort_unique_resize( name_list );

  SYS_T::print_exit_if( name_list.size() != IOdata.size(), "Error: In TET_T::write_triangle_grid, there are %d data in the IOdata that have the same name.\n", IOdata.size() - name_list.size() + 1 );

// We add the IOdata for VTK
  for( auto data : IOdata )
  {
    if( data.get_object() == AssociateObject::Node )
      VTK_T::add_int_PointData( grid_w, data.get_data(), data.get_name() );
    else if( data.get_object() == AssociateObject::Cell )
      VTK_T::add_int_CellData( grid_w, data.get_data(), data.get_name() );
    else
      SYS_T::print_exit( "Error: In TET_T::write_triangle_grid, there is an unknown object type in DataVecStr %s", data.get_name().c_str() );
  }

  // Write the prepared grid_w to a vtu file on disk
  VTK_T::write_vtkPointSet(filename, grid_w);

  grid_w->Delete();
}


void TET_T::gen_quadratic_triangle_grid( vtkUnstructuredGrid * const &grid_w,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array )
{
  // check the input data compatibility
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: TET_T::gen_quadratic_triangle_grid point vector size does not match the number of points. \n");

  if(int(ien_array.size()) != 6*numcels) SYS_T::print_fatal("Error: TET_T::gen_quadratic_quadratic_triangle_grid ien array size does not match the number of cells. \n");

  // 1. nodal points
  vtkPoints * ppt = vtkPoints::New(); 
  ppt->SetDataTypeToDouble();
  for(int ii=0; ii<numpts; ++ii)
  {
    const double coor[3] = { pt[3*ii], pt[3*ii+1], pt[3*ii+2] };
    ppt -> InsertPoint(ii, coor);
  }

  grid_w -> SetPoints(ppt);
  ppt -> Delete();

  // 2. Cell
  vtkCellArray * cl = vtkCellArray::New();

  for(int ii=0; ii<numcels; ++ii)
  {
    vtkQuadraticTriangle * tr = vtkQuadraticTriangle::New();

    tr->GetPointIds()->SetId( 0, ien_array[6*ii] );
    tr->GetPointIds()->SetId( 1, ien_array[6*ii+1] );
    tr->GetPointIds()->SetId( 2, ien_array[6*ii+2] );
    tr->GetPointIds()->SetId( 3, ien_array[6*ii+3] );
    tr->GetPointIds()->SetId( 4, ien_array[6*ii+4] );
    tr->GetPointIds()->SetId( 5, ien_array[6*ii+5] );
    cl -> InsertNextCell(tr);
    tr -> Delete();
  }

  grid_w -> SetCells(22, cl);
  cl -> Delete();
}

double TET_T::get_aspect_ratio( const std::vector<double> &pt )
{
  double edge[6];
  // e01
  edge[0] = MATH_T::norm2(pt[3]-pt[0], pt[4]-pt[1], pt[5]-pt[2]) ;
  // e12
  edge[1] = MATH_T::norm2(pt[6]-pt[3], pt[7]-pt[4], pt[8]-pt[5]) ;
  // e02
  edge[2] = MATH_T::norm2(pt[6]-pt[0], pt[7]-pt[1], pt[8]-pt[2]) ;
  // e03
  edge[3] = MATH_T::norm2(pt[9]-pt[0], pt[10]-pt[1], pt[11]-pt[2]) ;
  // e13
  edge[4] = MATH_T::norm2(pt[9]-pt[3], pt[10]-pt[4], pt[11]-pt[5]) ;
  // e23
  edge[5] = MATH_T::norm2(pt[9]-pt[6], pt[10]-pt[7], pt[11]-pt[8]) ;

  const double emax = *std::max_element(edge, edge+6);
  const double emin = *std::min_element(edge, edge+6);

  return emax / emin;
}

Vector_3 TET_T::get_out_normal( const std::string &file,
    const std::vector<double> &vol_ctrlPts,
    const IIEN * const &vol_ien )
{
  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien;

  // Analyze the file type
  std::string fend; fend.assign( file.end()-4 , file.end() );

  if( fend.compare(".vtp") == 0 )
    VTK_T::read_vtp_grid( file, numpts, numcels, pts, ien);
  else if( fend.compare(".vtu") == 0 )
    VTK_T::read_vtu_grid( file, numpts, numcels, pts, ien);
  else
    SYS_T::print_fatal("Error: get_out_normal unknown file type.\n");

  const std::vector<int> gnode = VTK_T::read_int_PointData(file, "GlobalNodeID");
  const std::vector<int> gelem = VTK_T::read_int_CellData(file, "GlobalElementID");
  
  // triangle nodes' global indices
  const std::vector<int> trn = { gnode[ ien[0] ], gnode[ ien[1] ], gnode[ ien[2] ] };

  const int tete0 = gelem[0]; // triangle's associated tet element indices

  SYS_T::print_fatal_if(tete0 == -1, "Error: TET_T::get_out_normal requires the element indices for the vtp file.\n");

  const int ten[4] { vol_ien->get_IEN(tete0, 0), vol_ien->get_IEN(tete0, 1), 
    vol_ien->get_IEN(tete0, 2), vol_ien->get_IEN(tete0, 3) };

  int inside_node = 0, node_check = 0;
  for(int ii=0; ii<4; ++ii)
  {
    if( !VEC_T::is_invec( trn, ten[ii] ) ) inside_node = ten[ii];
    else node_check += 1;
  }

  SYS_T::print_fatal_if(node_check!=3, "Error: TET_T::get_out_normal, the associated tet element is incompatible with the triangle element. \n");

  // make cross line-0-1 and line-0-2
  const Vector_3 l01( vol_ctrlPts[3*trn[1]] - vol_ctrlPts[3*trn[0]],
      vol_ctrlPts[3*trn[1]+1] - vol_ctrlPts[3*trn[0]+1],
      vol_ctrlPts[3*trn[1]+2] - vol_ctrlPts[3*trn[0]+2] );

  const Vector_3 l02( vol_ctrlPts[3*trn[2]] - vol_ctrlPts[3*trn[0]],
      vol_ctrlPts[3*trn[2]+1] - vol_ctrlPts[3*trn[0]+1],
      vol_ctrlPts[3*trn[2]+2] - vol_ctrlPts[3*trn[0]+2] );

  Vector_3 outVec = cross_product( l01, l02 );

  outVec.normalize();

  // obtain the line from tri node 0 to the out-of-surface node
  const Vector_3 inw( vol_ctrlPts[inside_node*3] - vol_ctrlPts[3*trn[0]],
      vol_ctrlPts[inside_node*3+1] - vol_ctrlPts[3*trn[0]+1],
      vol_ctrlPts[inside_node*3+2] - vol_ctrlPts[3*trn[0]+2] );

  // inner product outward with inward, and correct outVec
  if(inw.dot_product(outVec) > 0.0) outVec *= -1.0;

  return outVec;
}

void TET_T::tetgenio2vtu( const tetgenio &meshout, const std::string &fName )
{
  const int index_offset = meshout.firstnumber;
  const int numpts = meshout.numberofpoints;
  const int numcel = meshout.numberoftetrahedra;

  std::vector<double> pt;
  pt.resize( numpts * 3 );
  std::vector<int> ien;
  ien.resize( numcel * 4 );

  std::cout<<"Preparing for writing .vtu grid file ... \n";
  std::cout<<"  -- "<<numpts<<" mesh points detected. \n";
  std::cout<<"  -- "<<numcel<<" mesh tetrahedral cell detected. \n";

  for(int ii=0; ii<numpts*3; ++ii) pt[ii] = meshout.pointlist[ii];

  for(int ii=0; ii<numcel*4; ++ii) ien[ii] = meshout.tetrahedronlist[ii] - index_offset;

  std::vector<DataVecStr<int>> input_vtk_data {};

  std::vector<int> temp_nid(numpts, 0);
  for(int ii=0; ii<numpts; ++ii) temp_nid[ii] = ii;
  input_vtk_data.push_back({temp_nid, "GlobalNodeID", AssociateObject::Node});

  std::vector<int> temp_eid(numcel, 0);
  for(int ii=0; ii<numcel; ++ii) temp_eid[ii] = ii;
  input_vtk_data.push_back({temp_eid, "GlobalElementID", AssociateObject::Cell});

  write_tet_grid(fName, numpts, numcel, pt, ien, input_vtk_data);

  std::cout<<"Volumetric grid file "<<fName<<".vtu has been written on disk. \n";
}


void TET_T::tetgenio2vtp( const tetgenio &meshout, const std::string &fName,
    const int &bcmarker )
{
  const int index_offset = meshout.firstnumber;
  const int numpts = meshout.numberofpoints;
  const int numcel = meshout.numberoftetrahedra;

  // global point coordinates
  std::vector<double> pt;
  pt.resize( numpts * 3 );
  for(int ii=0; ii<numpts*3; ++ii) pt[ii] = meshout.pointlist[ii];

  // Get global volumetric IEN array
  std::vector<int> ien;
  ien.resize( numcel * 4 );
  for(int ii=0; ii<numcel*4; ++ii) ien[ii] = meshout.tetrahedronlist[ii] - index_offset;

  // Store the point on the boundary in bcpt
  std::cout<<"Preparing for writing boundary "<<bcmarker<<" mesh file \n";
  std::vector<int> bcpt;
  const int numface = meshout.numberoftrifaces;
  bcpt.clear();
  std::cout<<"  -- "<<numface<<" faces detected. \n";
  int bcnumcl = 0;
  std::vector<int> trien_global;
  trien_global.clear();
  for(int ee=0; ee<numface; ++ee)
  {
    if( meshout.trifacemarkerlist[ee] == bcmarker )
    {
      bcnumcl += 1;
      bcpt.push_back( meshout.trifacelist[ee*3] - index_offset );
      bcpt.push_back( meshout.trifacelist[ee*3+1] - index_offset );
      bcpt.push_back( meshout.trifacelist[ee*3+2] - index_offset );
    }
  }
  trien_global = bcpt;

  std::cout<<"  -- "<<bcnumcl<<" faces on BC "<<bcmarker<<" detected. \n";

  // list the bc node in unique ascending order
  VEC_T::sort_unique_resize(bcpt);
  const int bcnumpt = static_cast<int>( bcpt.size() );

  std::cout<<"  -- "<<bcnumpt<<" points on BC "<<bcmarker<<" detected. \n";

  // tript stores the coordinates of the boundary points 
  std::vector<double> tript; tript.clear(); tript.resize(3*bcpt.size());
  for( int ii=0; ii<bcnumpt; ++ii )
  {
    tript[ii*3]   = pt[bcpt[ii]*3] ;
    tript[ii*3+1] = pt[bcpt[ii]*3+1] ;
    tript[ii*3+2] = pt[bcpt[ii]*3+2] ;
  }

  // generate a mapper that maps the bc node to 1, other node to 0
  bool * bcmap = new bool [numpts];
  for(int ii=0; ii<numpts; ++ii) bcmap[ii] = 0;
  for(int ii=0; ii<bcnumpt; ++ii) bcmap[bcpt[ii]] = 1;

  std::vector<int> gelem; gelem.clear();
  for( int ee=0; ee<numcel; ++ee )
  {
    int total = 0;
    total += bcmap[ ien[4*ee] ];
    total += bcmap[ ien[4*ee+1] ];
    total += bcmap[ ien[4*ee+2] ];
    total += bcmap[ ien[4*ee+3] ];

    // get a list of volume element that has at least 3 poitns on
    // the boundary surface 
    if(total >= 3) gelem.push_back(ee);
  }
  delete [] bcmap;
  std::cout<<"  -- "<<gelem.size()<<" tets have face on this BC surface. \n";

  // generate the local triangle IEN array
  std::vector<int> trien; trien.clear(); 
  int node0, node1, node2;
  std::vector<int>::iterator it;
  for( int ee=0; ee<bcnumcl; ++ee)
  {
    node0 = trien_global[3*ee];
    node1 = trien_global[3*ee+1];
    node2 = trien_global[3*ee+2];

    // trien_global and bcpt stores the same items, just ordering is
    // changed. So, it should always find the node in bcpt list.
    it = find(bcpt.begin(), bcpt.end(), node0);
    trien.push_back( it - bcpt.begin() );

    it = find(bcpt.begin(), bcpt.end(), node1);
    trien.push_back( it - bcpt.begin() );

    it = find(bcpt.begin(), bcpt.end(), node2);
    trien.push_back( it - bcpt.begin() );
  }
  std::cout<<"  -- triangle IEN generated. \n"; 

  // determine the face-to-element mapping
  std::vector<int> face2elem; face2elem.resize(bcnumcl);
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
    while( !gotit )
    {
      ee += 1;
      vol_elem = gelem[ee];
      vnode[0] = ien[4*vol_elem];
      vnode[1] = ien[4*vol_elem+1];
      vnode[2] = ien[4*vol_elem+2];
      vnode[3] = ien[4*vol_elem+3];
      std::sort(vnode, vnode+4);
      got0 = ( std::find(vnode, vnode+4, node0) != vnode+4 );
      got1 = ( std::find(vnode, vnode+4, node1) != vnode+4 );
      got2 = ( std::find(vnode, vnode+4, node2) != vnode+4 );

      gotit = got0 && got1 && got2;
    }
    face2elem[ff] = gelem[ee];
  }

  // preparing for mesh file (.vtp) name
  std::string outname(fName);
  outname.append(".");
  std::string strbcindex = std::to_string(bcmarker);
  outname.append(strbcindex);

  std::vector<DataVecStr<int>> input_vtk_data {};
  input_vtk_data.push_back({bcpt, "GlobalNodeID", AssociateObject::Node});
  input_vtk_data.push_back({face2elem, "GlobalElementID", AssociateObject::Cell});

  write_triangle_grid( outname, bcnumpt, bcnumcl, tript, trien, input_vtk_data );
  std::cout<<"Surface file "<<outname<<".vtp has been written on disk. \n";  
}


namespace TET_T
{
  Tet4::Tet4()
  {
    pts[0] = 0.0; pts[1] = 0.0;  pts[2] = 0.0;
    pts[3] = 1.0; pts[4] = 0.0;  pts[5] = 0.0;
    pts[6] = 0.0; pts[7] = 1.0;  pts[8] = 0.0;
    pts[9] = 0.0; pts[10] = 0.0; pts[11] = 1.0;

    gindex[0] = 0; gindex[1] = 1;
    gindex[2] = 2; gindex[3] = 3;
  }

  Tet4::Tet4( const std::vector<double> &in_nodes )
  {
    SYS_T::print_exit_if( in_nodes.size() != 12, "Error: input nodal list shall have 4 nodes with xyz-coordinates. \n");

    for(int ii=0; ii<12; ++ii) pts[ii] = in_nodes[ii];

    gindex[0] = 0; gindex[1] = 1;
    gindex[2] = 2; gindex[3] = 3;
  }

  Tet4::Tet4( const std::vector<double> &ctrlPts,
      const int &ien0, const int &ien1,
      const int &ien2, const int &ien3 )
  {
    gindex[0] = ien0;
    gindex[1] = ien1;
    gindex[2] = ien2;
    gindex[3] = ien3;

    for(int ii=0; ii<4; ++ii)
    {
      pts[ii*3]   = ctrlPts[gindex[ii]*3];
      pts[ii*3+1] = ctrlPts[gindex[ii]*3+1];
      pts[ii*3+2] = ctrlPts[gindex[ii]*3+2];
    }
  }

  Tet4::~Tet4()
  {
  }

  void Tet4::reset( const std::vector<double> &in_nodes )
  {
    SYS_T::print_exit_if( in_nodes.size() != 12, "Error: input nodal list shall have 4 nodes with xyz-coordinates. \n");
    for(int ii=0; ii<12; ++ii) pts[ii] = in_nodes[ii];

    gindex[0] = 0; gindex[1] = 1;
    gindex[2] = 2; gindex[3] = 3;
  }

  void Tet4::reset( const std::vector<double> &ctrlPts,
      const int &ien0, const int &ien1,
      const int &ien2, const int &ien3 )
  {
    gindex[0] = ien0;
    gindex[1] = ien1;
    gindex[2] = ien2;
    gindex[3] = ien3;

    for(int ii=0; ii<4; ++ii)
    {
      pts[ii*3]   = ctrlPts[gindex[ii]*3];
      pts[ii*3+1] = ctrlPts[gindex[ii]*3+1];
      pts[ii*3+2] = ctrlPts[gindex[ii]*3+2];
    }
  }

  void Tet4::reset( const int &ien0, const int &ien1,
      const int &ien2, const int &ien3 )
  {
    gindex[0] = ien0;
    gindex[1] = ien1;
    gindex[2] = ien2;
    gindex[3] = ien3;
  }

  void Tet4::reset( const std::vector<double> &ctrlPts,
      const IIEN * const &ien_ptr, const int &ee )
  {
    gindex[0] = ien_ptr->get_IEN(ee, 0);
    gindex[1] = ien_ptr->get_IEN(ee, 1);
    gindex[2] = ien_ptr->get_IEN(ee, 2);
    gindex[3] = ien_ptr->get_IEN(ee, 3);

    for(int ii=0; ii<4; ++ii)
    {
      pts[ii*3]   = ctrlPts[gindex[ii]*3];
      pts[ii*3+1] = ctrlPts[gindex[ii]*3+1];
      pts[ii*3+2] = ctrlPts[gindex[ii]*3+2];
    }
  }

  void Tet4::print_info() const
  {
    std::cout<<"Tet4 object : \n";
    std::cout<<" -- global indices: "<<gindex[0]<<'\t';
    std::cout<<gindex[1]<<'\t'<<gindex[2]<<'\t'<<gindex[3]<<'\n';
    std::cout<<" -- pt0 : "<<pts[0]<<'\t'<<pts[1]<<'\t'<<pts[2]<<'\n';
    std::cout<<" -- pt1 : "<<pts[3]<<'\t'<<pts[4]<<'\t'<<pts[5]<<'\n';
    std::cout<<" -- pt2 : "<<pts[6]<<'\t'<<pts[7]<<'\t'<<pts[8]<<'\n';
    std::cout<<" -- pt3 : "<<pts[9]<<'\t'<<pts[10]<<'\t'<<pts[11]<<'\n';
  }

  double Tet4::get_aspect_ratio() const
  {
    double edge[6];
    // e01
    edge[0] = MATH_T::norm2(pts[3]-pts[0], pts[4]-pts[1], pts[5]-pts[2]) ;
    // e12
    edge[1] = MATH_T::norm2(pts[6]-pts[3], pts[7]-pts[4], pts[8]-pts[5]) ;
    // e02
    edge[2] = MATH_T::norm2(pts[6]-pts[0], pts[7]-pts[1], pts[8]-pts[2]) ;
    // e03
    edge[3] = MATH_T::norm2(pts[9]-pts[0], pts[10]-pts[1], pts[11]-pts[2]) ;
    // e13
    edge[4] = MATH_T::norm2(pts[9]-pts[3], pts[10]-pts[4], pts[11]-pts[5]) ;
    // e23
    edge[5] = MATH_T::norm2(pts[9]-pts[6], pts[10]-pts[7], pts[11]-pts[8]) ;

    const double emax = *std::max_element(edge, edge+6);
    const double emin = *std::min_element(edge, edge+6);

    return emax / emin;
  }

  int Tet4::get_face_id( const int &n0, const int &n1, const int &n2) const
  {
    const int temp[3] { n0, n1, n2 };

    const auto it0 = std::find(temp, temp+3, gindex[0]);
    const auto it1 = std::find(temp, temp+3, gindex[1]);
    const auto it2 = std::find(temp, temp+3, gindex[2]);
    const auto it3 = std::find(temp, temp+3, gindex[3]);

    // flg is true if the node is found in temp
    bool flg[4] { (it0 != temp+3), (it1 != temp+3), (it2 != temp+3), (it3 != temp+3) };

    int sum = flg[0] + flg[1] + flg[2] + flg[3];

    SYS_T::print_exit_if( sum != 3, "Error: Tet4::find_face_id input is not a proper face of this element. \n");

    int loc = 0;
    while( flg[loc] ) loc+=1;

    return loc;
  }

  double Tet4::get_diameter() const
  {
    double x, y, z, r;
    MATH_T::get_tet_sphere_info(
        pts[0], pts[3], pts[6], pts[9],
        pts[1], pts[4], pts[7], pts[10],
        pts[2], pts[5], pts[8], pts[11],
        x, y, z, r );

    return 2.0 * r;
  }

  double Tet4::get_volume() const
  {
    const double x10 = pts[3] - pts[0];
    const double y01 = pts[1] - pts[4];
    const double z01 = pts[2] - pts[5];

    const double x21 = pts[6] - pts[3];
    const double y12 = pts[4] - pts[7];
    const double z12 = pts[5] - pts[8];

    const double x32 = pts[9] - pts[6];
    const double z23 = pts[8] - pts[11];
    const double y23 = pts[7] - pts[10];

    const double sum = x10 * (y12*z23 - y23*z12)
      + x21*(y23*z01 - y01*z23) + x32*(y01*z12 - y12*z01);

    return sum / 6.0;
  }

  void Tet4::write_vtu( const std::string &fileName ) const
  {
    // Setup the VTK objects
    vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

    std::vector<double> input_pt(pts, pts+12);
    std::vector<int> input_ien{0, 1, 2, 3};

    gen_tet_grid( grid_w, 4, 1, input_pt, input_ien );

    std::vector<int> input_node_index(gindex, gindex+4);

    VTK_T::add_int_PointData( grid_w, input_node_index, "GlobalNodeID" );

    // write vtu
    VTK_T::write_vtkPointSet( fileName, grid_w, true );

    grid_w->Delete();
  }
}


void TET_T::tetmesh_check(const std::vector<double> &cpts,
    const IIEN * const &ienptr, const int &nelem,
    const double &crit_aspect_ratio )
{
  std::cout<<"\n======= Tet mesh quality =======\n";
  TET_T::Tet4 * teton = new TET_T::Tet4();
  teton -> reset( cpts, ienptr, 0 );
  double teton_max_vol = teton -> get_volume();
  double teton_min_vol = teton_max_vol;
  double teton_max_h = teton -> get_diameter();
  double teton_min_h = teton_max_h;
  int num_dist_elem = 0; // number of element that has negative volume
  int num_aspt_elem = 0; // number of element that has aspect ratio larger
  // than the given critical aspect ratio value
  for(int ee = 0; ee<nelem; ++ee)
  {
    teton->reset( cpts, ienptr, ee );

    if( teton->get_aspect_ratio() > crit_aspect_ratio )
      num_aspt_elem += 1;

    double teton_ee_vol = teton -> get_volume();
    if( teton_ee_vol < 0.0 ) num_dist_elem += 1;

    if( teton_max_vol < teton_ee_vol) teton_max_vol = teton_ee_vol;

    if( teton_min_vol > teton_ee_vol) teton_min_vol = teton_ee_vol;

    double teton_he = teton -> get_diameter();

    if( teton_max_h < teton_he ) teton_max_h = teton_he;

    if( teton_min_h > teton_he ) teton_min_h = teton_he;
  }
  delete teton;
  cout<<"- maximum tetrahedron volume : "<<teton_max_vol<<endl;
  cout<<"- minimum tetrahedron volume : "<<teton_min_vol<<endl;
  cout<<"- maximum tetrahedron diameter : "<<teton_max_h<<endl;
  cout<<"- minimum tetrahedron diameter : "<<teton_min_h<<endl;
  cout<<"- number of distorted element : "<<num_dist_elem<<endl;
  cout<<"- number of element with aspect ratio larger than "<<crit_aspect_ratio<<" : "<<num_aspt_elem<<endl;
  std::cout<<"==================================\n";
}

// EOF
