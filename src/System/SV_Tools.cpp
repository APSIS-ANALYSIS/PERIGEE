#include "SV_Tools.hpp"

void SV_T::update_sv_vtu( const std::string &filename,
    const std::string &write_name,
    int &sv_node_start, int &sv_elem_start )
{
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  VTK_T::read_vtu_grid(filename, nFunc, nElem, ctrlPts, vecIEN);

  // Read element index
  vtkXMLUnstructuredGridReader * reader = vtkXMLUnstructuredGridReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();
  vtkUnstructuredGrid * vtkugrid = reader -> GetOutput();

  vtkCellData * celldata = vtkugrid->GetCellData();
  vtkDataArray * cd = celldata->GetScalars("GlobalElementID");

  std::vector<int> eid; eid.clear();
  for(int ii=0; ii<nElem; ++ii)
    eid.push_back( static_cast<int>( cd->GetComponent(ii, 0) ) );

  vtkPointData * pointdata = vtkugrid->GetPointData();
  vtkDataArray * pd = pointdata->GetScalars("GlobalNodeID");

  std::vector<int> nid; nid.clear();
  for(int ii=0; ii<nFunc; ++ii)
    nid.push_back( static_cast<int>(pd->GetComponent(ii,0)) );

  reader->Delete();

  // Determine the starting indices as the minimum in nid / eid
  sv_node_start = *std::min_element( nid.begin(), nid.end() );
  sv_elem_start = *std::min_element( eid.begin(), eid.end() );

  // Write my file
  std::string fname(write_name);
  std::string fend;
  fend.assign( fname.end()-4 , fname.end() );

  // If the last four is .vtu, remove them for TET_T::write_tet_grid 
  if(fend.compare(".vtu") == 0) fname.erase( fname.end()-4, fname.end() );
  
  std::vector<DataVecStr<int>> input_vtk_data {};

  std::vector<int> temp_nid(nFunc, 0);
  for(int ii=0; ii<nFunc; ++ii) temp_nid[ii] = ii;
  input_vtk_data.push_back({temp_nid, "GlobalNodeID", AssociateObject::Node});

  std::vector<int> temp_eid(nElem, 0);
  for(int ii=0; ii<nElem; ++ii) temp_eid[ii] = ii;
  input_vtk_data.push_back({temp_eid, "GlobalElementID", AssociateObject::Cell});

  TET_T::write_tet_grid(fname, nFunc, nElem, ctrlPts, vecIEN, input_vtk_data);
}


void SV_T::gen_sv_fsi_vtus( const std::string &filename_f,
    const std::string &filename_s,
    const std::string &filename_f_wall,
    const std::string &writename_whole,
    const std::string &writename_solid,
    std::vector<int> &map_s_node, std::vector<int> &map_s_elem )
{
  // Merge data
  std::vector<int> wIEN; // whole mesh IEN
  std::vector<int> wtag; // whole mesh domain tag
  std::vector<double> ctrlPts;
  wIEN.clear(); wtag.clear(); ctrlPts.clear();

  // Read fluid mesh
  int nFunc_f, nElem_f;
  std::vector<int> vecIEN_f;
  std::vector<double> ctrlPts_f;

  VTK_T::read_vtu_grid(filename_f, nFunc_f, nElem_f, ctrlPts_f, vecIEN_f);

  // The FSI mesh IEN and control points start from those in fluid sub-domain
  VEC_T::insert_end( wIEN, vecIEN_f );
  VEC_T::clean( vecIEN_f );
  for(int ii=0; ii<nElem_f; ++ii) wtag.push_back(0);
  VEC_T::insert_end( ctrlPts, ctrlPts_f );
  VEC_T::clean( ctrlPts_f );

  // Read the interface mesh
  // This surface mesh is assumed to be the UPDATED fluid wall vtp mesh.
  int nFunc_i, nElem_i;
  std::vector<int> vecIEN_i;
  std::vector<double> ctrlPts_i;
  VTK_T::read_vtp_grid( filename_f_wall, nFunc_i, nElem_i, ctrlPts_i, vecIEN_i );
  std::vector<int> node_idx_i = VTK_T::read_int_PointData(filename_f_wall, "GlobalNodeID");

  // Read solid mesh
  int nFunc_s, nElem_s;
  std::vector<int> vecIEN_s;
  std::vector<double> ctrlPts_s;

  VTK_T::read_vtu_grid(filename_s, nFunc_s, nElem_s, ctrlPts_s, vecIEN_s);

  // Create a node mapping for the solid nodes that maps the solid own
  // node index to the whole FSI mesh nodal index.
  map_s_node.resize(nFunc_s);
  int counter = 0;
  for(int ii=0; ii<nFunc_s; ++ii)
  {
    const int loc = find_idx( ctrlPts_i, nFunc_i, ctrlPts_s[3*ii], 
        ctrlPts_s[3*ii+1], ctrlPts_s[3*ii+2], 1.0e-6 );

    if( loc == -1 )
    {
      map_s_node[ii] = nFunc_f + counter;
      ctrlPts.push_back( ctrlPts_s[3*ii] );
      ctrlPts.push_back( ctrlPts_s[3*ii+1] );
      ctrlPts.push_back( ctrlPts_s[3*ii+2] );
      counter += 1;
    }
    else
      map_s_node[ii] = node_idx_i[loc]; 
  }

  SYS_T::print_fatal_if( counter != nFunc_s - nFunc_i, "Error: SV_T::merge_sv_vtus, there are points in the interface not located in the solid domain. \n");

  // Clean the interface data
  VEC_T::clean(ctrlPts_i); VEC_T::clean(node_idx_i);

  // adjust the solid IEN array in the FSI IEN array
  std::vector<int> tempIEN; tempIEN.resize( vecIEN_s.size() ); 
  for(unsigned int ii=0; ii<vecIEN_s.size(); ++ii)
    tempIEN[ii] = map_s_node[ vecIEN_s[ii] ]; 

  VEC_T::insert_end( wIEN, tempIEN );
  VEC_T::clean(tempIEN);

  // Generate the physical tag array 
  for(int ii=0; ii<nElem_s; ++ii) wtag.push_back(1);

  // Write FSI mesh file
  std::string fname( writename_whole );
  std::string fend;
  fend.assign( fname.end()-4 , fname.end() );

  // If the last four is .vtu, remove them for TET_T::write_tet_grid 
  if(fend.compare(".vtu") == 0) fname.erase( fname.end()-4, fname.end() );

  std::vector<DataVecStr<int>> input_vtk_data {};
  input_vtk_data.push_back({wtag, "Physics_tag", AssociateObject::Cell});
  
  std::vector<int> temp_nid(nFunc_f + counter, 0);
  for(int ii=0; ii<nFunc_f + counter; ++ii) temp_nid[ii] = ii;
  input_vtk_data.push_back({temp_nid, "GlobalNodeID", AssociateObject::Node});

  std::vector<int> temp_eid(nElem_f + nElem_s, 0);
  for(int ii=0; ii<nElem_f + nElem_s; ++ii) temp_eid[ii] = ii;
  input_vtk_data.push_back({temp_eid, "GlobalElementID", AssociateObject::Cell});
  
  TET_T::write_tet_grid( fname, nFunc_f + counter, nElem_f + nElem_s,
      ctrlPts, wIEN, input_vtk_data, true );

  std::cout<<"Status: "<<writename_whole<<" is generated. \n";

  // Write solid vtu file
  fname = writename_solid;
  fend.assign( fname.end()-4 , fname.end() );
  if(fend.compare(".vtu") == 0) fname.erase( fname.end()-4, fname.end() );

  map_s_elem.resize(nElem_s);
  for(int ii=0; ii<nElem_s; ++ii) map_s_elem[ii] = nElem_f + ii;

  input_vtk_data.clear();
  input_vtk_data.push_back({map_s_node, "GlobalNodeID", AssociateObject::Node});
  input_vtk_data.push_back({map_s_elem, "GlobalElementID", AssociateObject::Cell});
  TET_T::write_tet_grid( fname, nFunc_s, nElem_s, ctrlPts_s, vecIEN_s, input_vtk_data );

  std::cout<<"Status: "<<filename_s<<" is updated to "<<writename_solid<<'\n';
}


void SV_T::update_sv_vtp( const std::string &filename,
    const std::string &writename,
    const int &nstart, const int &estart )
{
  // Read the vtp file
  vtkXMLPolyDataReader * reader = vtkXMLPolyDataReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();

  vtkPolyData * polydata = reader -> GetOutput();
  int numpts = static_cast<int>( polydata -> GetNumberOfPoints() );
  int numcels = static_cast<int>( polydata -> GetNumberOfPolys() );

  vtkCellData * celldata = polydata->GetCellData();
  vtkDataArray * cd = celldata->GetScalars("GlobalElementID");

  vtkPointData * pointdata = polydata->GetPointData();
  vtkDataArray * pd = pointdata->GetScalars("GlobalNodeID");

  std::vector<double> pt;
  pt.clear();
  std::vector<int> global_node_index; 
  global_node_index.clear();
  double pt_xyz[3];
  for(int ii=0; ii<numpts; ++ii)
  {
    polydata -> GetPoint(ii, pt_xyz);
    pt.push_back(pt_xyz[0]);
    pt.push_back(pt_xyz[1]);
    pt.push_back(pt_xyz[2]);

    global_node_index.push_back( static_cast<int>(pd->GetComponent(ii,0)) );
  }

  std::vector<int> ien_array, global_ele_index;
  ien_array.clear(); global_ele_index.clear();
  for(int ii=0; ii<numcels; ++ii)
  {
    vtkCell * cell = polydata -> GetCell(ii);

    if( cell->GetCellType() != 5 ) SYS_T::print_fatal("Error: read_vtp_grid read a mesh with non triangle elements. \n");

    ien_array.push_back( static_cast<int>( cell->GetPointId(0) ) );
    ien_array.push_back( static_cast<int>( cell->GetPointId(1) ) );
    ien_array.push_back( static_cast<int>( cell->GetPointId(2) ) );

    global_ele_index.push_back( static_cast<int>( cd->GetComponent(ii, 0)));
  }

  reader->Delete();

  // Update the nodal and elemental indices
  for(unsigned int ii=0; ii<global_node_index.size(); ++ii)
    global_node_index[ii] -= nstart;

  for(unsigned int ii=0; ii<global_ele_index.size(); ++ii)
    global_ele_index[ii] -= estart;

  // File name generate
  std::string fname(writename);
  std::string fend;
  fend.assign( fname.end()-4 , fname.end() );

  if(fend.compare(".vtp") == 0)
    fname.erase(fname.end()-4, fname.end());

  std::vector<DataVecStr<int>> input_vtk_data {};
  input_vtk_data.push_back({global_node_index, "GlobalNodeID", AssociateObject::Node});
  input_vtk_data.push_back({global_ele_index, "GlobalElementID", AssociateObject::Cell});
  TET_T::write_triangle_grid( fname, numpts, numcels, pt, ien_array, input_vtk_data );
}


void SV_T::update_sv_sur_vtu( const std::string &filename,
    const std::string &writename,
    const int &nstart, const int &estart )
{
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  VTK_T::read_vtu_grid( filename, nFunc, nElem, ctrlPts, vecIEN );

  vtkXMLUnstructuredGridReader * reader = vtkXMLUnstructuredGridReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();
  vtkUnstructuredGrid * vtkugrid = reader -> GetOutput();

  vtkCellData * celldata = vtkugrid->GetCellData();
  vtkDataArray * cd = celldata->GetScalars("GlobalElementID");

  std::vector<int> eid; eid.clear();
  for(int ii=0; ii<nElem; ++ii)
    eid.push_back( static_cast<int>( cd->GetComponent(ii, 0) ) - estart );

  vtkPointData * pointdata = vtkugrid->GetPointData();
  vtkDataArray * pd = pointdata->GetScalars("GlobalNodeID");

  std::vector<int> nid; nid.clear();
  for(int ii=0; ii<nFunc; ++ii)
    nid.push_back( static_cast<int>(pd->GetComponent(ii,0)) - nstart );

  reader->Delete();

  std::string fname(writename);
  std::string fend;
  fend.assign( fname.end()-4 , fname.end() );

  if(fend.compare(".vtu") == 0) fname.erase(fname.end()-4, fname.end());
  
  std::vector<DataVecStr<int>> input_vtk_data {};
  input_vtk_data.push_back({nid, "GlobalNodeID", AssociateObject::Node});
  input_vtk_data.push_back({eid, "GlobalElementID", AssociateObject::Cell});
  TET_T::write_quadratic_triangle_grid( fname, nFunc, nElem, ctrlPts, vecIEN, input_vtk_data );
}


void SV_T::update_sv_vtp( const std::string &filename,
    const std::string &writename,
    const int &nstart, const int &estart,
    const std::vector<int> &nmap, const std::vector<int> &emap )
{
  vtkXMLPolyDataReader * reader = vtkXMLPolyDataReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();

  vtkPolyData * polydata = reader -> GetOutput();
  const int numpts = static_cast<int>( polydata -> GetNumberOfPoints() );
  const int numcels = static_cast<int>( polydata -> GetNumberOfPolys() );

  vtkCellData * celldata = polydata->GetCellData();
  vtkDataArray * cd = celldata->GetScalars("GlobalElementID");

  vtkPointData * pointdata = polydata->GetPointData();
  vtkDataArray * pd = pointdata->GetScalars("GlobalNodeID");

  std::vector<double> pt;
  pt.clear();
  std::vector<int> global_node_index; 
  global_node_index.clear();
  for(int ii=0; ii<numpts; ++ii)
  {
    double pt_xyz[3];
    polydata -> GetPoint(ii, pt_xyz);
    pt.push_back(pt_xyz[0]);
    pt.push_back(pt_xyz[1]);
    pt.push_back(pt_xyz[2]);

    global_node_index.push_back( static_cast<int>(pd->GetComponent(ii,0)) );
  }

  std::vector<int> ien_array, global_ele_index;
  ien_array.clear(); global_ele_index.clear();
  for(int ii=0; ii<numcels; ++ii)
  {
    vtkCell * cell = polydata -> GetCell(ii);

    if( cell->GetCellType() != 5 ) SYS_T::print_fatal("Error: read_vtp_grid read a mesh with non triangle elements. \n");

    ien_array.push_back( static_cast<int>( cell->GetPointId(0) ) );
    ien_array.push_back( static_cast<int>( cell->GetPointId(1) ) );
    ien_array.push_back( static_cast<int>( cell->GetPointId(2) ) );

    global_ele_index.push_back( static_cast<int>( cd->GetComponent(ii, 0)));
  }

  reader->Delete();

  // Update the nodal and elemental indices
  // minus 1 because the files are from SimVascular
  for(unsigned int ii=0; ii<global_node_index.size(); ++ii)
    global_node_index[ii] = nmap[ global_node_index[ii] - nstart ];

  for(unsigned int ii=0; ii<global_ele_index.size(); ++ii)
    global_ele_index[ii] = emap[ global_ele_index[ii] - estart ];

  // File name generate
  std::string fname(writename);
  std::string fend;
  fend.assign( fname.end()-4 , fname.end() );

  if(fend.compare(".vtp") == 0)
    fname.erase(fname.end()-4, fname.end());

  std::vector<DataVecStr<int>> input_vtk_data {};
  input_vtk_data.push_back({global_node_index, "GlobalNodeID", AssociateObject::Node});
  input_vtk_data.push_back({global_ele_index, "GlobalElementID", AssociateObject::Cell});
  TET_T::write_triangle_grid( fname, numpts, numcels, pt, ien_array, input_vtk_data );
}


int SV_T::find_idx( const std::vector<double> &pt, const int &len,
    const double &x, const double &y, const double &z,
    const double tol )
{
  for(int ii=0; ii<len; ++ii)
  {
    if( MATH_T::equals(x, pt[3*ii], tol) && MATH_T::equals(y, pt[3*ii+1], tol)
        && MATH_T::equals(z, pt[3*ii+2], tol) ) return ii;
  }
  return -1;
}


void SV_T::compare_sv_vtp( const std::string &filename_1,
    const std::string &filename_2 )
{
  int numpts_1, numcels_1, numpts_2, numcels_2;
  std::vector<double> pt_1, pt_2;
  std::vector<int> ien_1, ien_2;

  VTK_T::read_vtp_grid( filename_1, numpts_1, numcels_1, pt_1, ien_1 );
  VTK_T::read_vtp_grid( filename_2, numpts_2, numcels_2, pt_2, ien_2 );

  SYS_T::print_fatal_if(numpts_1 != numpts_2, "Error: SV_T::compare_sv_vtp number of points does not match. \n");

  SYS_T::print_fatal_if(numcels_1 != numcels_2, "Error: SV_T::compare_sv_vtp number of cells does not match. \n");

  const int numpts = numpts_1;

  std::vector<double> pt2_x(numpts), pt2_y(numpts), pt2_z(numpts);
  for(int ii=0; ii<numpts; ++ii)
  {
    pt2_x[ii] = pt_2[3*ii];
    pt2_y[ii] = pt_2[3*ii+1];
    pt2_z[ii] = pt_2[3*ii+2];
  }

  // Check the points coordinates are matched.
  std::vector<int> map_idx; map_idx.resize(numpts);
  for(int ii=0; ii<numpts; ++ii)
  {
    const int loc = find_idx( pt_1, numpts, pt2_x[ii], pt2_y[ii], pt2_z[ii], 1.0e-15 );
    SYS_T::print_fatal_if( loc == -1, "Error: SV_T::compare_sv_vtp, There are points not found in the vtp file.\n" );
    map_idx[ii] = loc;
  }
}

// EOF
