#include "NodalBC_3D_vtu.hpp"

NodalBC_3D_vtu::NodalBC_3D_vtu( const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_dir_nodes = 0;
  num_per_nodes = 0;

  Create_ID( nFunc );

  std::cout<<"===> NodalBC_3D_vtu: No nodal BC is generated. \n";
}

NodalBC_3D_vtu::NodalBC_3D_vtu( const std::string &vtufilename, 
    const int &nFunc )
{
  SYS_T::file_check( vtufilename );

  // ----- Read the vtu file
  vtkXMLUnstructuredGridReader * reader = vtkXMLUnstructuredGridReader::New();
  reader -> SetFileName( vtufilename.c_str() );
  reader -> Update();
  vtkUnstructuredGrid * vtkugrid = reader -> GetOutput();

  const int numpts = static_cast<int>( vtkugrid -> GetNumberOfPoints() );

  vtkPointData * pointdata = vtkugrid->GetPointData();
  vtkDataArray * pd = pointdata->GetScalars("GlobalNodeID");

  std::vector<unsigned int> gnode; gnode.clear();

  for(int ii=0; ii<numpts; ++ii) gnode.push_back( static_cast<unsigned int>(pd->GetComponent(ii,0)) );

  reader->Delete();

  VEC_T::sort_unique_resize( gnode );
  // ----- Finish the read

  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;
  num_dir_nodes = gnode.size();

  VEC_T::insert_end(dir_nodes, gnode);
  VEC_T::sort_unique_resize(dir_nodes);
  num_dir_nodes = dir_nodes.size();

  // Generate the ID array
  Create_ID(nFunc);
  std::cout<<"===> NodalBC_3D_vtu specified by "<<vtufilename<<" is generated. \n";
}


NodalBC_3D_vtu::NodalBC_3D_vtu( const std::vector<std::string> &vtufileList,
    const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  const unsigned int num_file = vtufileList.size();
  for(unsigned int ii=0; ii<num_file; ++ii)
  {
    SYS_T::file_check( vtufileList[ii] );

    int numpts, numcels;
    std::vector<double> pts;
    std::vector<int> ien;

    VTK_T::read_vtu_grid( vtufileList[ii], numpts, numcels, pts, ien );
 
    const std::vector<int> gnode = VTK_T::read_int_PointData( vtufileList[ii], "GlobalNodeID");
    const std::vector<int> gelem = VTK_T::read_int_CellData( vtufileList[ii], "GlobalElementID");

    for(unsigned int jj=0; jj<gnode.size(); ++jj)
    {
      if(gnode[jj]<0) SYS_T::print_fatal("Error: there are negative nodal index! \n");

      dir_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
    }
  }

  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  std::cout<<"===> NodalBC_3D_vtu specified by \n";
  for(unsigned int ii=0; ii<num_file; ++ii)
    std::cout<<"     "<<vtufileList[ii]<<std::endl;
  std::cout<<"     is generated. \n";
}


NodalBC_3D_vtu::NodalBC_3D_vtu( const std::string &vtufilename,
    const std::vector<std::string> &vtpfileList, const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  // Read vtp files
  const unsigned int num_file = vtpfileList.size();

  for(unsigned int ii=0; ii<num_file; ++ii)
  {
    SYS_T::file_check( vtpfileList[ii] );

    int numpts, numcels;
    std::vector<double> pts;
    std::vector<int> ien;

    VTK_T::read_vtp_grid( vtpfileList[ii], numpts, numcels, pts, ien );

    const std::vector<int> gnode = VTK_T::read_int_PointData( vtpfileList[ii], "GlobalNodeID");
    const std::vector<int> gelem = VTK_T::read_int_CellData( vtpfileList[ii], "GlobalElementID");

    if( numpts != static_cast<int>(gnode.size()) )
      SYS_T::print_fatal("Error: the numpts != global_node.size()! \n");

    for(unsigned int jj=0; jj<gnode.size(); ++jj)
    {
      if(gnode[jj]<0) SYS_T::print_fatal("Error: there are negative nodal index! \n");

      dir_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
    }
  }

  // Read vtu file
  SYS_T::file_check( vtufilename );

  vtkXMLUnstructuredGridReader * reader = vtkXMLUnstructuredGridReader::New();
  reader -> SetFileName( vtufilename.c_str() );
  reader -> Update();
  vtkUnstructuredGrid * vtkugrid = reader -> GetOutput();

  const int numpts = static_cast<int>( vtkugrid -> GetNumberOfPoints() );

  vtkPointData * pointdata = vtkugrid->GetPointData();
  vtkDataArray * pd = pointdata->GetScalars("GlobalNodeID");

  std::vector<unsigned int> gnode; gnode.clear();

  for(int ii=0; ii<numpts; ++ii) gnode.push_back( static_cast<unsigned int>(pd->GetComponent(ii,0)) );

  reader->Delete();

  VEC_T::sort_unique_resize( gnode );

  VEC_T::insert_end(dir_nodes, gnode);

  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  std::cout<<"===> NodalBC_3D_vtu specified by \n";
  for(unsigned int ii=0; ii<num_file; ++ii)
    std::cout<<"     "<<vtpfileList[ii]<<std::endl;
  std::cout<<"     "<<vtufilename<<"     is generated. \n";
}


NodalBC_3D_vtu::~NodalBC_3D_vtu()
{}


NodalBC_3D_vtu::NodalBC_3D_vtu( const std::string &vtufilename,
    const std::vector<std::string> &vtpfileList, 
    const std::vector<std::string> &vtpminus, 
    const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  // Read the minus vtp files first
  std::vector<unsigned int> minus_nodes;
  const unsigned int num_minus_file = vtpminus.size();
  minus_nodes.clear();

  for(unsigned int ii=0; ii<num_minus_file; ++ii)
  {
    SYS_T::file_check( vtpminus[ii] );

    int m_numpts, m_numcels;
    std::vector<double> m_pts;
    std::vector<int> m_ien;

    VTK_T::read_vtp_grid( vtpminus[ii], m_numpts, m_numcels, 
        m_pts, m_ien );
    const std::vector<int> m_gnode = VTK_T::read_int_PointData( vtpminus[ii], "GlobalNodeID");
    const std::vector<int> m_gelem = VTK_T::read_int_CellData( vtpminus[ii], "GlobalElementID");

    if( m_numpts != static_cast<int>(m_gnode.size()) )
      SYS_T::print_fatal("Error: the m_numpts != m_global_node.size()! \n");

    for(unsigned int jj=0; jj<m_gnode.size(); ++jj)
    {
      if(m_gnode[jj]<0) 
        SYS_T::print_fatal("Error: there are negative nodal index!\n");

      minus_nodes.push_back( static_cast<unsigned int>(m_gnode[jj]) );
    }
  }

  VEC_T::sort_unique_resize(minus_nodes);

  // Read vtu file
  SYS_T::file_check( vtufilename );

  vtkXMLUnstructuredGridReader * reader = vtkXMLUnstructuredGridReader::New();
  reader -> SetFileName( vtufilename.c_str() );
  reader -> Update();
  vtkUnstructuredGrid * vtkugrid = reader -> GetOutput();

  const int numpts = static_cast<int>( vtkugrid -> GetNumberOfPoints() );

  vtkPointData * pointdata = vtkugrid->GetPointData();
  vtkDataArray * pd = pointdata->GetScalars("GlobalNodeID");

  std::vector<unsigned int> gnode; gnode.clear();

  for(int ii=0; ii<numpts; ++ii) gnode.push_back( static_cast<unsigned int>(pd->GetComponent(ii,0)) );

  reader->Delete();

  VEC_T::sort_unique_resize( gnode );

  // If the gnode[ii] is not in the minus node, add it to dir_nodes
  for(unsigned int ii=0; ii<gnode.size(); ++ii)
  {
    if( !VEC_T::is_invec(minus_nodes, gnode[ii]) )
      dir_nodes.push_back( gnode[ii] );
  }

  VEC_T::sort_unique_resize(dir_nodes);

  // Read vtp files
  const unsigned int num_file = vtpfileList.size();

  for(unsigned int ii=0; ii<num_file; ++ii)
  {
    SYS_T::file_check( vtpfileList[ii] );

    const std::vector<int> vtp_gnode = VTK_T::read_int_PointData( vtpfileList[ii], "GlobalNodeID");

    for(unsigned int jj=0; jj<vtp_gnode.size(); ++jj)
    {
      if(vtp_gnode[jj]<0) SYS_T::print_fatal("Error: there are negative nodal index! \n");

      dir_nodes.push_back( static_cast<unsigned int>( vtp_gnode[jj] ) );
    }
  }

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  std::cout<<"===> NodalBC_3D_vtu specified by \n";
  for(unsigned int ii=0; ii<num_file; ++ii)
    std::cout<<"     "<<vtpfileList[ii]<<std::endl;
  std::cout<<"     "<<vtufilename<<std::endl;
  std::cout<<"      minus \n";
  for(unsigned int ii=0; ii<num_minus_file; ++ii)
    std::cout<<"     "<<vtpminus[ii]<<std::endl;
  std::cout<<"     is generated. \n";
}

// EOF
