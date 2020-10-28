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


NodalBC_3D_vtu::NodalBC_3D_vtu( const std::string &inflow_vtu_file,
    const std::string &wall_vtu_file,
    const std::vector<std::string> &outflow_vtu_files,
    const int &nFunc )
{
  // TO BE IMPLEMENTED.
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
    std::vector<int> ien, gnode, gelem;

    TET_T::read_vtu_grid( vtufileList[ii], numpts, numcels, pts, ien, gnode, gelem );

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


NodalBC_3D_vtu::~NodalBC_3D_vtu()
{}

// EOF
