#include "NodalBC.hpp"

NodalBC::NodalBC( const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_dir_nodes = 0;
  num_per_nodes = 0;

  Create_ID( nFunc );
  
  std::cout<<"===> NodalBC: No nodal BC is generated. \n";
}

NodalBC::NodalBC( const std::vector<std::string> &vtkfileList, 
    const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  for( const auto &vtkfile : vtkfileList )
  {
    SYS_T::file_check( vtkfile );

    const auto gnode = VTK_T::read_int_PointData(vtkfile, "GlobalNodeID");
  
    for(unsigned int jj=0; jj<gnode.size(); ++jj)
    {
      SYS_T::print_exit_if(gnode[jj]<0, "Error: there are negative nodal index! \n");

      dir_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
    }
  }
  
  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  std::cout<<"===> NodalBC specified by \n";
  for( const auto &vtkfile : vtkfileList )
    std::cout<<"     "<<vtkfile<<"\n";
  std::cout<<"     is generated."<<std::endl;
}


NodalBC::NodalBC( const std::vector<std::string> &vtkfileList,
    const int &nFunc, const std::vector<int> &master_idx )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_dir_nodes = 0;

  const unsigned int num_file = vtkfileList.size();

  SYS_T::print_exit_if(vtkfileList.size() != master_idx.size(),
    "Error: the file size does not match the master idx size.\n");

  for (const auto &vtkfile : vtkfileList)
  {
    SYS_T::file_check( vtkfile );

    const auto gnode = VTK_T::read_int_PointData(vtkfile, "GlobalNodeID");

    for(unsigned int jj=0; jj<gnode.size(); ++jj)
    {
      SYS_T::print_exit_if(gnode[jj]<0, "Error: there are negative nodal index! \n");

      if( static_cast<int>(jj) != master_idx[ii] )
      {
        per_slave_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
        per_master_nodes.push_back( static_cast<unsigned int>(gnode[ master_idx[ii] ]) );
      }
      else
      {
        std::cout<<"jj = "<<jj<<'\t'<<gnode[jj]<<" matched master index "<<master_idx[ii]<<std::endl;
      }
    }
  }

  num_per_nodes = per_slave_nodes.size();

  Create_ID( nFunc );

  std::cout<<"===> NodalBC specified by \n";
  for(const auto &vtkfile : vtkfileList)
    std::cout<<"     "<<vtkfile<<" follows "<<master_idx[ii]<<std::endl;
  std::cout<<"     is generated. \n";
}


NodalBC::NodalBC( const std::vector<std::string> &vtkfileList,
    const int &nFunc, const int &type )
{
  clock_t log_time = clock();
 
  // Clean allocation first 
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();

  num_dir_nodes = 0;
  num_per_nodes = 0;

  switch( type )
  {
    case 1:
      BC_type_1( vtkfileList, nFunc );
      break;
    case 2:
      BC_type_2( vtkfileList, nFunc );
      break;
    case 3:
      BC_type_3( vtkfileList, nFunc );
      break;
    default:
      std::cerr<<"Error: NodalBC with bc type = "<<type<<" is not implemented. \n";
      exit(EXIT_FAILURE);
  }

  Create_ID( nFunc );

  log_time = clock() - log_time;
  std::cout<<"===> NodalBC, type = "<<type<<" is generated. ";
  std::cout<<"Time taken: "<<((float) log_time)/CLOCKS_PER_SEC<<" seconds."<<std::endl; 
}

NodalBC::NodalBC( const std::string &vtufile, 
        const std::vector<std::string> &vtpfileList, const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  // Read vtp files
  for( const auto &vtpfile : vtpfileList )
  {
    SYS_T::file_check( vtpfile );

    const std::vector<int> gnode = VTK_T::read_int_PointData( vtpfile , "GlobalNodeID");

    for(unsigned int jj=0; jj<gnode.size(); ++jj)
    {
      SYS_T::print_exit_if(gnode[jj]<0, "Error: there are negative nodal index! \n");

      dir_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
    }
  }

  // Read vtu file
  SYS_T::file_check( vtufile );

  // vtkXMLUnstructuredGridReader * reader = vtkXMLUnstructuredGridReader::New();
  // reader -> SetFileName( vtufile.c_str() );
  // reader -> Update();
  // vtkUnstructuredGrid * vtkugrid = reader -> GetOutput();

  // const int numpts = static_cast<int>( vtkugrid -> GetNumberOfPoints() );

  // vtkPointData * pointdata = vtkugrid->GetPointData();
  // vtkDataArray * pd = pointdata->GetScalars("GlobalNodeID");

  // std::vector<unsigned int> gnode; gnode.clear();

  // for(int ii=0; ii<numpts; ++ii) gnode.push_back( static_cast<unsigned int>(pd->GetComponent(ii,0)) );

  // reader->Delete();

  std::vector<int> gnode = VTK_T::read_int_PointData( vtufile, "GlobalNodeID" ); // Line 159 ~ 173

  VEC_T::sort_unique_resize( gnode );

  // VEC_T::insert_end(dir_nodes, gnode );

  VEC_T::insert_end(dir_nodes, std::vector<unsigned int> (gnode.begin(), gnode.end()) ); // Line 179

  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  std::cout<<"===> NodalBC_3D_vtu specified by \n";
  for(const auto &vtpfile : vtpfileList)
    std::cout<<"     "<<vtpfile<<std::endl;
  std::cout<<"     "<<vtufile<<"     is generated. \n";
}

NodalBC::~NodalBC()
{}

void NodalBC::BC_type_1( const std::vector<std::string> &vtkfileList,
    const int &nFunc )
{
  const unsigned int num_file = vtpfileList.size();

  for( const auto &vtkfile : vtkfileList )
  {
    SYS_T::file_check( vtkfile );

    const std::vector<int> gnode = VTK_T::read_int_PointData(vtkfile, "GlobalNodeID");

    for(unsigned int jj=1; jj<gnode.size(); ++jj)
    {
      SYS_T::print_exit_if(gnode[jj]<0, "Error: there are negative nodal index! \n");

      per_slave_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
      per_master_nodes.push_back( static_cast<unsigned int>(gnode[ 0 ]) );
    }
  }

  num_per_nodes = per_slave_nodes.size();
  std::cout<<"-----> Master slave relations: \n";
  for( const auto &vtkfile : vtkfileList )
    std::cout<<"     "<<vtkfile<<" follows 0th node in the file "<<std::endl;
}


void NodalBC::BC_type_2( const std::vector<std::string> &vtkfileList,
    const int &nFunc )
{
  SYS_T::print_exit_if(vtkfileList.size() != 2,
    "Error: NodalBC::BC_type_2 the number of vtp files is wrong. \n")

  SYS_T::file_check( vtpfileList[0] );

  std::vector<int> gnode = VTK_T::read_int_PointData(vtpfileList[0], "GlobalNodeID");

  dir_nodes.resize( gnode.size() );
  for(unsigned int ii=0; ii<gnode.size(); ++ii)
  {
    SYS_T::print_exit_if(gnode[ii]<0, "Error: there are negative nodal index! \n");

    dir_nodes[ii] = static_cast<unsigned int>( gnode[ii] );
  }

  SYS_T::file_check( vtpfileList[1] );

  gnode = VTK_T::read_int_PointData(vtpfileList[1], "GlobalNodeID");

  for(unsigned int jj=1; jj<gnode.size(); ++jj)
  {
    SYS_T::print_exit_if(gnode[jj]<0, "Error: there are negative nodal index! \n");

    per_slave_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
    per_master_nodes.push_back( static_cast<unsigned int>( gnode[ 0 ]) );
  }

  num_dir_nodes = dir_nodes.size();
  num_per_nodes = per_slave_nodes.size();

  std::cout<<"-----> Dirichlet nodes from "<<vtpfileList[0]<<std::endl;
  std::cout<<"       Master-slave from "<<vtpfileList[1]<<" with 0th node master."<<std::endl;
}


void NodalBC::BC_type_3( const std::vector<std::string> &vtpfileList,
    const int &nFunc  )
{
  const unsigned int num_file = vtpfileList.size();

  SYS_T::print_fatal_if( num_file != 1, 
      "Error: NodalBC::BC_type_3 the number of vtp files is wrong. \n" );

  SYS_T::file_check( vtpfileList[0] );

  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien;

  VTK_T::read_vtp_grid( vtpfileList[0], numpts, numcels, pts, ien );

  const std::vector<int> gnode = VTK_T::read_int_PointData(vtpfileList[0], "GlobalNodeID");

  SYS_T::print_fatal_if( numpts != static_cast<int>(gnode.size()),
      "Error: the numpts != global_node.size()! \n");

  SYS_T::print_fatal_if( numpts < 1,
      "Error: the numpts is less than 1 in the vtp file! \n");

  SYS_T::print_fatal_if( gnode[0]<0,
      "Error: there are negative nodal index! \n");

  dir_nodes.resize(1);
  dir_nodes[0] = static_cast<unsigned int>( gnode[0] );
  num_dir_nodes = dir_nodes.size();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;
  std::cout<<"-----> Dirichlet nodes "<<gnode[0]<<
    " from "<<vtpfileList[0]<<std::endl;
}

// EOF