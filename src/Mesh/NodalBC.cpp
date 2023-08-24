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

NodalBC::~NodalBC()
{}

void NodalBC::BC_type_1( const std::vector<std::string> &vtkfileList,
    const int &nFunc )
{
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

void NodalBC::BC_type_3( const std::vector<std::string> &vtkfileList,
    const int &nFunc  )
{
  SYS_T::print_exit_if( vtkfileList.size() != 1, 
      "Error: NodalBC::BC_type_3 the number of vtp files is wrong. \n" );

  SYS_T::file_check( vtpfileList[0] );

  const std::vector<int> gnode = VTK_T::read_int_PointData(vtkfileList[0], "GlobalNodeID");

  SYS_T::print_exit_if( gnode.size() < 1,
      "Error: the numpts is less than 1 in the vtp file! \n");

  SYS_T::print_exit_if( gnode[0]<0,
      "Error: there are negative nodal index! \n");

  dir_nodes.resize(1);
  dir_nodes[0] = static_cast<unsigned int>( gnode[0] );
  num_dir_nodes = 1;
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;
  std::cout<<"-----> Dirichlet node "<<gnode[0]<<
    " from "<<vtpfileList[0]<<std::endl;
}

// EOF
