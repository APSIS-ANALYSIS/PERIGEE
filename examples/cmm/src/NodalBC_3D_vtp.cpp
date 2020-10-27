#include "NodalBC_3D_vtp.hpp"

NodalBC_3D_vtp::NodalBC_3D_vtp( const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_dir_nodes = 0;
  num_per_nodes = 0;

  Create_ID( nFunc );
  
  std::cout<<"===> NodalBC_3D_vtp: No nodal BC is generated. \n";
}


NodalBC_3D_vtp::NodalBC_3D_vtp( const std::string &vtpfileName,
   const int &nFunc )
{
  SYS_T::file_check( vtpfileName );

  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien, gnode, gelem;

  TET_T::read_vtp_grid( vtpfileName, numpts, numcels, pts, ien, gnode, gelem );

  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;
  num_dir_nodes = numpts;

  if( numpts != static_cast<int>(gnode.size()) )
    SYS_T::print_fatal("Error: the numpts != global_node.size()! \n");

  dir_nodes.resize( gnode.size() );
  for(unsigned int ii=0; ii<gnode.size(); ++ii)
  {
    if(gnode[ii]<0) SYS_T::print_fatal("Error: there are negative nodal index! \n");
  
    dir_nodes[ii] = static_cast<unsigned int>( gnode[ii] ); 
  }

  // Generate the ID array
  Create_ID(nFunc);

  std::cout<<"===> NodalBC_3D_vtp specified by "<<vtpfileName<<" is generated. \n";
}


NodalBC_3D_vtp::NodalBC_3D_vtp( const std::vector<std::string> &vtpfileList, 
    const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  const unsigned int num_file = vtpfileList.size();

  for(unsigned int ii=0; ii<num_file; ++ii)
  {
    SYS_T::file_check( vtpfileList[ii] );

    int numpts, numcels;
    std::vector<double> pts;
    std::vector<int> ien, gnode, gelem;

    TET_T::read_vtp_grid( vtpfileList[ii], numpts, numcels, pts, ien, gnode, gelem );
  
    if( numpts != static_cast<int>(gnode.size()) )
      SYS_T::print_fatal("Error: the numpts != global_node.size()! \n");

    for(unsigned int jj=0; jj<gnode.size(); ++jj)
    {
      if(gnode[jj]<0) SYS_T::print_fatal("Error: there are negative nodal index! \n");

      dir_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
    }
  }
  
  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  std::cout<<"===> NodalBC_3D_vtp specified by \n";
  for(unsigned int ii=0; ii<num_file; ++ii)
    std::cout<<"     "<<vtpfileList[ii]<<std::endl;
  std::cout<<"     is generated. \n";
}


NodalBC_3D_vtp::NodalBC_3D_vtp( const std::vector<std::string> &vtpfileList,
    const int &nFunc, const std::vector<int> &master_idx )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_dir_nodes = 0;

  const unsigned int num_file = vtpfileList.size();

  if( num_file != master_idx.size() )
    SYS_T::print_fatal("Error: the file size does not match the master idx size.\n");

  for(unsigned int ii=0; ii<num_file; ++ii)
  {
    SYS_T::file_check( vtpfileList[ii] );

    int numpts, numcels;
    std::vector<double> pts;
    std::vector<int> ien, gnode, gelem;

    TET_T::read_vtp_grid( vtpfileList[ii], numpts, numcels, pts, ien, gnode, gelem );
  
    if( numpts != static_cast<int>(gnode.size()) )
      SYS_T::print_fatal("Error: the numpts != global_node.size()! \n");

    for(unsigned int jj=0; jj<gnode.size(); ++jj)
    {
      if(gnode[jj]<0) SYS_T::print_fatal("Error: there are negative nodal index! \n");

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

  std::cout<<"===> NodalBC_3D_vtp specified by \n";
  for(unsigned int ii=0; ii<num_file; ++ii)
    std::cout<<"     "<<vtpfileList[ii]<<" follows "<<master_idx[ii]<<std::endl;
  std::cout<<"     is generated. \n";
}


NodalBC_3D_vtp::NodalBC_3D_vtp( const std::vector<std::string> &vtpfileList,
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
      BC_type_1( vtpfileList, nFunc );
      break;
    case 2:
      BC_type_2( vtpfileList, nFunc );
      break;
    case 3:
      BC_type_3( vtpfileList, nFunc );
      break;
    default:
      std::cerr<<"Error: NodalBC_3D_vtp with bc type = "<<type<<" is not implemented. \n";
      exit(EXIT_FAILURE);
  }

  Create_ID( nFunc );

  log_time = clock() - log_time;
  std::cout<<"===> NodalBC_3D_vtp, type = "<<type<<" is generated. ";
  std::cout<<"Time taken: "<<((float) log_time)/CLOCKS_PER_SEC<<" seconds."<<std::endl; 
}


NodalBC_3D_vtp::~NodalBC_3D_vtp()
{}


void NodalBC_3D_vtp::BC_type_1( const std::vector<std::string> &vtpfileList,
    const int &nFunc  )
{
  const unsigned int num_file = vtpfileList.size();

  for(unsigned int ii=0; ii<num_file; ++ii)
  {
    SYS_T::file_check( vtpfileList[ii] );

    int numpts, numcels;
    std::vector<double> pts;
    std::vector<int> ien, gnode, gelem;

    TET_T::read_vtp_grid( vtpfileList[ii], numpts, numcels, pts, ien, gnode, gelem );

    if( numpts != static_cast<int>(gnode.size()) )
      SYS_T::print_fatal("Error: the numpts != global_node.size()! \n");

    for(unsigned int jj=1; jj<gnode.size(); ++jj)
    {
      if(gnode[jj]<0) SYS_T::print_fatal("Error: there are negative nodal index! \n");

      per_slave_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
      per_master_nodes.push_back( static_cast<unsigned int>(gnode[ 0 ]) );
    }
  }

  num_per_nodes = per_slave_nodes.size();
  std::cout<<"-----> Master slave relations: \n";
  for(unsigned int ii=0; ii<num_file; ++ii)
    std::cout<<"     "<<vtpfileList[ii]<<" follows 0th node in the file "<<std::endl;
}


void NodalBC_3D_vtp::BC_type_2( const std::vector<std::string> &vtpfileList,
    const int &nFunc  )
{
  const unsigned int num_file = vtpfileList.size();

  if( num_file != 2 ) 
    SYS_T::print_fatal("Error: NodalBC_3D_vtp::BC_type_2 the number of vtp files is wrong. \n");

  SYS_T::file_check( vtpfileList[0] );

  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien, gnode, gelem;

  TET_T::read_vtp_grid( vtpfileList[0], numpts, numcels, pts, ien, gnode, gelem );

  if( numpts != static_cast<int>(gnode.size()) )
    SYS_T::print_fatal("Error: the numpts != global_node.size()! \n");

  dir_nodes.resize( gnode.size() );
  for(unsigned int ii=0; ii<gnode.size(); ++ii)
  {
    if(gnode[ii]<0) SYS_T::print_fatal("Error: there are negative nodal index! \n");

    dir_nodes[ii] = static_cast<unsigned int>( gnode[ii] );
  }

  SYS_T::file_check( vtpfileList[1] );
  TET_T::read_vtp_grid( vtpfileList[1], numpts, numcels, pts, ien, gnode, gelem );

  if( numpts != static_cast<int>(gnode.size()) )
    SYS_T::print_fatal("Error: the numpts != global_node.size()! \n");

  for(unsigned int jj=1; jj<gnode.size(); ++jj)
  {
    if(gnode[jj]<0) SYS_T::print_fatal("Error: there are negative nodal index! \n");

    per_slave_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
    per_master_nodes.push_back( static_cast<unsigned int>( gnode[ 0 ]) );
  }

  num_dir_nodes = dir_nodes.size();
  num_per_nodes = per_slave_nodes.size();

  std::cout<<"-----> Dirichlet nodes from "<<vtpfileList[0]<<std::endl;
  std::cout<<"       Master-slave from "<<vtpfileList[1]<<" with 0th node master."<<std::endl;
}


void NodalBC_3D_vtp::BC_type_3( const std::vector<std::string> &vtpfileList,
    const int &nFunc  )
{
  const unsigned int num_file = vtpfileList.size();

  SYS_T::print_fatal_if( num_file != 1, 
      "Error: NodalBC_3D_vtp::BC_type_3 the number of vtp files is wrong. \n" );

  SYS_T::file_check( vtpfileList[0] );

  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien, gnode, gelem;

  TET_T::read_vtp_grid( vtpfileList[0], numpts, numcels, pts, ien, gnode, gelem );

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
