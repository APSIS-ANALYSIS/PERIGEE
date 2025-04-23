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

  for( int ii{0}; ii < VEC_T::get_size(vtkfileList); ++ii )
  {
    const auto vtkfile = vtkfileList[ii];
    SYS_T::file_check( vtkfile );

    const auto gnode = VTK_T::read_int_PointData(vtkfile, "GlobalNodeID");
  
    for(unsigned int jj=0; jj<gnode.size(); ++jj)
    {
      SYS_T::print_fatal_if( gnode[jj]<0 || gnode[jj]>=nFunc, "Error: the nodal index %d is not in the range [0, %d)! \n", gnode[jj], nFunc);

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
  const std::vector<std::string> &slafileList,
  const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  for( int ii{0}; ii < VEC_T::get_size(vtkfileList); ++ii )
  {
    const auto vtkfile = vtkfileList[ii];
    SYS_T::file_check( vtkfile );

    const auto gnode = VTK_T::read_int_PointData(vtkfile, "GlobalNodeID");

    for(unsigned int jj=0; jj<gnode.size(); ++jj)
    {
      SYS_T::print_fatal_if( gnode[jj]<0 || gnode[jj]>=nFunc, "Error: the nodal index %d is not in the range [0, %d)! \n", gnode[jj], nFunc);

      dir_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
    }
  }

  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  for( int ii{0}; ii < VEC_T::get_size(slafileList); ++ii )
  {
    const auto slafile = slafileList[ii];
    SYS_T::file_check( slafile );

    const auto sla_gnode = VTK_T::read_int_PointData(slafile, "GlobalNodeID");
    const auto mas_gnode = VTK_T::read_int_PointData(slafile, "MasterNodeID");

    for(unsigned int jj=0; jj<sla_gnode.size(); ++jj)
    {
      SYS_T::print_fatal_if( sla_gnode[jj]<0 || sla_gnode[jj]>=nFunc, "Error: the slave nodal index %d is not in the range [0, %d)! \n", sla_gnode[jj], nFunc);
      SYS_T::print_fatal_if( mas_gnode[jj]<0 || mas_gnode[jj]>=nFunc, "Error: the master nodal index %d is not in the range [0, %d)! \n", mas_gnode[jj], nFunc);

      per_slave_nodes.push_back( static_cast<unsigned int>( sla_gnode[jj]) );
      per_master_nodes.push_back( static_cast<unsigned int>( mas_gnode[jj]) );
    }
  }

  num_per_nodes = per_slave_nodes.size();

  Create_ID( nFunc );

  std::cout<<"===> NodalBC specified by \n";
  for( const auto &vtkfile : vtkfileList )
    std::cout<<"     "<<vtkfile<<"\n";
  std::cout<<"     is generated."<<std::endl;

  std::cout<<"===> Periodic NodalBC specified by \n";
  for( const auto &slavtkfile : slafileList )
    std::cout<<"     "<<slavtkfile<<"\n";
  std::cout<<"     is generated."<<std::endl;
}

NodalBC::NodalBC( const std::vector<std::string> &vtkfileList,
    const std::string &rotated_file, 
    const std::string &fixed_file,    
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
      SYS_T::print_fatal_if( gnode[jj]<0 || gnode[jj]>=nFunc, "Error: the nodal index %d is not in the range [0, %d)! \n", gnode[jj], nFunc);

      dir_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
    }
  }
  
  // Prepare the numbers that need to be shifted
  SYS_T::file_check( rotated_file ); 
  SYS_T::file_check( fixed_file );  

  const int fixed_nFunc = VTK_T::read_num_pt(fixed_file);

  auto rot_gnode = VTK_T::read_int_PointData(rotated_file, "GlobalNodeID");

  for(int &nodeid : rot_gnode)
    nodeid += fixed_nFunc;

  for(unsigned int jj=0; jj<rot_gnode.size(); ++jj)
  {
    SYS_T::print_fatal_if( rot_gnode[jj]<0 || rot_gnode[jj]>=nFunc, "Error: the nodal index %d is not in the range [0, %d)! \n", rot_gnode[jj], nFunc);

    dir_nodes.push_back( static_cast<unsigned int>( rot_gnode[jj]) );
  }

  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  std::cout<<"===> NodalBC specified by \n";
  for( const auto &vtkfile : vtkfileList )
    std::cout<<"     "<<vtkfile<<"\n";
  std::cout<<"     "<<rotated_file<<"\n";
  std::cout<<"     is generated."<<std::endl;
}

NodalBC::NodalBC( const std::vector<std::string> &vtkfileList,
    const std::string &rotated_file, 
    const std::string &rotated_wall_file, 
    const std::string &fixed_file,    
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
      SYS_T::print_fatal_if( gnode[jj]<0 || gnode[jj]>=nFunc, "Error: the nodal index %d is not in the range [0, %d)! \n", gnode[jj], nFunc);

      dir_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
    }
  }
  
  // Prepare the numbers that need to be shifted
  // SYS_T::file_check( rotated_file );
  SYS_T::file_check( rotated_wall_file ); 
  SYS_T::file_check( fixed_file );  

  const int fixed_nFunc = VTK_T::read_num_pt(fixed_file);

  auto rot_gnode = VTK_T::read_int_PointData(rotated_file, "GlobalNodeID");

  for(int &nodeid : rot_gnode)
    nodeid += fixed_nFunc;

  auto rot_wall_gnode = VTK_T::read_int_PointData(rotated_wall_file, "GlobalNodeID");

  for(int &nodeid : rot_wall_gnode)
    nodeid += fixed_nFunc;

  for(unsigned int jj=0; jj<rot_gnode.size(); ++jj)
  {
    SYS_T::print_fatal_if( rot_gnode[jj]<0 || rot_gnode[jj]>=nFunc, "Error: the nodal index %d is not in the range [0, %d)! \n", rot_gnode[jj], nFunc);

    dir_nodes.push_back( static_cast<unsigned int>( rot_gnode[jj]) );
  }

  for(unsigned int jj=0; jj<rot_wall_gnode.size(); ++jj)
  {
    SYS_T::print_fatal_if( rot_wall_gnode[jj]<0 || rot_wall_gnode[jj]>=nFunc, "Error: the nodal index %d is not in the range [0, %d)! \n", rot_wall_gnode[jj], nFunc);

    dir_nodes.push_back( static_cast<unsigned int>( rot_wall_gnode[jj]) );
  }

  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  std::cout<<"===> NodalBC specified by \n";
  for( const auto &vtkfile : vtkfileList )
    std::cout<<"     "<<vtkfile<<"\n";
  std::cout<<"     "<<rotated_file<<"\n";
  std::cout<<"     "<<rotated_wall_file<<"\n";
  std::cout<<"     is generated."<<std::endl;
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

void NodalBC::BC_type_1( const std::vector<std::string> &vtkfileList,
    const int &nFunc )
{
  for( const auto &vtkfile : vtkfileList )
  {
    SYS_T::file_check( vtkfile );

    const std::vector<int> gnode = VTK_T::read_int_PointData(vtkfile, "GlobalNodeID");

    for(unsigned int jj=1; jj<gnode.size(); ++jj)
    {
      SYS_T::print_fatal_if( gnode[jj]<0 || gnode[jj]>=nFunc, "Error: the nodal index %d is not in the range [0, %d)! \n", gnode[jj], nFunc);

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
  SYS_T::print_fatal_if(vtkfileList.size() != 2, "Error: NodalBC::BC_type_2 the number of vtk files is wrong. \n");

  SYS_T::file_check( vtkfileList[0] );

  std::vector<int> gnode = VTK_T::read_int_PointData(vtkfileList[0], "GlobalNodeID");

  dir_nodes.resize( gnode.size() );
  for(unsigned int ii=0; ii<gnode.size(); ++ii)
  {
    SYS_T::print_fatal_if( gnode[ii]<0 || gnode[ii]>=nFunc, "Error: the nodal index %d is not in the range [0, %d)! \n", gnode[ii], nFunc);

    dir_nodes[ii] = static_cast<unsigned int>( gnode[ii] );
  }

  SYS_T::file_check( vtkfileList[1] );

  gnode = VTK_T::read_int_PointData(vtkfileList[1], "GlobalNodeID");

  for(unsigned int jj=1; jj<gnode.size(); ++jj)
  {
    SYS_T::print_fatal_if( gnode[jj]<0 || gnode[jj]>=nFunc, "Error: the nodal index %d is not in the range [0, %d)! \n", gnode[jj], nFunc);

    per_slave_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
    per_master_nodes.push_back( static_cast<unsigned int>( gnode[ 0 ]) );
  }

  num_dir_nodes = dir_nodes.size();
  num_per_nodes = per_slave_nodes.size();

  std::cout<<"-----> Dirichlet nodes from "<<vtkfileList[0]<<std::endl;
  std::cout<<"       Master-slave from "<<vtkfileList[1]<<" with 0th node master."<<std::endl;
}

void NodalBC::BC_type_3( const std::vector<std::string> &vtkfileList,
    const int &nFunc  )
{
  SYS_T::print_fatal_if( vtkfileList.size() != 1, "Error: NodalBC::BC_type_3 the number of vtp files is wrong. \n" );

  SYS_T::file_check( vtkfileList[0] );

  const std::vector<int> gnode = VTK_T::read_int_PointData(vtkfileList[0], "GlobalNodeID");

  SYS_T::print_fatal_if( gnode.size() < 1, "Error: the numpts is less than 1 in the vtp file! \n");

  SYS_T::print_fatal_if( gnode[0]<0 || gnode[0]>=nFunc, "Error: the nodal index %d is not in the range [0, %d)! \n", gnode[0], nFunc);

  dir_nodes.resize(1);
  dir_nodes[0] = static_cast<unsigned int>( gnode[0] );
  num_dir_nodes = 1;
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;
  std::cout<<"-----> Dirichlet node "<<gnode[0]<<
    " from "<<vtkfileList[0]<<std::endl;
}

// EOF
