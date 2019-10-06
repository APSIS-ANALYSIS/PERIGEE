#include "NodalBC_3D_Linearelastic.hpp"

NodalBC_3D_Linearelastic::NodalBC_3D_Linearelastic( const int &nFunc, 
    const int &nFunc_x, const int &nFunc_y, const int &nFunc_z, 
    const int &bc_type)
{
  clock_t log_time = clock();

  // Clean the nodal BC data structure
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();

  num_dir_nodes = 0;
  num_per_nodes = 0;

  switch(bc_type)
  {
    case 1:
      BC_type_1(nFunc_x, nFunc_y, nFunc_z);
      break;
    case 2:
      BC_type_2();
      break;
    case 3:
      BC_type_3(nFunc_x, nFunc_y, nFunc_z);
      break;
    case 4:
      BC_type_4(nFunc_x, nFunc_y, nFunc_z);
      break;
    case 5:
      BC_type_5(nFunc_x, nFunc_y, nFunc_z);
      break;
    case 6:
      BC_type_6(nFunc_x, nFunc_y, nFunc_z);
      break;
    case 7:
      BC_type_7(nFunc_x, nFunc_y, nFunc_z);
      break;
    case 8:
      BC_type_8(nFunc_x, nFunc_y, nFunc_z);
      break;
    case 9:
      BC_type_9(nFunc_x, nFunc_y, nFunc_z);
      break;
    case 10:
      BC_type_10(nFunc_x, nFunc_y, nFunc_z);
      break;
    default:
      std::cerr<<"Error: NodalBC_3D_Linearelastic with bc_type = "
        <<bc_type<<" is not implemented. \n";
      exit(EXIT_FAILURE);
  }

  Create_ID(nFunc);

  log_time = clock() - log_time;

  std::cout<<"===> NodalBC_3D_Linearelastic, type = "<<bc_type<<" is generated. ";
  std::cout<<"Time taken: "<<((float) log_time)/CLOCKS_PER_SEC<<" seconds. \n"<<std::endl;
}



NodalBC_3D_Linearelastic::~NodalBC_3D_Linearelastic()
{}



void NodalBC_3D_Linearelastic::BC_type_1( const int &nFunc_x, 
    const int &nFunc_y, const int &nFunc_z )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  Generate_BCNode_3D_A(nFunc_x, nFunc_y, nFunc_z,
      nodes_front, nodes_back, nodes_left,
      nodes_right, nodes_top, nodes_bottom );
  dir_nodes.clear();
  dir_nodes.insert(dir_nodes.end(), nodes_front.begin(), nodes_front.end());
  dir_nodes.insert(dir_nodes.end(), nodes_back.begin(), nodes_back.end());
  dir_nodes.insert(dir_nodes.end(), nodes_left.begin(), nodes_left.end());
  dir_nodes.insert(dir_nodes.end(), nodes_right.begin(), nodes_right.end());
  dir_nodes.insert(dir_nodes.end(), nodes_top.begin(), nodes_top.end());
  dir_nodes.insert(dir_nodes.end(), nodes_bottom.begin(), nodes_bottom.end());

  num_dir_nodes = dir_nodes.size();
  std::cout<<"-----> Dirichlet BC on all 6 faces. \n";
}


void NodalBC_3D_Linearelastic::BC_type_2()
{
  std::cout<<"-----> No nodal BC. \n";
}



void NodalBC_3D_Linearelastic::BC_type_3( const int &nFunc_x, 
    const int &nFunc_y, const int &nFunc_z )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  Generate_BCNode_3D_A(nFunc_x, nFunc_y, nFunc_z,
      nodes_front, nodes_back, nodes_left,
      nodes_right, nodes_top, nodes_bottom );
  dir_nodes.clear();
  dir_nodes.insert(dir_nodes.end(), nodes_front.begin(), nodes_front.end());
  dir_nodes.insert(dir_nodes.end(), nodes_back.begin(), nodes_back.end());
  dir_nodes.insert(dir_nodes.end(), nodes_left.begin(), nodes_left.end());
  dir_nodes.insert(dir_nodes.end(), nodes_right.begin(), nodes_right.end());
  dir_nodes.insert(dir_nodes.end(), nodes_bottom.begin(), nodes_bottom.end());

  num_dir_nodes = dir_nodes.size();

  std::cout<<"-----> Dirchlet nodal BC on LEFT, RIGHT, FRONT, BACK, BOTTOM. \n";
}



void NodalBC_3D_Linearelastic::BC_type_4( const int &nFunc_x, 
    const int &nFunc_y, const int &nFunc_z )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  Generate_BCNode_3D_A(nFunc_x, nFunc_y, nFunc_z,
      nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom );
  
  dir_nodes.clear();
  dir_nodes.insert(dir_nodes.end(), nodes_back.begin(), nodes_back.end());

  num_dir_nodes = dir_nodes.size();

  std::cout<<"-----> Dirchlet nodal BC on BACK. \n";
}


void NodalBC_3D_Linearelastic::BC_type_5( const int &nFunc_x, 
    const int &nFunc_y, const int &nFunc_z )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  Generate_BCNode_3D_A(nFunc_x, nFunc_y, nFunc_z,
      nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom );
  
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_back.begin(), nodes_back.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_front.begin(), nodes_front.end());

  num_per_nodes = per_slave_nodes.size();

  std::cout<<"-----> C0 periodic pair for Front-Back nodes. \n";
}



void NodalBC_3D_Linearelastic::BC_type_6( const int &nFunc_x, 
    const int &nFunc_y, const int &nFunc_z )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  std::vector<int> e01,e02, e13, e23, e45, e46, e57, e67, e15, e37, e04, e26, cor;

  Generate_BCNode_3D_B(nFunc_x, nFunc_y, nFunc_z,
      nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom,
     e01, e02, e13, e23, e45, e46, e57, e67, e15, e37, e04, e26, cor );
  
  dir_nodes.clear();
  dir_nodes.insert(dir_nodes.end(), nodes_bottom.begin(), nodes_bottom.end());
  dir_nodes.insert(dir_nodes.end(), e01.begin(), e01.end());
  dir_nodes.insert(dir_nodes.end(), e02.begin(), e02.end());
  dir_nodes.insert(dir_nodes.end(), e13.begin(), e13.end());
  dir_nodes.insert(dir_nodes.end(), e23.begin(), e23.end());
  dir_nodes.push_back(cor[0]);
  dir_nodes.push_back(cor[1]);
  dir_nodes.push_back(cor[2]);
  dir_nodes.push_back(cor[3]);

  num_dir_nodes = dir_nodes.size();

  std::cout<<"-----> Dirchlet nodal BC on BOTTOM. \n";
}


void NodalBC_3D_Linearelastic::BC_type_7( const int &nFunc_x, 
    const int &nFunc_y, const int &nFunc_z )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  std::vector<int> e01,e02, e13, e23, e45, e46, e57, e67, e15, e37, e04, e26, cor;

  Generate_BCNode_3D_B(nFunc_x, nFunc_y, nFunc_z,
      nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom,
     e01, e02, e13, e23, e45, e46, e57, e67, e15, e37, e04, e26, cor );
  
  dir_nodes.clear();
  dir_nodes.insert(dir_nodes.end(), nodes_bottom.begin(), nodes_bottom.end());
  dir_nodes.insert(dir_nodes.end(), nodes_top.begin(), nodes_top.end());

  dir_nodes.insert(dir_nodes.end(), e01.begin(), e01.end());
  dir_nodes.insert(dir_nodes.end(), e02.begin(), e02.end());
  dir_nodes.insert(dir_nodes.end(), e13.begin(), e13.end());
  dir_nodes.insert(dir_nodes.end(), e23.begin(), e23.end());
  
  dir_nodes.insert(dir_nodes.end(), e45.begin(), e45.end());
  dir_nodes.insert(dir_nodes.end(), e46.begin(), e46.end());
  dir_nodes.insert(dir_nodes.end(), e67.begin(), e67.end());
  dir_nodes.insert(dir_nodes.end(), e57.begin(), e57.end());
  
  dir_nodes.insert(dir_nodes.end(), cor.begin(), cor.end());
  
  num_dir_nodes = dir_nodes.size();

  std::cout<<"-----> Dirchlet nodal BC on TOP & BOTTOM. \n";
}


void NodalBC_3D_Linearelastic::BC_type_8( const int &nFunc_x, 
    const int &nFunc_y, const int &nFunc_z )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  std::vector<int> e01,e02, e13, e23, e45, e46, e57, e67, e15, e37, e04, e26, cor;

  Generate_BCNode_3D_B(nFunc_x, nFunc_y, nFunc_z,
      nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom,
     e01, e02, e13, e23, e45, e46, e57, e67, e15, e37, e04, e26, cor );
  
  dir_nodes.clear();
  dir_nodes.insert(dir_nodes.end(), nodes_top.begin(), nodes_top.end());
  dir_nodes.insert(dir_nodes.end(), nodes_back.begin(), nodes_back.end());

  dir_nodes.insert(dir_nodes.end(), e02.begin(), e02.end());
  
  dir_nodes.insert(dir_nodes.end(), e45.begin(), e45.end());
  dir_nodes.insert(dir_nodes.end(), e46.begin(), e46.end());
  dir_nodes.insert(dir_nodes.end(), e67.begin(), e67.end());
  dir_nodes.insert(dir_nodes.end(), e57.begin(), e57.end());
  
  dir_nodes.insert(dir_nodes.end(), e04.begin(), e04.end());
  dir_nodes.insert(dir_nodes.end(), e26.begin(), e26.end());
  
  dir_nodes.push_back(cor[0]);
  dir_nodes.push_back(cor[2]);
  dir_nodes.push_back(cor[4]);
  dir_nodes.push_back(cor[5]);
  dir_nodes.push_back(cor[6]);
  dir_nodes.push_back(cor[7]);
  
  num_dir_nodes = dir_nodes.size();

  std::cout<<"-----> Dirchlet nodal BC on TOP & BACK. \n";
}



void NodalBC_3D_Linearelastic::BC_type_9( const int &nFunc_x, 
    const int &nFunc_y, const int &nFunc_z )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  std::vector<int> e01,e02, e13, e23, e45, e46, e57, e67, e15, e37, e04, e26, cor;

  Generate_BCNode_3D_B(nFunc_x, nFunc_y, nFunc_z,
      nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom,
     e01, e02, e13, e23, e45, e46, e57, e67, e15, e37, e04, e26, cor );
  
  dir_nodes.clear();
  dir_nodes.insert(dir_nodes.end(), nodes_top.begin(), nodes_top.end());
  dir_nodes.insert(dir_nodes.end(), nodes_left.begin(), nodes_left.end());

  dir_nodes.insert(dir_nodes.end(), e01.begin(), e01.end());
  
  dir_nodes.insert(dir_nodes.end(), e45.begin(), e45.end());
  dir_nodes.insert(dir_nodes.end(), e46.begin(), e46.end());
  dir_nodes.insert(dir_nodes.end(), e67.begin(), e67.end());
  dir_nodes.insert(dir_nodes.end(), e57.begin(), e57.end());
  
  dir_nodes.insert(dir_nodes.end(), e04.begin(), e04.end());
  dir_nodes.insert(dir_nodes.end(), e15.begin(), e15.end());
  
  dir_nodes.push_back(cor[0]);
  dir_nodes.push_back(cor[1]);
  dir_nodes.push_back(cor[4]);
  dir_nodes.push_back(cor[5]);
  dir_nodes.push_back(cor[6]);
  dir_nodes.push_back(cor[7]);

  num_dir_nodes = dir_nodes.size();

  std::cout<<"-----> Dirchlet nodal BC on TOP & LEFT. \n";
}


void NodalBC_3D_Linearelastic::BC_type_10( const int &nFunc_x,
    const int &nFunc_y, const int &nFunc_z )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  std::vector<int> e01,e02, e13, e23, e45, e46, e57, e67, e15, e37, e04, e26, cor;

  Generate_BCNode_3D_B(nFunc_x, nFunc_y, nFunc_z,
      nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom,
      e01, e02, e13, e23, e45, e46, e57, e67, e15, e37, e04, e26, cor );

  dir_nodes.clear();
  dir_nodes.insert(dir_nodes.end(), nodes_bottom.begin(), nodes_bottom.end());
  dir_nodes.insert(dir_nodes.end(), e01.begin(), e01.end());
  dir_nodes.insert(dir_nodes.end(), e02.begin(), e02.end());
  dir_nodes.insert(dir_nodes.end(), e13.begin(), e13.end());
  dir_nodes.insert(dir_nodes.end(), e23.begin(), e23.end());
  dir_nodes.push_back( cor[0] );
  dir_nodes.push_back( cor[1] );
  dir_nodes.push_back( cor[2] );
  dir_nodes.push_back( cor[3] );

  num_dir_nodes = dir_nodes.size();

  for( unsigned int ii=0; ii<nodes_top.size(); ++ii )
  {
    per_slave_nodes.push_back( nodes_top[ii] );
    per_master_nodes.push_back( cor[4] );
  }

  for(unsigned int ii=0; ii<e45.size(); ++ii)
  {
    per_slave_nodes.push_back( e45[ii] );
    per_master_nodes.push_back( cor[4] );
  }

  for(unsigned int ii=0; ii<e46.size(); ++ii)
  {
    per_slave_nodes.push_back( e46[ii] );
    per_master_nodes.push_back( cor[4] );
  }


  for(unsigned int ii=0; ii<e57.size(); ++ii)
  {
    per_slave_nodes.push_back( e57[ii] );
    per_master_nodes.push_back( cor[4] );
  }

  for(unsigned int ii=0; ii<e67.size(); ++ii)
  {
    per_slave_nodes.push_back( e67[ii] );
    per_master_nodes.push_back( cor[4] );
  }

  per_slave_nodes.push_back( cor[5] );
  per_slave_nodes.push_back( cor[6] );
  per_slave_nodes.push_back( cor[7] );

  per_master_nodes.push_back( cor[4] );
  per_master_nodes.push_back( cor[4] );
  per_master_nodes.push_back( cor[4] );

  num_per_nodes = per_slave_nodes.size();

  std::cout<<"-----> Dirchlet nodal BC on Bottom; Top nodes follow cor[4]. \n";
}


// EOF
