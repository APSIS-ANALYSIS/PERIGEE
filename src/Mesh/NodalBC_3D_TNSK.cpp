#include "NodalBC_3D_TNSK.hpp"

NodalBC_3D_TNSK::NodalBC_3D_TNSK(const int &nFunc, const int &nFunc_x, const int &nFunc_y,
            const int &nFunc_z, const int &bc_type)
{
  clock_t log_time = clock();

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
      BC_type_2(nFunc_x, nFunc_y, nFunc_z);
      break;
    case 3:
      BC_type_3();
      break;
    case 4:
      BC_type_4(nFunc_x, nFunc_y, nFunc_z);
      break;
    case 5:
      BC_type_5(nFunc_x, nFunc_y, nFunc_z);
      break;
    default:
      std::cerr<<"Error: NodalBC_3D_TNSK with bc_type = "<<bc_type<<" is not implemented. \n";
      exit(EXIT_FAILURE);
  }

  Create_ID(nFunc);

  log_time = clock() - log_time;
  std::cout<<"===> NodalBC_3D_TNSK, type = "<<bc_type<<" is generated. ";
  std::cout<<"Time taken: "<<((float) log_time)/CLOCKS_PER_SEC<<" seconds. \n"<<std::endl;
}


NodalBC_3D_TNSK::~NodalBC_3D_TNSK()
{}


void NodalBC_3D_TNSK::BC_type_1( const int &nFunc_x, const int &nFunc_y, const int &nFunc_z )
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
  std::cout<<"-----> Dirichlet BC. \n";
}


void NodalBC_3D_TNSK::BC_type_2( const int &nFunc_x, const int &nFunc_y, const int &nFunc_z )
{
  // full periodic. master - slave: front - back, left - right, bottom - top
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  std::vector<int> nodes_edge01, nodes_edge02, nodes_edge13, nodes_edge23,
    nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67, nodes_edge15,
    nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner;

  Generate_BCNode_3D_B(nFunc_x, nFunc_y, nFunc_z, 
      nodes_front, nodes_back, nodes_left, nodes_right,
      nodes_top, nodes_bottom, nodes_edge01, nodes_edge02, nodes_edge13,
      nodes_edge23, nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67,
      nodes_edge15, nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner);  

  unsigned int master_size, slave_size;
  master_size = 0; slave_size = 0;

  master_size += nodes_front.size();  slave_size += nodes_back.size();
  master_size += nodes_left.size();   slave_size += nodes_right.size();
  master_size += nodes_bottom.size(); slave_size += nodes_top.size();

  master_size += nodes_edge01.size(); slave_size += nodes_edge23.size();
  master_size += nodes_edge01.size(); slave_size += nodes_edge67.size();
  master_size += nodes_edge01.size(); slave_size += nodes_edge45.size();

  master_size += nodes_edge15.size(); slave_size += nodes_edge37.size();
  master_size += nodes_edge15.size(); slave_size += nodes_edge26.size();
  master_size += nodes_edge15.size(); slave_size += nodes_edge04.size();

  master_size += nodes_edge13.size(); slave_size += nodes_edge02.size();
  master_size += nodes_edge13.size(); slave_size += nodes_edge46.size();
  master_size += nodes_edge13.size(); slave_size += nodes_edge57.size();

  master_size += (nodes_corner.size() - 1);
  slave_size  += (nodes_corner.size() - 1);

  if(master_size != slave_size)
  {
    std::cerr<<"ERROR: The periodic boundary slave nodes does not match master nodes. \n";
    exit(1);
  }

  num_per_nodes = master_size;

  per_master_nodes.resize(num_per_nodes);
  per_slave_nodes.resize(num_per_nodes);

  unsigned int offset = 0;

  // 1. set face master slave relations
  for(unsigned int ii = 0; ii<nodes_back.size(); ++ii)
  {
    per_master_nodes[offset] = nodes_front[ii];
    per_slave_nodes[offset]  = nodes_back[ii];
    offset += 1;
  }

  for(unsigned int ii=0; ii<nodes_left.size(); ++ii)
  {
    per_master_nodes[offset] = nodes_left[ii];
    per_slave_nodes[offset]  = nodes_right[ii];
    offset += 1;
  }

  for(unsigned int ii=0; ii<nodes_bottom.size(); ++ii)
  {
    per_master_nodes[offset] = nodes_bottom[ii];
    per_slave_nodes[offset]  = nodes_top[ii];
    offset += 1;
  }

  // 2. set edge master slave relations
  // 2.1 master edge01, slave edge 23 edge 67 edge 45
  for(unsigned int ii=0; ii<nodes_edge01.size(); ++ii)
  {
    per_master_nodes[offset] = nodes_edge01[ii];
    per_slave_nodes[offset]  = nodes_edge23[ii];
    offset += 1;
  }

  for(unsigned int ii=0; ii<nodes_edge01.size(); ++ii)
  {
    per_master_nodes[offset] = nodes_edge01[ii];
    per_slave_nodes[offset]  = nodes_edge67[ii];
    offset += 1;
  }

  for(unsigned int ii=0; ii<nodes_edge01.size(); ++ii)
  {
    per_master_nodes[offset] = nodes_edge01[ii];
    per_slave_nodes[offset]  = nodes_edge45[ii];
    offset += 1;
  }

  // 2.2 master edge15 slave edge37 edge26 edge04
  for(unsigned int ii=0; ii<nodes_edge15.size(); ++ii)
  {
    per_master_nodes[offset] = nodes_edge15[ii];
    per_slave_nodes[offset]  = nodes_edge37[ii];
    offset += 1;
  }

  for(unsigned int ii=0; ii<nodes_edge15.size(); ++ii)
  {
    per_master_nodes[offset] = nodes_edge15[ii];
    per_slave_nodes[offset]  = nodes_edge26[ii];
    offset += 1;
  }

  for(unsigned int ii=0; ii<nodes_edge15.size(); ++ii)
  {
    per_master_nodes[offset] = nodes_edge15[ii];
    per_slave_nodes[offset]  = nodes_edge04[ii];
    offset += 1;
  }

  // 2.3 master edge13 slave edge02 edge46 edge57
  for(unsigned int ii=0; ii<nodes_edge13.size(); ++ii)
  {
    per_master_nodes[offset] = nodes_edge13[ii];
    per_slave_nodes[offset]  = nodes_edge02[ii];
    offset += 1;
  }

  for(unsigned int ii=0; ii<nodes_edge13.size(); ++ii)
  {
    per_master_nodes[offset] = nodes_edge13[ii];
    per_slave_nodes[offset]  = nodes_edge46[ii];
    offset += 1;
  }

  for(unsigned int ii=0; ii<nodes_edge13.size(); ++ii)
  {
    per_master_nodes[offset] = nodes_edge13[ii];
    per_slave_nodes[offset]  = nodes_edge57[ii];
    offset += 1;
  }

  // 3 set corner nodes master slave relations
  per_master_nodes[offset] = nodes_corner[1];
  per_slave_nodes[offset]  = nodes_corner[0];
  offset += 1;

  for(unsigned int ii=2; ii<nodes_corner.size(); ++ii)
  {
    per_master_nodes[offset] = nodes_corner[1];
    per_slave_nodes[offset]  = nodes_corner[ii];
    offset += 1;
  }

  VEC_T::shrink2fit(per_master_nodes);
  VEC_T::shrink2fit(per_slave_nodes);
  std::cout<<"-----> full C0 periodic. \n";
}


void NodalBC_3D_TNSK::BC_type_3()
{
  std::cout<<"-----> No strong boundary nodal imposition. \n";
}


void NodalBC_3D_TNSK::BC_type_4( const int &nFunc_x, const int &nFunc_y,
   const int &nFunc_z )
{
  // Nodes for outside faces
  std::vector<int> n_front, n_back, n_left, n_right, n_top, n_bottom;
  std::vector<int> n_edge01, n_edge02, n_edge13, n_edge23,
    n_edge45, n_edge46, n_edge57, n_edge67, n_edge15,
    n_edge37, n_edge04, n_edge26, n_corner;

  Generate_BCNode_3D_B( nFunc_x, nFunc_y, nFunc_z, 
      n_front, n_back, n_left, n_right,
      n_top, n_bottom, n_edge01, n_edge02, n_edge13,
      n_edge23, n_edge45, n_edge46, n_edge57, n_edge67,
      n_edge15, n_edge37, n_edge04, n_edge26, n_corner);  

  // Nodes for interior cube
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  std::vector<int> nodes_edge01, nodes_edge02, nodes_edge13, nodes_edge23,
    nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67, nodes_edge15,
    nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner;

  Generate_BCNode_3D_C( nFunc_x, nFunc_y, nFunc_z, 
      nodes_front, nodes_back, nodes_left, nodes_right,
      nodes_top, nodes_bottom, nodes_edge01, nodes_edge02, nodes_edge13,
      nodes_edge23, nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67,
      nodes_edge15, nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner);  

  // interial master, outside slave
  per_master_nodes.insert(per_master_nodes.end(), nodes_front.begin(), nodes_front.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_front.begin(), n_front.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_back.begin(), nodes_back.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_back.begin(), n_back.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_left.begin(), nodes_left.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_left.begin(), n_left.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_right.begin(), nodes_right.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_right.begin(), n_right.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_top.begin(), nodes_top.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_top.begin(), n_top.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_bottom.begin(), nodes_bottom.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_bottom.begin(), n_bottom.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge01.begin(), nodes_edge01.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge01.begin(), n_edge01.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge02.begin(), nodes_edge02.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge02.begin(), n_edge02.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge13.begin(), nodes_edge13.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge13.begin(), n_edge13.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge23.begin(), nodes_edge23.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge23.begin(), n_edge23.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge45.begin(), nodes_edge45.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge45.begin(), n_edge45.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge46.begin(), nodes_edge46.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge46.begin(), n_edge46.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge57.begin(), nodes_edge57.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge57.begin(), n_edge57.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge67.begin(), nodes_edge67.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge67.begin(), n_edge67.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge15.begin(), nodes_edge15.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge15.begin(), n_edge15.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge37.begin(), nodes_edge37.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge37.begin(), n_edge37.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge04.begin(), nodes_edge04.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge04.begin(), n_edge04.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge26.begin(), nodes_edge26.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge26.begin(), n_edge26.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_corner.begin(), nodes_corner.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_corner.begin(), n_corner.end());


  num_per_nodes = per_master_nodes.size();
  assert(num_per_nodes == per_slave_nodes.size());

  num_dir_nodes = 0;
  std::cout<<"-----> Grad c dot n = 0 for C1 function. \n";
}


void NodalBC_3D_TNSK::BC_type_5( const int &nFunc_x, const int &nFunc_y,
   const int &nFunc_z )
{
  std::vector<int> n_front, n_back, n_left, n_right, n_top, n_bottom;
  std::vector<int> n_edge01, n_edge02, n_edge13, n_edge23,
    n_edge45, n_edge46, n_edge57, n_edge67, n_edge15,
    n_edge37, n_edge04, n_edge26, n_corner;

  Generate_BCNode_3D_B(nFunc_x, nFunc_y, nFunc_z, 
      n_front, n_back, n_left, n_right,
      n_top, n_bottom, n_edge01, n_edge02, n_edge13,
      n_edge23, n_edge45, n_edge46, n_edge57, n_edge67,
      n_edge15, n_edge37, n_edge04, n_edge26, n_corner);  


  dir_nodes.insert(dir_nodes.end(), n_top.begin(), n_top.end());
  dir_nodes.insert(dir_nodes.end(), n_bottom.begin(), n_bottom.end());

  dir_nodes.insert(dir_nodes.end(), n_edge01.begin(), n_edge01.end());
  dir_nodes.insert(dir_nodes.end(), n_edge02.begin(), n_edge02.end());
  dir_nodes.insert(dir_nodes.end(), n_edge13.begin(), n_edge13.end());
  dir_nodes.insert(dir_nodes.end(), n_edge23.begin(), n_edge23.end());

  dir_nodes.insert(dir_nodes.end(), n_edge46.begin(), n_edge46.end());
  dir_nodes.insert(dir_nodes.end(), n_edge45.begin(), n_edge45.end());
  dir_nodes.insert(dir_nodes.end(), n_edge57.begin(), n_edge57.end());
  dir_nodes.insert(dir_nodes.end(), n_edge67.begin(), n_edge67.end());

  dir_nodes.insert(dir_nodes.end(), n_corner.begin(), n_corner.end());

  num_dir_nodes = dir_nodes.size();
  num_per_nodes = 0;
  std::cout<<"-----> Dirichlet on Top and Bottom surfaces. \n";
}



// EOF
