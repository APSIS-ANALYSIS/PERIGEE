#include "BoundaryCond2D.hpp"

BoundaryCond2D::BoundaryCond2D( const IMesh * const &mesh, const int &bc_type )
{
  clock_t log_time = clock();

  // Initialize data
  clear_nodes(); clear_elems();

  // Generate boundary nodes/elems according to the user definition
  switch(bc_type)
  {
    case 1:
      BC_type_1(mesh);
      break;
    case 2:
      BC_type_2(mesh);
      break;
    case 3:
      BC_type_3(mesh);
      break;
    case 4:
      BC_type_4(mesh);
      break;
    case 5:
      BC_type_5(mesh);
      break;
    case 6:
      BC_type_6(mesh);
      break;
    case 7:
      BC_type_7(mesh);
      break;
    case 8:
      BC_type_8(mesh);
      break;
    case 9:
      BC_type_9(mesh);
      break;
    default:
      std::cerr<<"Error: BoundaryCond2D the bc_type "<<bc_type<<" is not implemented!"<<std::endl;
      exit(1);
  }

  // Generate ID array
  Create_ID(mesh);

  // Calculate time and memory usage by this constructor
  std::cout<<"=== Boundary condition with type "<<bc_type<<" generated. ";
  std::cout<<"Time taken: "<<((float) log_time)/CLOCKS_PER_SEC<<" seconds"<<std::endl;
}



BoundaryCond2D::~BoundaryCond2D()
{}



void BoundaryCond2D::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"======== BC info ======="<<std::endl;
  if(num_dir_nodes > 0)
  {
    std::cout<<"Dirichlet nodes: "<<std::endl;
    for(unsigned int ii=0; ii<num_dir_nodes; ++ii)
      std::cout<<dir_nodes[ii]<<'\t';
    std::cout<<std::endl;
  }

  if(num_per_nodes>0)
  {
    std::cout<<"Periodic master - slave nodes: "<<std::endl;
    for(unsigned int ii=0; ii<num_per_nodes; ++ii)
      std::cout<<per_master_nodes[ii]<<'\t'<<per_slave_nodes[ii]<<std::endl;
  }

  std::cout<<std::endl<<"ID array: "<<std::endl;
  for(unsigned int ii=0; ii<ID.size(); ++ii)
    std::cout<<ID[ii]<<'\n';
  std::cout<<std::endl;

  std::cout<<std::endl<<"bottom elem: \t";
  for(unsigned int ii=0; ii<get_num_bottom_elem(); ++ii)
    std::cout<<get_bottom_elem(ii)<<'\t';
  std::cout<<std::endl<<"top elem: \t";
  for(unsigned int ii=0; ii<get_num_top_elem(); ++ii)
    std::cout<<get_top_elem(ii)<<'\t';
  std::cout<<std::endl<<"left elem: \t";
  for(unsigned int ii=0; ii<get_num_left_elem(); ++ii)
    std::cout<<get_left_elem(ii)<<'\t';
  std::cout<<std::endl<<"right elem: \t";
  for(unsigned int ii=0; ii<get_num_right_elem(); ++ii)
    std::cout<<get_right_elem(ii)<<'\t';
  std::cout<<std::endl<<"========================"<<std::endl;
}



void BoundaryCond2D::print_edge_cornder_nodes(const IMesh * const &mesh) const
{}


// clear_nodes(): 
// clear the std::std::vector: dir_nodes, per_slave_nodes, per_master_nodes
// set num_dir_nodes, num_per_nodes = 0               
void BoundaryCond2D::clear_nodes()
{
  dir_nodes.clear();
  num_dir_nodes = 0;

  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;
}


// clear_elems():
void BoundaryCond2D::clear_elems()
{
  left_elem.clear();   num_left_elem = 0;
  right_elem.clear();  num_right_elem = 0;
  top_elem.clear();    num_top_elem = 0;
  bottom_elem.clear(); num_bottom_elem = 0;
}


void BoundaryCond2D::Create_ID(const IMesh * const &mesh)
{
  int num_node = mesh->get_nFunc();
  ID.resize(num_node);
  
  std::vector<int>::iterator it;
  std::vector<int> temp_dir_nodes;
  temp_dir_nodes.clear();
  temp_dir_nodes.insert(temp_dir_nodes.end(), dir_nodes.begin(), dir_nodes.end());
  sort(temp_dir_nodes.begin(), temp_dir_nodes.end());
  for(int ii=0; ii<num_node; ++ii)
  {
    if( binary_search(temp_dir_nodes.begin(), temp_dir_nodes.end(), ii) )
      ID[ii] = -1;
    else
    {
      it = find(per_slave_nodes.begin(), per_slave_nodes.end(), ii);
      if(it == per_slave_nodes.end())
        ID[ii] = ii;
      else
        ID[ii] = per_master_nodes[it - per_slave_nodes.begin()];
    }
  }
  VEC_T::shrink2fit(ID);
}


void BoundaryCond2D::Generate_BCNodes_A( const IMesh * const &mesh,
    std::vector<int> &left, std::vector<int> &right, std::vector<int> &top, 
    std::vector<int> &bottom, std::vector<int> &corner ) const
{
  left.clear(); right.clear(); top.clear(); bottom.clear();

  int nFunc_x = mesh->get_nFunc_x();
  int nFunc_y = mesh->get_nFunc_y();

  for(int xx=1; xx<nFunc_x -1; ++xx)
  {
    bottom.push_back(xx);
    top.push_back(xx + nFunc_x * (nFunc_y - 1));
  }

  for(int yy=1; yy<nFunc_y-1; ++yy)
  {
    left.push_back(yy*nFunc_x);
    right.push_back((yy+1)*nFunc_x - 1);
  }

  corner.push_back(0); corner.push_back(nFunc_x - 1);
  corner.push_back(nFunc_x * (nFunc_y-1));
  corner.push_back(nFunc_x * nFunc_y - 1);

  VEC_T::shrink2fit(left);
  VEC_T::shrink2fit(right);
  VEC_T::shrink2fit(top);
  VEC_T::shrink2fit(bottom);
  VEC_T::shrink2fit(corner);
}


void BoundaryCond2D::Generate_BCNodes_B( const IMesh * const &mesh,
    std::vector<int> &lef1, std::vector<int> &rig1, std::vector<int> &top1, std::vector<int> &bot1,
    std::vector<int> &lef2, std::vector<int> &rig2, std::vector<int> &top2, std::vector<int> &bot2,
    std::vector<int> &corner1, std::vector<int> &corner2,
    std::vector<int> &corner3, std::vector<int> &corner4 ) const
{
  lef1.clear(); rig1.clear(); top1.clear(); bot1.clear();
  lef2.clear(); rig2.clear(); top2.clear(); bot2.clear();
  corner1.clear(); corner2.clear(); corner3.clear(); corner4.clear();

  int nFunc_x = mesh->get_nFunc_x();
  int nFunc_y = mesh->get_nFunc_y();

  for(int xx=2; xx<nFunc_x - 2; ++xx)
  {
    bot1.push_back(xx);
    bot2.push_back(xx + nFunc_x);
    top1.push_back(xx + nFunc_x * (nFunc_y - 1));
    top2.push_back(xx + nFunc_x * (nFunc_y - 2));
  }

  for(int yy=2; yy<nFunc_y - 2; ++yy)
  {
    lef1.push_back(yy*nFunc_x);
    lef2.push_back(yy*nFunc_x + 1);
    rig1.push_back((yy+1)*nFunc_x - 1);
    rig2.push_back((yy+1)*nFunc_x - 2);
  }

  corner1.push_back(0); 
  corner1.push_back(nFunc_x - 1);
  corner1.push_back(nFunc_x * (nFunc_y-1));
  corner1.push_back(nFunc_x * nFunc_y - 1);

  corner2.push_back(1);
  corner2.push_back(nFunc_x - 2); 
  corner2.push_back(nFunc_x * (nFunc_y - 1) + 1);
  corner2.push_back(nFunc_x * nFunc_y - 2);

  corner3.push_back(nFunc_x);
  corner3.push_back(nFunc_x * (nFunc_y-2));
  corner3.push_back(2*nFunc_x - 1);
  corner3.push_back(nFunc_x*(nFunc_y-1)-1);

  corner4.push_back(nFunc_x  + 1);
  corner4.push_back(2*nFunc_x - 2);
  corner4.push_back(nFunc_x * (nFunc_y-2)+1);
  corner4.push_back(nFunc_x * (nFunc_y-1)-2);

  VEC_T::shrink2fit(lef1); VEC_T::shrink2fit(lef2);
  VEC_T::shrink2fit(rig1); VEC_T::shrink2fit(rig2);
  VEC_T::shrink2fit(top1); VEC_T::shrink2fit(top2);
  VEC_T::shrink2fit(bot1); VEC_T::shrink2fit(bot2);

  VEC_T::shrink2fit(corner1); VEC_T::shrink2fit(corner2);
  VEC_T::shrink2fit(corner3); VEC_T::shrink2fit(corner4);
}


void BoundaryCond2D::Generate_BCElems_A( const IMesh * const &mesh,
    std::vector<int> &lef, std::vector<int> &rig, std::vector<int> &bot,
    std::vector<int> &top ) const
{
  lef.clear(); rig.clear(); bot.clear(); top.clear();

  int nElem_x = mesh->get_nElem_x();
  int nElem_y = mesh->get_nElem_y();

  for(int ii=0; ii<nElem_y; ++ii)
  {
    lef.push_back( ii * nElem_x );
    rig.push_back( (ii + 1)*nElem_x - 1  );
  }

  for(int ii=0; ii<nElem_x; ++ii)
  {
    bot.push_back(ii);
    top.push_back(ii + nElem_x * (nElem_y-1));
  }
}


void BoundaryCond2D::BC_type_1(const IMesh * const &mesh)
{
  std::vector<int> lef, rig, top, bot, cor;
  Generate_BCNodes_A(mesh, lef, rig, top, bot, cor);
  dir_nodes.clear();
  dir_nodes.insert(dir_nodes.end(), cor.begin(), cor.end());
  dir_nodes.insert(dir_nodes.end(), bot.begin(), bot.end());
  dir_nodes.insert(dir_nodes.end(), lef.begin(), lef.end());
  dir_nodes.insert(dir_nodes.end(), rig.begin(), rig.end());
  dir_nodes.insert(dir_nodes.end(), top.begin(), top.end());

  num_dir_nodes = dir_nodes.size();
  
  num_per_nodes = 0;

  std::cout<<"-----> C0 Dirichlet BC. \n";
}


void BoundaryCond2D::BC_type_2(const IMesh * const &mesh)
{
  std::vector<int> lef1, lef2, rig1, rig2, top1, top2, bot1, bot2;
  std::vector<int> cor1, cor2, cor3, cor4;

  Generate_BCNodes_B(mesh, lef1, rig1, top1, bot1, lef2, 
      rig2, top2, bot2, cor1, cor2, cor3, cor4);

  dir_nodes.clear(); num_dir_nodes = 0;

  per_slave_nodes.clear(); per_master_nodes.clear();

  // the inner nodes are set as master, the outer bc nodes are slave
  per_slave_nodes.insert(per_slave_nodes.end(), lef1.begin(), lef1.end());
  per_master_nodes.insert(per_master_nodes.end(), lef2.begin(), lef2.end());
  per_slave_nodes.insert(per_slave_nodes.end(), rig1.begin(), rig1.end());
  per_master_nodes.insert(per_master_nodes.end(), rig2.begin(), rig2.end());
  per_slave_nodes.insert(per_slave_nodes.end(), bot1.begin(), bot1.end());
  per_master_nodes.insert(per_master_nodes.end(), bot2.begin(), bot2.end());
  per_slave_nodes.insert(per_slave_nodes.end(), top1.begin(), top1.end());
  per_master_nodes.insert(per_master_nodes.end(), top2.begin(), top2.end());

  per_master_nodes.push_back(cor4[0]);
  per_master_nodes.push_back(cor4[0]);
  per_master_nodes.push_back(cor4[0]);
  per_slave_nodes.push_back(cor1[0]);
  per_slave_nodes.push_back(cor2[0]);
  per_slave_nodes.push_back(cor3[0]);
  
  per_master_nodes.push_back(cor4[1]);
  per_master_nodes.push_back(cor4[1]);
  per_master_nodes.push_back(cor4[1]);
  per_slave_nodes.push_back(cor1[1]);
  per_slave_nodes.push_back(cor2[1]);
  per_slave_nodes.push_back(cor3[2]);
  
  per_master_nodes.push_back(cor4[2]);
  per_master_nodes.push_back(cor4[2]);
  per_master_nodes.push_back(cor4[2]);
  per_slave_nodes.push_back(cor1[2]);
  per_slave_nodes.push_back(cor2[2]);
  per_slave_nodes.push_back(cor3[1]);
  
  per_master_nodes.push_back(cor4[3]);
  per_master_nodes.push_back(cor4[3]);
  per_master_nodes.push_back(cor4[3]);
  per_slave_nodes.push_back(cor1[3]);
  per_slave_nodes.push_back(cor2[3]);
  per_slave_nodes.push_back(cor3[3]);

  num_per_nodes = per_slave_nodes.size();
  std::cout<<"-----> C1 grad c dot n = 0 \n";
}


void BoundaryCond2D::BC_type_3(const IMesh * const &mesh)
{
  std::vector<int> lef, rig, bot, top, cor;
  Generate_BCNodes_A(mesh, lef, rig, top, bot, cor);
  
  dir_nodes.clear(); num_dir_nodes = 0;

  per_master_nodes.insert(per_master_nodes.end(), bot.begin(), bot.end());
  per_slave_nodes.insert(per_slave_nodes.end(), top.begin(), top.end());

  per_master_nodes.insert(per_master_nodes.end(), lef.begin(), lef.end());
  per_slave_nodes.insert(per_slave_nodes.end(), rig.begin(), rig.end());

  per_master_nodes.push_back(cor[0]);
  per_master_nodes.push_back(cor[0]);
  per_master_nodes.push_back(cor[0]);

  per_slave_nodes.push_back(cor[1]);
  per_slave_nodes.push_back(cor[2]);
  per_slave_nodes.push_back(cor[3]);

  num_per_nodes = per_slave_nodes.size();

  std::cout<<"-----> C0 periodic strong nodal imposition. \n";
}


void BoundaryCond2D::BC_type_4(const IMesh * const &mesh)
{
  std::vector<int> lef1, lef2, rig1, rig2, top1, top2, bot1, bot2;
  std::vector<int> cor1, cor2, cor3, cor4;

  Generate_BCNodes_B(mesh, lef1, rig1, top1, bot1, lef2, 
      rig2, top2, bot2, cor1, cor2, cor3, cor4);

  dir_nodes.clear(); num_dir_nodes = 0;

  per_slave_nodes.clear(); per_master_nodes.clear();

  per_slave_nodes.insert(per_slave_nodes.end(), rig1.begin(), rig1.end());
  per_master_nodes.insert(per_master_nodes.end(), lef2.begin(), lef2.end());
  
  per_slave_nodes.insert(per_slave_nodes.end(), rig2.begin(), rig2.end());
  per_master_nodes.insert(per_master_nodes.end(), lef1.begin(), lef1.end());
  
  per_slave_nodes.insert(per_slave_nodes.end(), top1.begin(), top1.end());
  per_master_nodes.insert(per_master_nodes.end(), bot2.begin(), bot2.end());
  
  per_slave_nodes.insert(per_slave_nodes.end(), top2.begin(), top2.end());
  per_master_nodes.insert(per_master_nodes.end(), bot1.begin(), bot1.end());

  per_master_nodes.push_back(cor1[0]); // 0
  per_master_nodes.push_back(cor1[0]); 
  per_master_nodes.push_back(cor1[0]);
  per_slave_nodes.push_back(cor2[1]);  // 3
  per_slave_nodes.push_back(cor3[1]);  //15
  per_slave_nodes.push_back(cor4[3]);  //18

  per_master_nodes.push_back(cor2[0]); // 1
  per_master_nodes.push_back(cor2[0]);
  per_master_nodes.push_back(cor2[0]);
  per_slave_nodes.push_back(cor1[1]);  // 4 
  per_slave_nodes.push_back(cor4[2]);  // 16
  per_slave_nodes.push_back(cor3[3]);  // 19

  per_master_nodes.push_back(cor3[0]); // 5
  per_master_nodes.push_back(cor3[0]);
  per_master_nodes.push_back(cor3[0]);
  per_slave_nodes.push_back(cor1[2]);  // 20
  per_slave_nodes.push_back(cor4[1]);  // 8  
  per_slave_nodes.push_back(cor2[3]);  // 23

  per_master_nodes.push_back(cor4[0]); // 6
  per_master_nodes.push_back(cor4[0]);
  per_master_nodes.push_back(cor4[0]);
  per_slave_nodes.push_back(cor3[2]);  // 9
  per_slave_nodes.push_back(cor2[2]);  // 21  
  per_slave_nodes.push_back(cor1[3]);  // 24

  num_per_nodes = per_slave_nodes.size();

  std::cout<<"-----> C1 periodic strong nodal imposition. \n";
  std::cout<<"       Periodic splines is needed here for this bc. \n";
}


void BoundaryCond2D::BC_type_5(const IMesh * const &mesh)
{
  Generate_BCElems_A(mesh, left_elem, right_elem, bottom_elem, top_elem);

  num_left_elem = left_elem.size();
  num_right_elem = right_elem.size();
  num_bottom_elem = bottom_elem.size();
  num_top_elem = top_elem.size();
  std::cout<<"-----> Weak element bc imposition with boundary faces specified. \n";
}


void BoundaryCond2D::BC_type_6(const IMesh * const &mesh)
{
  num_dir_nodes = 0;
  dir_nodes.clear();

  num_per_nodes = 0;
  per_slave_nodes.clear();
  per_master_nodes.clear();

  std::cout<<"-----> adiabatic bc. \n";
}


void BoundaryCond2D::BC_type_7(const IMesh * const &mesh)
{
  std::vector<int> lef, rig, top, bot, cor;
  Generate_BCNodes_A(mesh, lef, rig, top, bot, cor);
  dir_nodes.clear();
  dir_nodes.insert(dir_nodes.end(), cor.begin(), cor.end());
  dir_nodes.insert(dir_nodes.end(), bot.begin(), bot.end());
  dir_nodes.insert(dir_nodes.end(), top.begin(), top.end());

  num_dir_nodes = dir_nodes.size();
  
  num_per_nodes = 0;

  std::cout<<"-----> C0 Dirichlet BC on Top, Bottom, Corner; ";
  std::cout<<"No BC on Left / Right. \n ";
}

void BoundaryCond2D::BC_type_8(const IMesh * const &mesh)
{
  std::vector<int> lef, rig, top, bot, cor;
  Generate_BCNodes_A(mesh, lef, rig, top, bot, cor);
  dir_nodes.clear();
  dir_nodes.insert(dir_nodes.end(), cor.begin(), cor.end());
  dir_nodes.insert(dir_nodes.end(), lef.begin(), lef.end());
  dir_nodes.insert(dir_nodes.end(), rig.begin(), rig.end());

  num_dir_nodes = dir_nodes.size();
  
  num_per_nodes = 0;

  std::cout<<"-----> C0 Dirichlet BC on Left, Right, and Corner; ";
  std::cout<<"No BC on Top / Bottom. \n ";
}

void BoundaryCond2D::BC_type_9(const IMesh * const &mesh)
{
  std::vector<int> lef, rig, top, bot, cor;
  Generate_BCNodes_A(mesh, lef, rig, top, bot, cor);
  dir_nodes.clear();
  dir_nodes.insert(dir_nodes.end(), bot.begin(), bot.end());
  dir_nodes.push_back(cor[0]);
  dir_nodes.push_back(cor[1]);

  num_dir_nodes = dir_nodes.size();
  
  num_per_nodes = 0;

  std::cout<<"-----> C0 Dirichlet BC on Bottom with corner; ";
  std::cout<<"No BC on Top / Left /Right. \n ";
}

//EOF
