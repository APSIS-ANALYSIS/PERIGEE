#include "BoundaryCond.hpp"

BoundaryCond::BoundaryCond( const IMesh * const &mesh, const int &bc_type)
: bctype(bc_type)
{
  clock_t log_time = clock();

  // Initialize the data
  Clear_nodes(); Clear_elems();

  // Generate nodes/elems according to the bc_type
  switch(bc_type)
  {
    case -2:
      BC_test_2();
      break;
    case -1:
      BC_test_1();
      break;
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
    case 10:
      BC_type_10(mesh);
      break;
    case 11:
      BC_type_11(mesh);
      break;
    case 12:
      BC_type_12(mesh);
      break;
    case 13:
      BC_type_13(mesh);
      break;
    case 14:
      BC_type_14(mesh);
      break;
    case 15:
      BC_type_15(mesh);
      break;
    case 16:
      BC_type_16(mesh);
      break;
    case 17:
      BC_type_17(mesh);
      break;
    case 18:
      BC_type_18(mesh);
      break;
    case 19:
      BC_type_19(mesh);
      break;
    case 20:
      BC_type_20(mesh);
      break;
    case 21:
      BC_type_21(mesh);
      break;
    case 22:
      BC_type_22(mesh);
      break;
    case 23:
      BC_type_23(mesh);
      break;
    case 101:
      BC_type_101(mesh);
      break;
    case 102:
      BC_type_102(mesh);
      break;
    case 103:
      BC_type_103(mesh);
      break;
    default:
      std::cerr<<"WARNING: the bc_type "<<bc_type<<" is not implemented!"<<std::endl;
      exit(1);
  }

  // Generate ID array
  Create_ID(mesh->get_nFunc());


  log_time = clock() - log_time;

  // Print info on screen
  std::cout<<"=== Boundary condition with type "<<bc_type<<" generated, ";
  std::cout<<"taking "<<((double) log_time)/CLOCKS_PER_SEC<<" seconds"<<std::endl<<std::endl;
}



BoundaryCond::BoundaryCond( const IMesh * const &mesh, const int &bc_type,
   const std::vector<double> &cPts )
: bctype(bc_type)
{
  clock_t log_time = clock();

  // Initialize the data
  Clear_nodes(); Clear_elems();

  // Generate nodes/elems according to the bc_type
  switch(bc_type)
  {
    case 104:
      BC_type_104(mesh, cPts);
      break;
    case 105:
      BC_type_105(mesh, cPts);
      break;
    case 200:
      BC_type_200(mesh, cPts);
      break;
    case 201:
      BC_type_201(mesh, cPts);
      break;
    default:
      std::cerr<<"WARNING: the bc_type "<<bc_type<<" is not implemented!"<<std::endl;
      exit(1);
  }

  // Generate ID array
  Create_ID(mesh->get_nFunc());

  log_time = clock() - log_time;

  // Print info on screen
  std::cout<<"=== Boundary condition with type "<<bc_type<<" generated, ";
  std::cout<<"taking "<<((double) log_time)/CLOCKS_PER_SEC<<" seconds"<<std::endl<<std::endl;
}


BoundaryCond::~BoundaryCond()
{ 
  VEC_T::clean(dir_nodes);
  VEC_T::clean(per_slave_nodes);
  VEC_T::clean(per_master_nodes);
  VEC_T::clean(front_elem);
  VEC_T::clean(back_elem);
  VEC_T::clean(left_elem);
  VEC_T::clean(right_elem);
  VEC_T::clean(top_elem);
  VEC_T::clean(bottom_elem);
  std::cout<<"-- BoundaryCond "<<bctype<<" deleted. \n";
}


void BoundaryCond::Create_ID(const IMesh * const &mesh)
{
  // ID array records the boundary nodes that has strong imposition of
  // conditions, i.e., dirichlet nodes and slave periodic nodes.
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
      {
        ID[ii] = per_master_nodes[it - per_slave_nodes.begin()];
      }
    }
  }
  VEC_T::shrink2fit(ID);
}


void BoundaryCond::Create_ID(const int &in_nfunc)
{
  ID.clear();
  ID.resize(in_nfunc);
  VEC_T::shrink2fit(ID);
  for(unsigned int ii = 0; ii<ID.size(); ++ii)
    ID[ii] = ii;
  for(unsigned int ii = 0; ii<per_slave_nodes.size(); ++ii)
    ID[per_slave_nodes[ii]] = per_master_nodes[ii];
  for(unsigned int ii = 0; ii<dir_nodes.size(); ++ii)
    ID[dir_nodes[ii]] = -1;
}


void BoundaryCond::print_info() const
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
    std::cout<<ID[ii]<<'\t';
  std::cout<<std::endl;

  std::cout<<std::endl<<"front elem: \t";
  for(unsigned int ii=0; ii<get_num_front_elem(); ++ii)
    std::cout<<get_front_elem(ii)<<'\t';
  std::cout<<std::endl<<"back elem: \t";
  for(unsigned int ii=0; ii<get_num_back_elem(); ++ii)
    std::cout<<get_back_elem(ii)<<'\t';
  std::cout<<std::endl<<"left elem: \t";
  for(unsigned int ii=0; ii<get_num_left_elem(); ++ii)
    std::cout<<get_left_elem(ii)<<'\t';
  std::cout<<std::endl<<"right elem: \t";
  for(unsigned int ii=0; ii<get_num_right_elem(); ++ii)
    std::cout<<get_right_elem(ii)<<'\t';
  std::cout<<std::endl<<"top elem: \t";
  for(unsigned int ii=0; ii<get_num_top_elem(); ++ii)
    std::cout<<get_top_elem(ii)<<'\t';
  std::cout<<std::endl<<"bottom elem: \t";
  for(unsigned int ii=0; ii<get_num_bottom_elem(); ++ii)
    std::cout<<get_bottom_elem(ii)<<'\t';
  std::cout<<std::endl<<"========================"<<std::endl;
}


void BoundaryCond::Clear_nodes()
{
  dir_nodes.clear(); 
  num_dir_nodes = 0;
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;
}


void BoundaryCond::Clear_elems()
{
  front_elem.clear();
  back_elem.clear();
  left_elem.clear();
  right_elem.clear();
  top_elem.clear();
  bottom_elem.clear();

  num_front_elem = 0; num_back_elem = 0;
  num_left_elem = 0; num_right_elem = 0;
  num_top_elem = 0; num_bottom_elem = 0;
}


void BoundaryCond::Generate_BCNodes_A( const IMesh * const &mesh,
    std::vector<int> &front, std::vector<int> &back, std::vector<int> &left,
    std::vector<int> &right, std::vector<int> &top, std::vector<int> &bottom ) const
{
  // This function reads in basic mesh object, based on the nfunc,
  // nFunc_x, ... , nElem, etc. to determine the boundary nodes/elements.
  // This approach is restricted to NURBS single patch mesh.

  front.clear(); back.clear(); left.clear(); right.clear();
  top.clear(); bottom.clear();

  const int nFunc_x = mesh->get_nFunc_x();
  const int nFunc_y = mesh->get_nFunc_y();
  const int nFunc_z = mesh->get_nFunc_z();

  for(int zz=0; zz<nFunc_z; ++zz)
  {
    for(int yy=0; yy<nFunc_y; ++yy)
    {
      int node = zz * nFunc_x * nFunc_y+ yy * nFunc_x;
      back.push_back(node);
      front.push_back(node + nFunc_x - 1);
    }
  }
  VEC_T::shrink2fit( front );
  VEC_T::shrink2fit( back );

  for(int zz=0; zz<nFunc_z; ++zz)
  {
    for(int xx=1; xx<nFunc_x-1; ++xx)
    {
      int node = zz * nFunc_x*nFunc_y + xx;
      left.push_back(node);
      right.push_back(node + nFunc_x*(nFunc_y-1));
    }
  }
  VEC_T::shrink2fit(left);
  VEC_T::shrink2fit(right);

  for(int yy = 1; yy<nFunc_y - 1; ++yy)
  {
    for(int xx = 1; xx<nFunc_x - 1; ++xx)
    {
      int node = yy * nFunc_x + xx;
      bottom.push_back(node);
      top.push_back(node + (nFunc_z-1)*nFunc_x*nFunc_y );
    }
  }
  VEC_T::shrink2fit(top);
  VEC_T::shrink2fit(bottom);
}

void BoundaryCond::Generate_BCNodes_B( const IMesh * const &mesh,
    std::vector<int> &front, std::vector<int> &back, std::vector<int> &left,
    std::vector<int> &right, std::vector<int> &top, std::vector<int> &bottom,
    std::vector<int> &edge01, std::vector<int> &edge02, std::vector<int> &edge13,
    std::vector<int> &edge23, std::vector<int> &edge45, std::vector<int> &edge46,
    std::vector<int> &edge57, std::vector<int> &edge67, std::vector<int> &edge15,
    std::vector<int> &edge37, std::vector<int> &edge04, std::vector<int> &edge26,
    std::vector<int> &corner ) const
{
  front.clear(); back.clear(); left.clear(); right.clear(); top.clear(); bottom.clear();
  edge01.clear(); edge02.clear(); edge13.clear(); edge23.clear();
  edge45.clear(); edge46.clear(); edge57.clear(); edge67.clear();
  edge15.clear(); edge37.clear(); edge04.clear(); edge26.clear();
  corner.clear();

  const int nFunc_x = mesh->get_nFunc_x();
  const int nFunc_y = mesh->get_nFunc_y();
  const int nFunc_z = mesh->get_nFunc_z();

  for( int zz=1; zz<nFunc_z-1; ++zz )
  {
    for( int yy=1; yy<nFunc_y-1; ++yy)
    {
      int node = zz * nFunc_x * nFunc_y + yy * nFunc_x;
      back.push_back(node);
      front.push_back(node + nFunc_x - 1);
    }
    for(int xx=1; xx<nFunc_x-1; ++xx)
    {
      int node = zz * nFunc_x * nFunc_y + xx;
      left.push_back(node);
      right.push_back(node + nFunc_x * (nFunc_y - 1));
    }
  }

  for(int yy=1; yy<nFunc_y-1; ++yy)
  {
    for(int xx=1; xx<nFunc_x - 1; ++xx)
    {
      int node = yy * nFunc_x + xx;
      bottom.push_back(node);
      top.push_back(node + (nFunc_z-1)*nFunc_x * nFunc_y);
    }
  }
  VEC_T::shrink2fit(front); VEC_T::shrink2fit(back); VEC_T::shrink2fit(left);
  VEC_T::shrink2fit(right); VEC_T::shrink2fit(top); VEC_T::shrink2fit(bottom);

  for(int xx = 1; xx<nFunc_x-1; ++xx)
  {
    edge01.push_back(xx);
    edge23.push_back((nFunc_y-1)*nFunc_x + xx);
    edge45.push_back((nFunc_z-1)*nFunc_x*nFunc_y + xx);
    edge67.push_back((nFunc_z-1)*nFunc_x*nFunc_y + (nFunc_y-1)*nFunc_x + xx);
  }
  VEC_T::shrink2fit(edge01); VEC_T::shrink2fit(edge23); VEC_T::shrink2fit(edge45); 
  VEC_T::shrink2fit(edge67);

  for(int yy=1; yy<nFunc_y-1; ++yy)
  {
    edge02.push_back(yy*nFunc_x);
    edge13.push_back(yy*nFunc_x + nFunc_x -1);
    edge46.push_back((nFunc_z-1)*nFunc_x*nFunc_y + yy*nFunc_x);
    edge57.push_back((nFunc_z-1)*nFunc_x*nFunc_y + yy*nFunc_x + nFunc_x - 1);
  }
  VEC_T::shrink2fit(edge02); VEC_T::shrink2fit(edge13); 
  VEC_T::shrink2fit(edge46); VEC_T::shrink2fit(edge57);

  for(int zz=1; zz<nFunc_z-1; ++zz)
  {
    edge04.push_back(zz*nFunc_x*nFunc_y);
    edge15.push_back(zz*nFunc_x*nFunc_y + nFunc_x - 1);
    edge26.push_back(zz*nFunc_x*nFunc_y + (nFunc_y-1)*nFunc_x);
    edge37.push_back(zz*nFunc_x*nFunc_y + (nFunc_y-1)*nFunc_x + nFunc_x - 1);
  }
  VEC_T::shrink2fit(edge04); VEC_T::shrink2fit(edge26); 
  VEC_T::shrink2fit(edge15); VEC_T::shrink2fit(edge37);

  corner.push_back(0);
  corner.push_back(nFunc_x-1);
  corner.push_back((nFunc_y-1)*nFunc_x);
  corner.push_back((nFunc_y-1)*nFunc_x + nFunc_x - 1);
  corner.push_back((nFunc_z-1)*nFunc_x*nFunc_y);
  corner.push_back((nFunc_z-1)*nFunc_x*nFunc_y + nFunc_x -1);
  corner.push_back((nFunc_z-1)*nFunc_x*nFunc_y + (nFunc_y-1)*nFunc_x);
  corner.push_back((nFunc_z-1)*nFunc_x*nFunc_y + (nFunc_y-1)*nFunc_x + nFunc_x - 1);
  VEC_T::shrink2fit(corner);
}


void BoundaryCond::Generate_BCNodes_C( const IMesh * const &mesh,
    std::vector<int> &front, std::vector<int> &back, std::vector<int> &left,
    std::vector<int> &right, std::vector<int> &top, std::vector<int> &bottom,
    std::vector<int> &edge01, std::vector<int> &edge02, std::vector<int> &edge13,
    std::vector<int> &edge23, std::vector<int> &edge45, std::vector<int> &edge46,
    std::vector<int> &edge57, std::vector<int> &edge67, std::vector<int> &edge15,
    std::vector<int> &edge37, std::vector<int> &edge04, std::vector<int> &edge26,
    std::vector<int> &corner    ) const
{
  front.clear(); back.clear(); left.clear(); right.clear(); top.clear(); bottom.clear();
  edge01.clear(); edge02.clear(); edge13.clear(); edge23.clear();
  edge45.clear(); edge46.clear(); edge57.clear(); edge67.clear();
  edge15.clear(); edge37.clear(); edge04.clear(); edge26.clear();
  corner.clear();

  const int nFunc_x = mesh->get_nFunc_x();
  const int nFunc_y = mesh->get_nFunc_y();
  const int nFunc_z = mesh->get_nFunc_z();

  const int nFxy = nFunc_x * nFunc_y;

  // Genearte the face nodes, note the nodes here includes the edge and corner
  // nodes. I will specify the edge and corner later as well.
  int node;
  for( int zz = 1; zz<nFunc_z-1; ++zz )
  {
    for( int yy=1; yy<nFunc_y-1; ++yy )
    {
      node = zz * nFxy + yy * nFunc_x;
      back.push_back(node+1);
      front.push_back(node + nFunc_x - 2);
    }
    for( int xx=1; xx<nFunc_x - 1; ++xx )
    {
      node = zz * nFxy + xx;
      left.push_back(node + nFunc_x);
      right.push_back(node + nFunc_x * (nFunc_y - 2));
    }
  }

  for( int yy=1; yy<nFunc_y-1; ++yy )
  {
    for(int xx=1; xx<nFunc_x-1; ++xx )
    {
      node = yy * nFunc_x + xx;
      bottom.push_back(node + nFxy);
      top.push_back(node + (nFunc_z-2)*nFxy);
    }
  }

  VEC_T::shrink2fit(front);  VEC_T::shrink2fit(back);
  VEC_T::shrink2fit(left);   VEC_T::shrink2fit(right);
  VEC_T::shrink2fit(top);    VEC_T::shrink2fit(bottom);

  // Now sepcify the interior edges
  for(int xx=1; xx<nFunc_x-1; ++xx)
  {
    edge01.push_back(xx + nFxy + nFunc_x);
    edge23.push_back(xx + nFxy + (nFunc_y-2) * nFunc_x );
    edge45.push_back(xx + (nFunc_z-2)*nFxy + nFunc_x);
    edge67.push_back(xx + (nFunc_z-2)*nFxy + (nFunc_y-2)*nFunc_x );
  }
  VEC_T::shrink2fit(edge01); VEC_T::shrink2fit(edge23);
  VEC_T::shrink2fit(edge45); VEC_T::shrink2fit(edge67);

  for(int yy=1; yy<nFunc_y-1; ++yy)
  {
    edge02.push_back(yy*nFunc_x + nFxy + 1);
    edge13.push_back(yy*nFunc_x + nFxy + nFunc_x -2);
    edge46.push_back(yy*nFunc_x + (nFunc_z-2)*nFxy + 1);
    edge57.push_back(yy*nFunc_x + (nFunc_z-2)*nFxy + nFunc_x - 2);
  }
  VEC_T::shrink2fit(edge02); VEC_T::shrink2fit(edge13);
  VEC_T::shrink2fit(edge46); VEC_T::shrink2fit(edge57);

  for(int zz=1; zz<nFunc_z-1; ++zz)
  {
    edge04.push_back(zz*nFxy + nFunc_x + 1);
    edge15.push_back(zz*nFxy + nFunc_x + nFunc_x - 2);
    edge26.push_back(zz*nFxy + (nFunc_y-2)*nFunc_x + 1);
    edge37.push_back(zz*nFxy + (nFunc_y-2)*nFunc_x + nFunc_x - 2);
  }
  VEC_T::shrink2fit(edge04); VEC_T::shrink2fit(edge15);
  VEC_T::shrink2fit(edge26); VEC_T::shrink2fit(edge37);

  // corner nodes
  corner.push_back(nFxy + nFunc_x + 1);
  corner.push_back(nFxy + nFunc_x + nFunc_x - 2);
  corner.push_back(nFxy + (nFunc_y-2)*nFunc_x + 1);
  corner.push_back(nFxy + (nFunc_y-2)*nFunc_x + nFunc_x - 2);

  corner.push_back((nFunc_z-2)*nFxy + nFunc_x + 1);
  corner.push_back((nFunc_z-2)*nFxy + nFunc_x + nFunc_x - 2);
  corner.push_back((nFunc_z-2)*nFxy + (nFunc_y-2)*nFunc_x + 1);
  corner.push_back((nFunc_z-2)*nFxy + (nFunc_y-2)*nFunc_x + nFunc_x - 2);

  VEC_T::shrink2fit(corner);
}


void BoundaryCond::Generate_BCNodes_D( const int &nFunc_x, const int &nFunc_y,
    const int &nFunc_z, const int &f_start,
    int * &back_front, int &nbf, 
    int * &left_right, int &nlr,
    int * &bottom_top, int &nbt,
    int * &edge01, int &nex, 
    int * &edge02, int &ney,
    int * &edge04, int &nez,
    int * &corner ) const
{
  if(nFunc_z <= 2)
  {
    std::cerr<<"Error: in Generate_BCNodes_D, nFunc_z has to be greater than 2. \n";
    exit(EXIT_FAILURE);
  }
  if(nFunc_y <= 2)
  {
    std::cerr<<"Error: in Generate_BCNodes_D, nFunc_y has to be greater than 2. \n";
    exit(EXIT_FAILURE);
  }
  if(nFunc_x <= 2)
  {
    std::cerr<<"Error: in Generate_BCNodes_D, nFunc_x has to be greater than 2. \n";
    exit(EXIT_FAILURE);
  }

  nbf = (nFunc_z-2) * (nFunc_y-2);
  nlr = (nFunc_z-2) * (nFunc_x-2);
  nbt = (nFunc_y-2) * (nFunc_x-2); 
  
  back_front   = new int [2*nbf];
  left_right   = new int [2*nlr];
  bottom_top   = new int [2*nbt];

  nex = nFunc_x - 2;
  ney = nFunc_y - 2;
  nez = nFunc_z - 2;

  edge01 = new int [4*nex]; 
  edge02 = new int [4*ney];
  edge04 = new int [4*nez];
  corner = new int [8];
 
  const int nfxm1     = nFunc_x - 1; 
  const int nfxy      = nFunc_x * nFunc_y;
  const int nfxfym1   = nFunc_x * (nFunc_y - 1);
  const int nfzm1fxfy = (nFunc_z - 1) * nFunc_x * nFunc_y;

  int node, temp;
  int index_bf = 0;
  int index_lr = 0;
  for( int zz=1; zz<nFunc_z-1; ++zz )
  {
    temp = f_start + zz * nfxy;
    for( int yy=1; yy<nFunc_y-1; ++yy )
    {
      node = temp + yy * nFunc_x;
      back_front[index_bf] = node;
      back_front[index_bf + nbf] = node + nfxm1;
      ++index_bf;
    }

    for( int xx=1; xx<nFunc_x-1; ++xx )
    {
      left_right[index_lr] = temp + xx;
      left_right[index_lr+nlr] = temp + xx + nfxfym1; 
      ++index_lr;
    }
  }

  int index_bt = 0;
  for(int yy=1; yy<nFunc_y-1; ++yy)
  {
    temp = f_start + yy * nFunc_x;
    for(int xx=1; xx<nFunc_x-1; ++xx)
    {
      bottom_top[index_bt] = temp + xx;
      bottom_top[index_bt + nbt] = temp + xx + nfzm1fxfy;
      ++index_bt;
    }
  }

  int counter = 0;
  const int m23 = (nFunc_y-1) * nFunc_x;
  const int m45 = (nFunc_z-1) * nFunc_x * nFunc_y;
  const int m67 = m45 + (nFunc_y-1) * nFunc_x;
  for(int xx=1; xx<nFunc_x-1; ++xx)
  {
    edge01[counter] = xx + f_start;
    edge01[counter + nex] = xx + m23 + f_start;
    edge01[counter + 2*nex] = xx + m45 + f_start;
    edge01[counter + 3*nex] = xx + m67 + f_start; 
    ++counter;
  }

  counter = 0;
  int m02;
  const int m13 = nFunc_x - 1;
  const int m46 = (nFunc_z-1)*nFunc_x*nFunc_y;
  const int m57 = m46 + m13;
  for(int yy=1; yy<nFunc_y-1; ++yy)
  {
    m02 = yy * nFunc_x + f_start;
    edge02[counter] = m02;
    edge02[counter + ney] = m02 + m13;
    edge02[counter + 2*ney] = m02 + m46;
    edge02[counter + 3*ney] = m02 + m57;
    ++counter;
  }

  counter = 0;
  int m04;
  const int m15 = nFunc_x - 1;
  const int m26 = (nFunc_y-1) * nFunc_x;
  const int m37 = m26 + m15;
  for(int zz=1; zz<nFunc_z-1; ++zz)
  {
    m04 = zz * nfxy + f_start;
    edge04[counter] = m04;
    edge04[counter + nez] = m04 + m15;
    edge04[counter + 2*nez] = m04 + m26;
    edge04[counter + 3*nez] = m04 + m37;
    ++counter;
  }

  corner[0] = f_start;
  corner[1] = f_start + nFunc_x - 1;
  corner[2] = f_start + (nFunc_y-1) * nFunc_x;
  corner[3] = f_start + (nFunc_y-1)*nFunc_x + nFunc_x - 1;
  corner[4] = f_start + (nFunc_z-1)*nFunc_x*nFunc_y;
  corner[5] = f_start + (nFunc_z-1)*nFunc_x*nFunc_y + nFunc_x -1;
  corner[6] = f_start + (nFunc_z-1)*nFunc_x*nFunc_y + (nFunc_y-1)*nFunc_x;
  corner[7] = f_start + (nFunc_z-1)*nFunc_x*nFunc_y + (nFunc_y-1)*nFunc_x + nFunc_x - 1;
}


void BoundaryCond::Generate_BCElems_A( const IMesh * const &mesh,
    std::vector<int> &front, std::vector<int> &back, std::vector<int> &left,
    std::vector<int> &right, std::vector<int> &top, std::vector<int> &bottom ) const
{
  // Reads single patch NURBS information from mesh, then determine the 
  // elements that have faces on the front, back, ..., bottom boundary faces.
  front.clear(); back.clear(); left.clear(); right.clear();
  top.clear(); bottom.clear();

  const int nElem_x = mesh->get_nElem_x();
  const int nElem_y = mesh->get_nElem_y();
  const int nElem_z = mesh->get_nElem_z();

  for(int zz=0; zz<nElem_z; ++zz)
  {
    for(int yy=0; yy<nElem_y; ++yy)
    {
      int elem = zz * nElem_x * nElem_y + yy * nElem_x;
      back.push_back(elem);
      front.push_back(elem + nElem_x - 1);
    }
  }
  VEC_T::shrink2fit(front); VEC_T::shrink2fit(back);

  for(int zz=0; zz<nElem_z; ++zz)
  {
    for(int xx=0; xx<nElem_x; ++xx)
    {
      int elem = zz * nElem_x * nElem_y + xx;
      left.push_back(elem);
      right.push_back(elem + nElem_x * (nElem_y-1));
    }
  }
  VEC_T::shrink2fit(left); VEC_T::shrink2fit(right);

  for(int yy=0; yy<nElem_y; ++yy)
  {
    for(int xx=0 ;xx<nElem_x; ++xx)
    {
      int elem = yy*nElem_x + xx;
      bottom.push_back(elem);
      top.push_back(elem+(nElem_z-1)*nElem_x*nElem_y);
    }
  }
  VEC_T::shrink2fit(bottom); VEC_T::shrink2fit(top);
}



void BoundaryCond::Generate_BCElems_B( const IMesh * const &mesh,
    std::vector<int> &front, std::vector<int> &back, std::vector<int> &left,
    std::vector<int> &right, std::vector<int> &top, std::vector<int> &bottom ) const
{
  front.clear(); back.clear(); left.clear(); right.clear();
  top.clear(); bottom.clear();

  const int nElem_x = mesh->get_nElem_x();
  const int nElem_y = mesh->get_nElem_y();
  const int nElem_z = mesh->get_nElem_z();

  const int es = mesh->get_nElem_start();
  
  const int nExEy = nElem_x * nElem_y;

  int temp;

  for(int zz=0; zz<nElem_z; ++zz)
  {
    temp = zz * nExEy + es;
    for(int yy=0; yy<nElem_y; ++yy)
    {
      int elem = temp + yy * nElem_x;
      back.push_back(elem);
      front.push_back(elem + nElem_x - 1);
    }
  }
  VEC_T::shrink2fit(front); VEC_T::shrink2fit(back);

  for(int zz=0; zz<nElem_z; ++zz)
  {
    temp = zz * nExEy + es;
    for(int xx=0; xx<nElem_x; ++xx)
    {
      int elem = temp + xx;
      left.push_back(elem);
      right.push_back(elem + nElem_x * (nElem_y-1));
    }
  }
  VEC_T::shrink2fit(left); VEC_T::shrink2fit(right);

  for(int yy=0; yy<nElem_y; ++yy)
  {
    temp = yy * nElem_x + es;
    for(int xx=0 ;xx<nElem_x; ++xx)
    {
      int elem = temp + xx;
      bottom.push_back(elem);
      top.push_back(elem+(nElem_z-1)*nExEy);
    }
  }
  VEC_T::shrink2fit(bottom); VEC_T::shrink2fit(top);
}


void BoundaryCond::BC_test_1() const
{
  const int insdeg  = 1;
  const int intdeg  = 2;
  const int inudeg  = 3;
  const int inElemx = 1200 - insdeg;
  const int inElemy = 1100 - intdeg;
  const int inElemz = 1000 - inudeg;
  IMesh * testmesh = new Mesh_NURBS_1Patch_3D(insdeg, intdeg, inudeg, 
      0.2, 0.3, 0.5, 0.1, 0.2, 0.4, 
      inElemx, inElemy, inElemz);
  
  std::cout<<"Test mesh has "<<testmesh->get_nFunc_x()<<" * "<<testmesh->get_nFunc_y();
  std::cout<<" * "<<testmesh->get_nFunc_z()<<" = "<<testmesh->get_nFunc()<<" basis functions. \n";
  
  clock_t log_time = clock();
  std::vector<int> bfr, bba, ble, bri, bto, bbo;
  std::vector<int> be01, be02, be13, be23, be45, be46, be57, be67, be15, be37, be04, be26, bcorner;
  
  Generate_BCNodes_B(testmesh, bfr, bba, ble, bri, bto, bbo, be01, be02, be13,
      be23, be45, be46, be57, be67, be15, be37, be04, be26, bcorner);

  log_time = clock() - log_time;
  std::cout<<"Generate_BCNodes_B function takes ";
  std::cout<<((double) log_time)/CLOCKS_PER_SEC<<" seconds"<<std::endl;

  
  log_time = clock();
  const int fun_start = 2;

  int * dbf; int dnbf;
  int * dlr; int dnlr;
  int * dbt; int dnbt;
  int * de01; int dnex;
  int * de02; int dney;
  int * de04; int dnez;
  int * dcorner;

  Generate_BCNodes_D(testmesh->get_nFunc_x(), testmesh->get_nFunc_y(),
      testmesh->get_nFunc_z(), fun_start,
      dbf, dnbf, dlr, dnlr, dbt, dnbt, de01, dnex, de02, dney, de04, dnez, dcorner );
  
  log_time = clock() - log_time;
  std::cout<<"Generate_BCNodes_D function takes ";
  std::cout<<((double) log_time)/CLOCKS_PER_SEC<<" seconds"<<std::endl;
  

  // Start comparing the two function results
  assert(dnbf == int(bba.size())); assert(dnbf == int(bfr.size()));
  for(int ii=0; ii<dnbf; ++ii)
  {
    assert(dbf[ii] = bba[ii] + fun_start);
    assert(dbf[ii+dnbf] = bfr[ii] + fun_start);
  }

  assert(dnlr == int(ble.size())); assert(dnlr == int(bri.size()));
  for(int ii=0; ii<dnlr; ++ii)
  {
    assert(dlr[ii] == ble[ii] + fun_start);
    assert(dlr[dnlr+ii] == bri[ii] + fun_start);
  }

  assert(dnbt == int(bto.size())); assert(dnbt == int(bbo.size()));
  for(int ii=0; ii<dnbt; ++ii)
  {
    assert(dbt[ii] == bbo[ii] + fun_start);
    assert(dbt[dnbt+ii] = bto[ii] + fun_start);
  }

  assert(dnex == int(be01.size())); assert(dnex == int(be23.size()));
  assert(dnex == int(be45.size())); assert(dnex == int(be67.size()));
  for(int ii=0; ii<dnex; ++ii)
  {
    assert(de01[ii] = be01[ii] + fun_start);
    assert(de01[ii+dnex] = be23[ii] + fun_start);
    assert(de01[ii+2*dnex] = be45[ii] + fun_start);
    assert(de01[ii+3*dnex] = be67[ii] + fun_start);
  }

  assert(dney == int(be02.size())); assert(dney == int(be13.size()));
  assert(dney == int(be46.size())); assert(dney == int(be57.size()));
  for(int ii=0; ii<dney; ++ii)
  {
    assert(de02[ii] = be02[ii] + fun_start);
    assert(de02[ii+dney] = be13[ii] + fun_start);
    assert(de02[ii+2*dney] = be46[ii] + fun_start);
    assert(de02[ii+3*dney] = be57[ii] + fun_start);
  }

  assert(dnez == int(be04.size())); assert(dnez == int(be15.size()));
  assert(dnez == int(be26.size())); assert(dnez == int(be37.size()));
  for(int ii=0; ii<dnez; ++ii)
  {
    assert(de04[ii] = be04[ii] + fun_start);
    assert(de04[ii+dnez] = be15[ii] + fun_start);
    assert(de04[ii+2*dnez] = be26[ii] + fun_start);
    assert(de04[ii+3*dnez] = be37[ii] + fun_start);
  }
  
  for(int ii=0; ii<8; ++ii)
    assert(dcorner[ii] = bcorner[ii] + fun_start);


  delete [] dbf; delete [] dlr; delete [] dbt;
  delete [] de01; delete [] de02; delete [] de04; delete [] dcorner;
  delete testmesh;

  std::cout<<"-- Test 1 passed. The Generate_BCNodes_D gives identical results to _B version. \n";
}

void BoundaryCond::BC_test_2() const
{
  const int insdeg  = 1;
  const int intdeg  = 2;
  const int inudeg  = 3;
  const int inElemx = 1232 - insdeg;
  const int inElemy = 1152 - intdeg;
  const int inElemz = 125 - inudeg;
  IMesh * atestmesh = new Mesh_NURBS_1Patch_3D(insdeg, intdeg, inudeg, 
      0.2, 0.3, 0.5, 0.1, 0.2, 0.4, 
      inElemx, inElemy, inElemz);

  std::vector<int> afro, abac, alef, arig, atop, abot;

  Generate_BCElems_A(atestmesh, afro, abac, alef, arig, atop, abot);

  const int patch_index = 1;
  const int patch_elem  = 102;
  const int patch_func  = 152;
  IMesh * btestmesh = new Mesh_NURBS_1ofnPatch_3D(insdeg, intdeg, inudeg,
      0.2, 0.3, 0.5, 0.1, 0.2, 0.4, inElemx, inElemy, inElemz,
      patch_index, patch_elem, patch_func );

  std::vector<int> bfro, bbac, blef, brig, btop, bbot;

  Generate_BCElems_B(btestmesh, bfro, bbac, blef, brig, btop, bbot);

  assert(afro.size() == bfro.size());
  assert(abac.size() == bbac.size());
  assert(alef.size() == blef.size());
  assert(arig.size() == brig.size());
  assert(atop.size() == btop.size());
  assert(abot.size() == bbot.size());

  for(unsigned int ii=0; ii<afro.size(); ++ii)
    assert(afro[ii] == bfro[ii] - patch_elem);

  for(unsigned int ii=0; ii<abac.size(); ++ii)
    assert(abac[ii] == bbac[ii] - patch_elem);

  for(unsigned int ii=0; ii<alef.size(); ++ii)
    assert(alef[ii] == blef[ii] - patch_elem);

  for(unsigned int ii=0; ii<arig.size(); ++ii)
    assert(arig[ii] == brig[ii] - patch_elem);

  for(unsigned int ii=0; ii<atop.size(); ++ii)
    assert(atop[ii] == btop[ii] - patch_elem);

  for(unsigned int ii=0; ii<abot.size(); ++ii)
    assert(abot[ii] == bbot[ii] - patch_elem);

  std::cout<<"--- Test 2 passed. The Generate_BCElems_B gives identical results to _A version. \n";

  delete atestmesh;
  delete btestmesh;
}


void BoundaryCond::BC_type_1( const IMesh * const &mesh )
{
  // all nodes are dirichlet
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  Generate_BCNodes_A(mesh, nodes_front, nodes_back, nodes_left,
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


void BoundaryCond::BC_type_2( const IMesh * const &mesh )
{
  std::cout<<"Warning: This Periodic boundary imposition function is slow. Use BC_type_8. \n";

  // full periodic. master - slave: front - back, left - right, bottom - top
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  std::vector<int> nodes_edge01, nodes_edge02, nodes_edge13, nodes_edge23,
    nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67, nodes_edge15,
    nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner;

  Generate_BCNodes_B(mesh, nodes_front, nodes_back, nodes_left, nodes_right,
      nodes_top, nodes_bottom, nodes_edge01, nodes_edge02, nodes_edge13,
      nodes_edge23, nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67,
      nodes_edge15, nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner);  

  // 1. set face master slave relations
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_back.begin(), nodes_back.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_right.begin(), nodes_right.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_top.begin(), nodes_top.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_front.begin(), nodes_front.end());
  per_master_nodes.insert(per_master_nodes.end(), nodes_left.begin(), nodes_left.end());
  per_master_nodes.insert(per_master_nodes.end(), nodes_bottom.begin(), nodes_bottom.end());

  // 2 set edge master slave relations
  // 2.1 master edge01, slave edge23 edge67 edge45
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge23.begin(), nodes_edge23.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge67.begin(), nodes_edge67.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge45.begin(), nodes_edge45.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge01.begin(), nodes_edge01.end());
  per_master_nodes.insert(per_master_nodes.end(), nodes_edge01.begin(), nodes_edge01.end());
  per_master_nodes.insert(per_master_nodes.end(), nodes_edge01.begin(), nodes_edge01.end());

  // 2.2 master edge15 slave edge37 edge26 edge04
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge37.begin(), nodes_edge37.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge26.begin(), nodes_edge26.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge04.begin(), nodes_edge04.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge15.begin(), nodes_edge15.end());
  per_master_nodes.insert(per_master_nodes.end(), nodes_edge15.begin(), nodes_edge15.end());
  per_master_nodes.insert(per_master_nodes.end(), nodes_edge15.begin(), nodes_edge15.end());

  // 2.3 master edge13 slave edge02 edge46 edge57
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge02.begin(), nodes_edge02.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge46.begin(), nodes_edge46.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge57.begin(), nodes_edge57.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge13.begin(), nodes_edge13.end());
  per_master_nodes.insert(per_master_nodes.end(), nodes_edge13.begin(), nodes_edge13.end());
  per_master_nodes.insert(per_master_nodes.end(), nodes_edge13.begin(), nodes_edge13.end());

  // 3 set corner nodes master slave relations
  per_master_nodes.push_back(nodes_corner[1]);
  per_slave_nodes.push_back(nodes_corner[0]);
  for(unsigned int ii=2; ii<nodes_corner.size(); ++ii)
  {
    per_master_nodes.push_back(nodes_corner[1]);
    per_slave_nodes.push_back(nodes_corner[ii]);
  }

  num_per_nodes = per_master_nodes.size();

  assert(num_per_nodes == per_slave_nodes.size()); 
}

void BoundaryCond::BC_type_3( const IMesh * const &mesh )
{ 
  std::cout<<"-----> Nothing on the boundary is enforced. \n"; 
}

void BoundaryCond::BC_type_4( const IMesh * const &mesh )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  Generate_BCNodes_A(mesh, nodes_front, nodes_back, nodes_left,
      nodes_right, nodes_top, nodes_bottom );

  // front master, back slave
  per_master_nodes.insert(per_master_nodes.end(), nodes_front.begin(), nodes_front.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_back.begin(), nodes_back.end());

  num_per_nodes = per_master_nodes.size();
  assert(num_per_nodes == per_slave_nodes.size());

  // rest are dirichlet nodes
  dir_nodes.insert(dir_nodes.end(), nodes_left.begin(), nodes_left.end());
  dir_nodes.insert(dir_nodes.end(), nodes_right.begin(), nodes_right.end());
  dir_nodes.insert(dir_nodes.end(), nodes_top.begin(), nodes_top.end());
  dir_nodes.insert(dir_nodes.end(), nodes_bottom.begin(), nodes_bottom.end());

  num_dir_nodes = dir_nodes.size();
}

void BoundaryCond::BC_type_5( const IMesh * const &mesh )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  Generate_BCNodes_A(mesh, nodes_front, nodes_back, nodes_left,
      nodes_right, nodes_top, nodes_bottom );

  // front master, back slave
  per_master_nodes.insert(per_master_nodes.end(), nodes_front.begin(), nodes_front.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_back.begin(), nodes_back.end());

  num_per_nodes = per_master_nodes.size();
  assert(num_per_nodes == per_slave_nodes.size());
}

void BoundaryCond::BC_type_6( const IMesh * const &mesh )
{
  Generate_BCElems_A(mesh, front_elem, back_elem, left_elem, right_elem,
      top_elem, bottom_elem);
  num_front_elem = front_elem.size();
  num_back_elem = back_elem.size();
  num_left_elem = left_elem.size();
  num_right_elem = right_elem.size();
  num_top_elem = top_elem.size();
  num_bottom_elem = bottom_elem.size();
}

void BoundaryCond::BC_type_7( const IMesh * const &mesh )
{
  Generate_BCElems_A(mesh, front_elem, back_elem, left_elem, right_elem,
      top_elem, bottom_elem);

  front_elem.clear();
  back_elem.clear();
  left_elem.clear();
  right_elem.clear();
  top_elem.clear();

  num_front_elem = front_elem.size();
  num_back_elem = back_elem.size();
  num_left_elem = left_elem.size();
  num_right_elem = right_elem.size();
  num_top_elem = top_elem.size();
  num_bottom_elem = bottom_elem.size();
}

void BoundaryCond::BC_type_8( const IMesh * const &mesh )
{
  // full periodic. master - slave: front - back, left - right, bottom - top
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  std::vector<int> nodes_edge01, nodes_edge02, nodes_edge13, nodes_edge23,
    nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67, nodes_edge15,
    nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner;

  Generate_BCNodes_B(mesh, nodes_front, nodes_back, nodes_left, nodes_right,
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
  std::cout<<"-----> full periodic. \n";
}


void BoundaryCond::BC_type_9( const IMesh * const &mesh )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  std::vector<int> nodes_edge01, nodes_edge02, nodes_edge13, nodes_edge23,
    nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67, nodes_edge15,
    nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner;

  Generate_BCNodes_B(mesh, nodes_front, nodes_back, nodes_left, nodes_right,
      nodes_top, nodes_bottom, nodes_edge01, nodes_edge02, nodes_edge13,
      nodes_edge23, nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67,
      nodes_edge15, nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner);  

  // front master, back slave
  per_master_nodes.insert(per_master_nodes.end(), nodes_front.begin(), nodes_front.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_back.begin(), nodes_back.end());

  num_per_nodes = per_master_nodes.size();
  assert(num_per_nodes == per_slave_nodes.size());

  std::cout<<"-----> front-back periodic. \n";
}


void BoundaryCond::BC_type_10( const IMesh * const &mesh )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  std::vector<int> nodes_edge01, nodes_edge02, nodes_edge13, nodes_edge23,
    nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67, nodes_edge15,
    nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner;

  Generate_BCNodes_B(mesh, nodes_front, nodes_back, nodes_left, nodes_right,
      nodes_top, nodes_bottom, nodes_edge01, nodes_edge02, nodes_edge13,
      nodes_edge23, nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67,
      nodes_edge15, nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner);  

  // front master, back slave
  per_master_nodes.insert(per_master_nodes.end(), nodes_front.begin(), nodes_front.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_back.begin(), nodes_back.end());

  num_per_nodes = per_master_nodes.size();
  assert(num_per_nodes == per_slave_nodes.size());

  // rest are dirichlet nodes
  dir_nodes.insert(dir_nodes.end(), nodes_left.begin(), nodes_left.end());
  dir_nodes.insert(dir_nodes.end(), nodes_right.begin(), nodes_right.end());
  dir_nodes.insert(dir_nodes.end(), nodes_top.begin(), nodes_top.end());
  dir_nodes.insert(dir_nodes.end(), nodes_bottom.begin(), nodes_bottom.end());

  dir_nodes.insert(dir_nodes.end(), nodes_edge01.begin(), nodes_edge01.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge02.begin(), nodes_edge02.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge13.begin(), nodes_edge13.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge23.begin(), nodes_edge23.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge45.begin(), nodes_edge45.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge46.begin(), nodes_edge46.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge57.begin(), nodes_edge57.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge67.begin(), nodes_edge67.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge15.begin(), nodes_edge15.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge37.begin(), nodes_edge37.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge04.begin(), nodes_edge04.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge26.begin(), nodes_edge26.end());
  dir_nodes.insert(dir_nodes.end(), nodes_corner.begin(), nodes_corner.end());

  num_dir_nodes = dir_nodes.size();
  std::cout<<"-----> front-back periodic, rest dirichlet \n";
}

void BoundaryCond::BC_type_11( const IMesh * const &mesh )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  std::vector<int> nodes_edge01, nodes_edge02, nodes_edge13, nodes_edge23,
    nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67, nodes_edge15,
    nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner;

  Generate_BCNodes_B(mesh, nodes_front, nodes_back, nodes_left, nodes_right,
      nodes_top, nodes_bottom, nodes_edge01, nodes_edge02, nodes_edge13,
      nodes_edge23, nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67,
      nodes_edge15, nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner);  

  // front master, back slave
  per_master_nodes.insert(per_master_nodes.end(), nodes_front.begin(), nodes_front.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_back.begin(), nodes_back.end());

  // bottom master, top slave
  per_master_nodes.insert(per_master_nodes.end(), nodes_bottom.begin(), nodes_bottom.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_top.begin(), nodes_top.end());

  // edge 13 master, 57 46 02 slave
  per_master_nodes.insert(per_master_nodes.end(), nodes_edge13.begin(), nodes_edge13.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge57.begin(), nodes_edge57.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge13.begin(), nodes_edge13.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge46.begin(), nodes_edge46.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge13.begin(), nodes_edge13.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge02.begin(), nodes_edge02.end());

  num_per_nodes = per_master_nodes.size();
  assert(num_per_nodes == per_slave_nodes.size());

  std::cout<<"-----> front-back, top-bottom periodic. \n";
}

void BoundaryCond::BC_type_12( const IMesh * const &mesh )
{
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  std::vector<int> nodes_edge01, nodes_edge02, nodes_edge13, nodes_edge23,
    nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67, nodes_edge15,
    nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner;

  Generate_BCNodes_B(mesh, nodes_front, nodes_back, nodes_left, nodes_right,
      nodes_top, nodes_bottom, nodes_edge01, nodes_edge02, nodes_edge13,
      nodes_edge23, nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67,
      nodes_edge15, nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner);  

  // front master, back slave
  per_master_nodes.insert(per_master_nodes.end(), nodes_front.begin(), nodes_front.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_back.begin(), nodes_back.end());

  // bottom master, top slave
  per_master_nodes.insert(per_master_nodes.end(), nodes_bottom.begin(), nodes_bottom.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_top.begin(), nodes_top.end());

  // edge 13 master, 57 46 02 slave
  per_master_nodes.insert(per_master_nodes.end(), nodes_edge13.begin(), nodes_edge13.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge57.begin(), nodes_edge57.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge13.begin(), nodes_edge13.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge46.begin(), nodes_edge46.end());

  per_master_nodes.insert(per_master_nodes.end(), nodes_edge13.begin(), nodes_edge13.end());
  per_slave_nodes.insert(per_slave_nodes.end(), nodes_edge02.begin(), nodes_edge02.end());

  num_per_nodes = per_master_nodes.size();
  assert(num_per_nodes == per_slave_nodes.size());

  // rest are dirichlet nodes
  dir_nodes.insert(dir_nodes.end(), nodes_left.begin(), nodes_left.end());
  dir_nodes.insert(dir_nodes.end(), nodes_right.begin(), nodes_right.end());

  dir_nodes.insert(dir_nodes.end(), nodes_edge01.begin(), nodes_edge01.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge23.begin(), nodes_edge23.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge45.begin(), nodes_edge45.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge67.begin(), nodes_edge67.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge15.begin(), nodes_edge15.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge37.begin(), nodes_edge37.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge04.begin(), nodes_edge04.end());
  dir_nodes.insert(dir_nodes.end(), nodes_edge26.begin(), nodes_edge26.end());
  dir_nodes.insert(dir_nodes.end(), nodes_corner.begin(), nodes_corner.end());

  num_dir_nodes = dir_nodes.size();
  std::cout<<"-----> front-back, top-bottom periodic, rest dirichlet \n";
}


void BoundaryCond::BC_type_13( const IMesh * const &mesh )
{
  // Nodes for outside faces
  std::vector<int> n_front, n_back, n_left, n_right, n_top, n_bottom;
  std::vector<int> n_edge01, n_edge02, n_edge13, n_edge23,
    n_edge45, n_edge46, n_edge57, n_edge67, n_edge15,
    n_edge37, n_edge04, n_edge26, n_corner;

  Generate_BCNodes_B(mesh, n_front, n_back, n_left, n_right,
      n_top, n_bottom, n_edge01, n_edge02, n_edge13,
      n_edge23, n_edge45, n_edge46, n_edge57, n_edge67,
      n_edge15, n_edge37, n_edge04, n_edge26, n_corner);  

  // Nodes for interior cube
  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;
  std::vector<int> nodes_edge01, nodes_edge02, nodes_edge13, nodes_edge23,
    nodes_edge45, nodes_edge46, nodes_edge57, nodes_edge67, nodes_edge15,
    nodes_edge37, nodes_edge04, nodes_edge26, nodes_corner;

  Generate_BCNodes_C(mesh, nodes_front, nodes_back, nodes_left, nodes_right,
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


void BoundaryCond::BC_type_14( const IMesh * const &mesh )
{
  std::vector<int> n_front, n_back, n_left, n_right, n_top, n_bottom;
  std::vector<int> n_edge01, n_edge02, n_edge13, n_edge23,
    n_edge45, n_edge46, n_edge57, n_edge67, n_edge15,
    n_edge37, n_edge04, n_edge26, n_corner;

  Generate_BCNodes_B(mesh, n_front, n_back, n_left, n_right,
      n_top, n_bottom, n_edge01, n_edge02, n_edge13,
      n_edge23, n_edge45, n_edge46, n_edge57, n_edge67,
      n_edge15, n_edge37, n_edge04, n_edge26, n_corner);  


  dir_nodes.insert(dir_nodes.end(), n_front.begin(), n_front.end());
  dir_nodes.insert(dir_nodes.end(), n_back.begin(), n_back.end());

  dir_nodes.insert(dir_nodes.end(), n_edge13.begin(), n_edge13.end());
  dir_nodes.insert(dir_nodes.end(), n_edge15.begin(), n_edge15.end());
  dir_nodes.insert(dir_nodes.end(), n_edge37.begin(), n_edge37.end());
  dir_nodes.insert(dir_nodes.end(), n_edge57.begin(), n_edge57.end());

  dir_nodes.insert(dir_nodes.end(), n_edge02.begin(), n_edge02.end());
  dir_nodes.insert(dir_nodes.end(), n_edge04.begin(), n_edge04.end());
  dir_nodes.insert(dir_nodes.end(), n_edge26.begin(), n_edge26.end());
  dir_nodes.insert(dir_nodes.end(), n_edge46.begin(), n_edge46.end());

  dir_nodes.insert(dir_nodes.end(), n_corner.begin(), n_corner.end());

  num_dir_nodes = dir_nodes.size();
  num_per_nodes = 0;
  std::cout<<"-----> Slip for velocity in x direction. \n";
}

void BoundaryCond::BC_type_15( const IMesh * const &mesh )
{
  std::vector<int> n_front, n_back, n_left, n_right, n_top, n_bottom;
  std::vector<int> n_edge01, n_edge02, n_edge13, n_edge23,
    n_edge45, n_edge46, n_edge57, n_edge67, n_edge15,
    n_edge37, n_edge04, n_edge26, n_corner;

  Generate_BCNodes_B(mesh, n_front, n_back, n_left, n_right,
      n_top, n_bottom, n_edge01, n_edge02, n_edge13,
      n_edge23, n_edge45, n_edge46, n_edge57, n_edge67,
      n_edge15, n_edge37, n_edge04, n_edge26, n_corner);  


  dir_nodes.insert(dir_nodes.end(), n_left.begin(), n_left.end());
  dir_nodes.insert(dir_nodes.end(), n_right.begin(), n_right.end());

  dir_nodes.insert(dir_nodes.end(), n_edge01.begin(), n_edge01.end());
  dir_nodes.insert(dir_nodes.end(), n_edge04.begin(), n_edge04.end());
  dir_nodes.insert(dir_nodes.end(), n_edge15.begin(), n_edge15.end());
  dir_nodes.insert(dir_nodes.end(), n_edge45.begin(), n_edge45.end());

  dir_nodes.insert(dir_nodes.end(), n_edge26.begin(), n_edge26.end());
  dir_nodes.insert(dir_nodes.end(), n_edge23.begin(), n_edge23.end());
  dir_nodes.insert(dir_nodes.end(), n_edge37.begin(), n_edge37.end());
  dir_nodes.insert(dir_nodes.end(), n_edge67.begin(), n_edge67.end());

  dir_nodes.insert(dir_nodes.end(), n_corner.begin(), n_corner.end());

  num_dir_nodes = dir_nodes.size();
  num_per_nodes = 0;
  std::cout<<"-----> Slip for velocity in y direction. \n";
}


void BoundaryCond::BC_type_16( const IMesh * const &mesh )
{
  std::vector<int> n_front, n_back, n_left, n_right, n_top, n_bottom;
  std::vector<int> n_edge01, n_edge02, n_edge13, n_edge23,
    n_edge45, n_edge46, n_edge57, n_edge67, n_edge15,
    n_edge37, n_edge04, n_edge26, n_corner;

  Generate_BCNodes_B(mesh, n_front, n_back, n_left, n_right,
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
  std::cout<<"-----> Slip for velocity in z direction. \n";
}


void BoundaryCond::BC_type_17( const IMesh * const &mesh )
{
  std::vector<int> n_front, n_back, n_left, n_right, n_top, n_bottom;
  std::vector<int> n_edge01, n_edge02, n_edge13, n_edge23,
    n_edge45, n_edge46, n_edge57, n_edge67, n_edge15,
    n_edge37, n_edge04, n_edge26, n_corner;

  Generate_BCNodes_B(mesh, n_front, n_back, n_left, n_right,
      n_top, n_bottom, n_edge01, n_edge02, n_edge13,
      n_edge23, n_edge45, n_edge46, n_edge57, n_edge67,
      n_edge15, n_edge37, n_edge04, n_edge26, n_corner);  

  dir_nodes.insert(dir_nodes.end(), n_corner[0]);

  num_dir_nodes = 1;
  num_per_nodes = 0;

  std::cout<<"-----> Pressure fixed at (0,0,0). \n";
}


void BoundaryCond::BC_type_18( const IMesh * const &mesh )
{
  std::vector<int> n_front, n_back, n_left, n_right, n_top, n_bottom;

  Generate_BCNodes_A( mesh, n_front, n_back, n_left, n_right,
      n_top, n_bottom );

  dir_nodes.insert(dir_nodes.end(), n_top.begin(), n_top.end());
  dir_nodes.insert(dir_nodes.end(), n_left.begin(), n_left.end());
  dir_nodes.insert(dir_nodes.end(), n_right.begin(), n_right.end());
  dir_nodes.insert(dir_nodes.end(), n_front.begin(), n_front.end());
  dir_nodes.insert(dir_nodes.end(), n_back.begin(), n_back.end());

  num_dir_nodes = dir_nodes.size();
  num_per_nodes = 0;

  std::cout<<"-----> Dirichlet on top, left, right, front, back. Do nothing on bottom. \n";
}


void BoundaryCond::BC_type_19( const IMesh * const &mesh )
{
  std::vector<int> n_front, n_back, n_left, n_right, n_top, n_bottom;
  std::vector<int> n_edge01, n_edge02, n_edge13, n_edge23,
    n_edge45, n_edge46, n_edge57, n_edge67, n_edge15,
    n_edge37, n_edge04, n_edge26, n_corner;

  Generate_BCNodes_B(mesh, n_front, n_back, n_left, n_right,
      n_top, n_bottom, n_edge01, n_edge02, n_edge13,
      n_edge23, n_edge45, n_edge46, n_edge57, n_edge67,
      n_edge15, n_edge37, n_edge04, n_edge26, n_corner);  

  dir_nodes.insert(dir_nodes.end(), n_top.begin(), n_top.end());
  dir_nodes.insert(dir_nodes.end(), n_right.begin(), n_right.end());

  dir_nodes.insert(dir_nodes.end(), n_edge45.begin(), n_edge45.end());
  dir_nodes.insert(dir_nodes.end(), n_edge67.begin(), n_edge67.end());
  dir_nodes.insert(dir_nodes.end(), n_edge57.begin(), n_edge57.end());
  dir_nodes.insert(dir_nodes.end(), n_edge46.begin(), n_edge46.end());

  dir_nodes.insert(dir_nodes.end(), n_edge23.begin(), n_edge23.end());
  dir_nodes.insert(dir_nodes.end(), n_edge26.begin(), n_edge26.end());
  dir_nodes.insert(dir_nodes.end(), n_edge37.begin(), n_edge37.end());

  dir_nodes.insert(dir_nodes.end(), n_corner[2]);
  dir_nodes.insert(dir_nodes.end(), n_corner[3]);
  dir_nodes.insert(dir_nodes.end(), n_corner[4]);
  dir_nodes.insert(dir_nodes.end(), n_corner[5]);
  dir_nodes.insert(dir_nodes.end(), n_corner[6]);
  dir_nodes.insert(dir_nodes.end(), n_corner[7]);

  per_master_nodes.insert(per_master_nodes.end(), n_front.begin(), n_front.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_back.begin(), n_back.end());

  per_master_nodes.insert(per_master_nodes.end(), n_edge15.begin(), n_edge15.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge04.begin(), n_edge04.end());

  per_master_nodes.insert(per_master_nodes.end(), n_edge13.begin(), n_edge13.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge02.begin(), n_edge02.end());

  per_master_nodes.insert(per_master_nodes.end(), n_corner[1]);
  per_slave_nodes.insert(per_slave_nodes.end(), n_corner[0]);

  num_dir_nodes = dir_nodes.size();
  num_per_nodes = per_master_nodes.size();

  std::cout<<"-----> Top: Dirichlet; Bottom: Nothing; \n";
  std::cout<<"       Front -- Back C0 periodic; \n";
  std::cout<<"       Left: Nothing; Right: Dirichlet"<<std::endl;
}


void BoundaryCond::BC_type_20( const IMesh * const &mesh )
{
  std::vector<int> n_front, n_back, n_left, n_right, n_top, n_bottom;
  std::vector<int> n_edge01, n_edge02, n_edge13, n_edge23,
    n_edge45, n_edge46, n_edge57, n_edge67, n_edge15,
    n_edge37, n_edge04, n_edge26, n_corner;

  Generate_BCNodes_B(mesh, n_front, n_back, n_left, n_right,
      n_top, n_bottom, n_edge01, n_edge02, n_edge13,
      n_edge23, n_edge45, n_edge46, n_edge57, n_edge67,
      n_edge15, n_edge37, n_edge04, n_edge26, n_corner);  

  per_master_nodes.insert(per_master_nodes.end(), n_front.begin(), n_front.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_back.begin(), n_back.end());

  per_master_nodes.insert(per_master_nodes.end(), n_edge15.begin(), n_edge15.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge04.begin(), n_edge04.end());

  per_master_nodes.insert(per_master_nodes.end(), n_edge13.begin(), n_edge13.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge02.begin(), n_edge02.end());

  per_master_nodes.insert(per_master_nodes.end(), n_corner[1]);
  per_slave_nodes.insert(per_slave_nodes.end(), n_corner[0]);

  num_dir_nodes = dir_nodes.size();
  num_per_nodes = per_master_nodes.size();

  std::cout<<"-----> Top, Bottom, Left, Right: Nothing; \n";
  std::cout<<"       Front -- Back C0 periodic."<<std::endl;
}


void BoundaryCond::print_face_edge_corner_nodes( const IMesh * const &mesh) const
{
  std::vector<int> front, back, left, right, top, bottom;
  std::vector<int> edge01, edge02, edge13, edge23, edge45, 
    edge46, edge57, edge67, edge15, edge37, edge04, edge26, corner;

  Generate_BCNodes_B(mesh, front, back, left, right, top, bottom, edge01,
      edge02, edge13, edge23, edge45, edge46, edge57, edge67, edge15, edge37,
      edge04, edge26, corner);
  std::cout<<std::endl;
  std::cout<<"=============================================== "<<std::endl;
  std::cout<<"front: ";
  VEC_T::print(front);
  std::cout<<"back: ";
  VEC_T::print(back);
  std::cout<<"left: ";
  VEC_T::print(left);
  std::cout<<"right: ";
  VEC_T::print(right);
  std::cout<<"top: ";
  VEC_T::print(top);
  std::cout<<"bottom: ";
  VEC_T::print(bottom);

  std::cout<<"Edge01: ";
  VEC_T::print(edge01);
  std::cout<<"Edge02: ";
  VEC_T::print(edge02);
  std::cout<<"Edge13: ";
  VEC_T::print(edge13);
  std::cout<<"Edge23: ";
  VEC_T::print(edge23);
  std::cout<<"Edge45: ";
  VEC_T::print(edge45);
  std::cout<<"Edge46: ";
  VEC_T::print(edge46);
  std::cout<<"Edge57: ";
  VEC_T::print(edge57);
  std::cout<<"Edge67: ";
  VEC_T::print(edge67);
  std::cout<<"Edge15: ";
  VEC_T::print(edge15);
  std::cout<<"Edge37: ";
  VEC_T::print(edge37);
  std::cout<<"Edge04: ";
  VEC_T::print(edge04);
  std::cout<<"Edge26: ";
  VEC_T::print(edge26);
  std::cout<<std::endl<<"corner: ";
  VEC_T::print(corner);
  std::cout<<"=============================================== "<<std::endl;
}

void BoundaryCond::print_face_edge_corner_nodes_2( const IMesh * const &mesh) const
{
  std::vector<int> front, back, left, right, top, bottom;
  std::vector<int> edge01, edge02, edge13, edge23, edge45, 
    edge46, edge57, edge67, edge15, edge37, edge04, edge26, corner;

  Generate_BCNodes_C(mesh, front, back, left, right, top, bottom, edge01,
      edge02, edge13, edge23, edge45, edge46, edge57, edge67, edge15, edge37,
      edge04, edge26, corner);
  std::cout<<std::endl;
  std::cout<<"=============================================== "<<std::endl;
  std::cout<<"front: ";
  VEC_T::print(front);
  std::cout<<"back: ";
  VEC_T::print(back);
  std::cout<<"left: ";
  VEC_T::print(left);
  std::cout<<"right: ";
  VEC_T::print(right);
  std::cout<<"top: ";
  VEC_T::print(top);
  std::cout<<"bottom: ";
  VEC_T::print(bottom);

  std::cout<<"Edge01: ";
  VEC_T::print(edge01);
  std::cout<<"Edge02: ";
  VEC_T::print(edge02);
  std::cout<<"Edge13: ";
  VEC_T::print(edge13);
  std::cout<<"Edge23: ";
  VEC_T::print(edge23);
  std::cout<<"Edge45: ";
  VEC_T::print(edge45);
  std::cout<<"Edge46: ";
  VEC_T::print(edge46);
  std::cout<<"Edge57: ";
  VEC_T::print(edge57);
  std::cout<<"Edge67: ";
  VEC_T::print(edge67);
  std::cout<<"Edge15: ";
  VEC_T::print(edge15);
  std::cout<<"Edge37: ";
  VEC_T::print(edge37);
  std::cout<<"Edge04: ";
  VEC_T::print(edge04);
  std::cout<<"Edge26: ";
  VEC_T::print(edge26);
  std::cout<<std::endl<<"corner: ";
  VEC_T::print(corner);
  std::cout<<"=============================================== "<<std::endl;
}


void BoundaryCond::BC_type_21( const IMesh * const &mesh )
{
  Generate_BCElems_A(mesh, front_elem, back_elem, left_elem, right_elem,
      top_elem, bottom_elem);

  front_elem.clear();
  back_elem.clear();
  left_elem.clear();
  right_elem.clear();

  num_front_elem  = front_elem.size();
  num_back_elem   = back_elem.size();
  num_left_elem   = left_elem.size();
  num_right_elem  = right_elem.size();
  num_top_elem    = top_elem.size();
  num_bottom_elem = bottom_elem.size();

  std::vector<int> nodes_front, nodes_back, nodes_left, nodes_right, nodes_top, nodes_bottom;

  Generate_BCNodes_A(mesh, nodes_front, nodes_back, nodes_left,
      nodes_right, nodes_top, nodes_bottom );
  dir_nodes.clear();
  dir_nodes.insert(dir_nodes.end(), nodes_front.begin(), nodes_front.end());
  dir_nodes.insert(dir_nodes.end(), nodes_back.begin(), nodes_back.end());
  dir_nodes.insert(dir_nodes.end(), nodes_left.begin(), nodes_left.end());
  dir_nodes.insert(dir_nodes.end(), nodes_right.begin(), nodes_right.end());

  num_dir_nodes = dir_nodes.size();

  std::cout<<"-----> NBC on top and bottom, rest face impose Dirichlet BC. \n";
}

void BoundaryCond::BC_type_22( const IMesh * const &mesh )
{
  // Specify Top and Bottom as NBC face
  Generate_BCElems_A(mesh, front_elem, back_elem, left_elem, right_elem,
      top_elem, bottom_elem);

  front_elem.clear();
  back_elem.clear();
  left_elem.clear();
  right_elem.clear();

  num_front_elem  = front_elem.size();
  num_back_elem   = back_elem.size();
  num_left_elem   = left_elem.size();
  num_right_elem  = right_elem.size();
  num_top_elem    = top_elem.size();
  num_bottom_elem = bottom_elem.size();

  // Extract nodal info 
  std::vector<int> n_front, n_back, n_left, n_right, n_top, n_bottom;
  std::vector<int> n_edge01, n_edge02, n_edge13, n_edge23,
    n_edge45, n_edge46, n_edge57, n_edge67, n_edge15,
    n_edge37, n_edge04, n_edge26, n_corner;

  Generate_BCNodes_B(mesh, n_front, n_back, n_left, n_right,
      n_top, n_bottom, n_edge01, n_edge02, n_edge13,
      n_edge23, n_edge45, n_edge46, n_edge57, n_edge67,
      n_edge15, n_edge37, n_edge04, n_edge26, n_corner);  

  // Insert Dir nodes
  dir_nodes.insert(dir_nodes.end(), n_left.begin(), n_left.end());

  dir_nodes.insert(dir_nodes.end(), n_edge04.begin(), n_edge04.end());
  dir_nodes.insert(dir_nodes.end(), n_edge15.begin(), n_edge15.end());
  dir_nodes.insert(dir_nodes.end(), n_edge01.begin(), n_edge01.end());
  dir_nodes.insert(dir_nodes.end(), n_edge45.begin(), n_edge45.end());

  dir_nodes.insert(dir_nodes.end(), n_corner[0]);
  dir_nodes.insert(dir_nodes.end(), n_corner[1]);
  dir_nodes.insert(dir_nodes.end(), n_corner[4]);
  dir_nodes.insert(dir_nodes.end(), n_corner[5]);

  // Insert per nodes
  per_master_nodes.insert(per_master_nodes.end(), n_front.begin(), n_front.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_back.begin(), n_back.end());

  per_master_nodes.insert(per_master_nodes.end(), n_edge37.begin(), n_edge37.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge26.begin(), n_edge26.end());

  per_master_nodes.insert(per_master_nodes.end(), n_edge13.begin(), n_edge13.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge02.begin(), n_edge02.end());

  per_master_nodes.insert(per_master_nodes.end(), n_edge57.begin(), n_edge57.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge46.begin(), n_edge46.end());

  per_master_nodes.insert(per_master_nodes.end(), n_corner[3]);
  per_slave_nodes.insert(per_slave_nodes.end(), n_corner[2]);

  per_master_nodes.insert(per_master_nodes.end(), n_corner[7]);
  per_slave_nodes.insert(per_slave_nodes.end(), n_corner[6]);

  num_dir_nodes = dir_nodes.size();
  num_per_nodes = per_master_nodes.size();

  std::cout<<"-----> NBC on top and bottom. Front-Back C0 periodic.";
  std::cout<<" Right nothing. Left: Dirichlet. \n";
}

void BoundaryCond::BC_type_23( const IMesh * const &mesh )
{
  // ----- From here is a direct copy from BC_type_22 -----
  // Specify Top and Bottom as NBC face
  Generate_BCElems_A(mesh, front_elem, back_elem, left_elem, right_elem,
      top_elem, bottom_elem);

  front_elem.clear();
  back_elem.clear();
  left_elem.clear();
  right_elem.clear();

  num_front_elem  = front_elem.size();
  num_back_elem   = back_elem.size();
  num_left_elem   = left_elem.size();
  num_right_elem  = right_elem.size();
  num_top_elem    = top_elem.size();
  num_bottom_elem = bottom_elem.size();

  // Extract nodal info 
  std::vector<int> n_front, n_back, n_left, n_right, n_top, n_bottom;
  std::vector<int> n_edge01, n_edge02, n_edge13, n_edge23,
    n_edge45, n_edge46, n_edge57, n_edge67, n_edge15,
    n_edge37, n_edge04, n_edge26, n_corner;

  Generate_BCNodes_B(mesh, n_front, n_back, n_left, n_right,
      n_top, n_bottom, n_edge01, n_edge02, n_edge13,
      n_edge23, n_edge45, n_edge46, n_edge57, n_edge67,
      n_edge15, n_edge37, n_edge04, n_edge26, n_corner);  

  // Insert Dir nodes
  dir_nodes.insert(dir_nodes.end(), n_left.begin(), n_left.end());

  dir_nodes.insert(dir_nodes.end(), n_edge04.begin(), n_edge04.end());
  dir_nodes.insert(dir_nodes.end(), n_edge15.begin(), n_edge15.end());
  dir_nodes.insert(dir_nodes.end(), n_edge01.begin(), n_edge01.end());
  dir_nodes.insert(dir_nodes.end(), n_edge45.begin(), n_edge45.end());

  dir_nodes.insert(dir_nodes.end(), n_corner[0]);
  dir_nodes.insert(dir_nodes.end(), n_corner[1]);
  dir_nodes.insert(dir_nodes.end(), n_corner[4]);
  dir_nodes.insert(dir_nodes.end(), n_corner[5]);

  // Insert per nodes
  per_master_nodes.insert(per_master_nodes.end(), n_front.begin(), n_front.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_back.begin(), n_back.end());

  per_master_nodes.insert(per_master_nodes.end(), n_edge37.begin(), n_edge37.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge26.begin(), n_edge26.end());

  per_master_nodes.insert(per_master_nodes.end(), n_edge13.begin(), n_edge13.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge02.begin(), n_edge02.end());

  per_master_nodes.insert(per_master_nodes.end(), n_edge57.begin(), n_edge57.end());
  per_slave_nodes.insert(per_slave_nodes.end(), n_edge46.begin(), n_edge46.end());

  per_master_nodes.insert(per_master_nodes.end(), n_corner[3]);
  per_slave_nodes.insert(per_slave_nodes.end(), n_corner[2]);

  per_master_nodes.insert(per_master_nodes.end(), n_corner[7]);
  per_slave_nodes.insert(per_slave_nodes.end(), n_corner[6]);
  // ----- Finish the BC_type_22 part -----

  // ----- Manually impose master-slave on right face -----
  const int nFunc_x = mesh->get_nFunc_x();
  const int nFunc_y = mesh->get_nFunc_y();
  const int nFunc_z = mesh->get_nFunc_z();

  for( int zz = 0; zz < nFunc_z; ++zz )
  {
    int mnode = zz * nFunc_x * nFunc_y + (nFunc_y - 1) * nFunc_x + nFunc_x - 1;
    for( int xx = 1; xx < nFunc_x - 1; ++xx )
    {
      int snode = zz * nFunc_x * nFunc_y + (nFunc_y-1) * nFunc_x + xx;

      per_master_nodes.push_back(mnode);
      per_slave_nodes.push_back(snode);
    }
  }
  // ----- Finish the right face imposition -----

  num_dir_nodes = dir_nodes.size();
  num_per_nodes = per_master_nodes.size();

  std::cout<<"-----> NBC on top and bottom; Front-Back C0 periodic;";
  std::cout<<" Right: Following edge 37; Left: Dirichlet. \n";
}

void BoundaryCond::BC_type_101( const IMesh * const &mesh )
{
  const int numPat = mesh->get_num_patch();

  if(numPat != 3)
  {
    std::cerr<<"Error: BC_type_101 is designed for 3-patch coronary geometry only. \n";
    exit(EXIT_FAILURE);
  }

  // We extract the boundary nodes for each patch
  // Patch 0
  int nFunc_x_0, nFunc_y_0, nFunc_z_0, f_start_0;
  nFunc_x_0 = mesh->get_patch_mesh(0)->get_nFunc_x();
  nFunc_y_0 = mesh->get_patch_mesh(0)->get_nFunc_y();
  nFunc_z_0 = mesh->get_patch_mesh(0)->get_nFunc_z();
  f_start_0 = mesh->get_patch_mesh(0)->get_nFunc_start();
  
  int * back_front_0, * left_right_0, * bottom_top_0, * edge01_0, 
      * edge02_0, * edge04_0, * corner_0; 

  int nbf_0, nlr_0, nbt_0, nex_0, ney_0, nez_0;

  Generate_BCNodes_D(nFunc_x_0, nFunc_y_0, nFunc_z_0, f_start_0,
      back_front_0, nbf_0, left_right_0, nlr_0, bottom_top_0, nbt_0,
      edge01_0, nex_0, edge02_0, ney_0, edge04_0, nez_0, corner_0);
  
  // Patch 1
  int nFunc_x_1, nFunc_y_1, nFunc_z_1, f_start_1;
  nFunc_x_1 = mesh->get_patch_mesh(1)->get_nFunc_x();
  nFunc_y_1 = mesh->get_patch_mesh(1)->get_nFunc_y();
  nFunc_z_1 = mesh->get_patch_mesh(1)->get_nFunc_z();
  f_start_1 = mesh->get_patch_mesh(1)->get_nFunc_start();
 
  int * back_front_1, * left_right_1, * bottom_top_1, * edge01_1, 
      * edge02_1, * edge04_1, * corner_1; 

  int nbf_1, nlr_1, nbt_1, nex_1, ney_1, nez_1;

  Generate_BCNodes_D(nFunc_x_1, nFunc_y_1, nFunc_z_1, f_start_1,
      back_front_1, nbf_1, left_right_1, nlr_1, bottom_top_1, nbt_1,
      edge01_1, nex_1, edge02_1, ney_1, edge04_1, nez_1, corner_1);


  // Patch 2
  int nFunc_x_2, nFunc_y_2, nFunc_z_2, f_start_2;
  nFunc_x_2 = mesh->get_patch_mesh(2)->get_nFunc_x();
  nFunc_y_2 = mesh->get_patch_mesh(2)->get_nFunc_y();
  nFunc_z_2 = mesh->get_patch_mesh(2)->get_nFunc_z();
  f_start_2 = mesh->get_patch_mesh(2)->get_nFunc_start();
 
  int * back_front_2, * left_right_2, * bottom_top_2, * edge01_2, 
      * edge02_2, * edge04_2, * corner_2; 

  int nbf_2, nlr_2, nbt_2, nex_2, ney_2, nez_2;

  Generate_BCNodes_D(nFunc_x_2, nFunc_y_2, nFunc_z_2, f_start_2,
      back_front_2, nbf_2, left_right_2, nlr_2, bottom_top_2, nbt_2,
      edge01_2, nex_2, edge02_2, ney_2, edge04_2, nez_2, corner_2);


  // This is a 3-patch coronary geometry, with top/bottom faces of each patch as
  // the interface
  if( (nFunc_x_0 != nFunc_x_1) || (nFunc_y_0 != nFunc_y_1) )
  {
    std::cerr<<"Error: Patch 0 and Patch 1 interface function does not match. \n";
    exit(EXIT_FAILURE);
  }

  if( (nFunc_x_0 != nFunc_x_2) || (nFunc_y_0 != nFunc_y_2) )
  {
    std::cerr<<"Error: Patch 0 and Patch 2 interface function does not match. \n";
    exit(EXIT_FAILURE);
  }
  const int nFunc_x = nFunc_x_0;
  const int nFunc_y = nFunc_y_0;
  
  const int nFxFy = nFunc_x * nFunc_y;

  if( (nFunc_x_0-1) % 4 != 0 )
  {
    std::cerr<<"Error: nFunc_x = 4 (2+inserted_s) + 1. \n ";
    exit(EXIT_FAILURE);
  }

  // num_ins = number of basis functions between two C0 knots
  // The structure is
  // 0 , ... ,  2+num_ins, ..., 2(2+num_ins), ..., 3(2+num_ins), ... , 4(2+num_ins) 
  // e.g. num_ins = 0 (aka no knot insertion or degree elevation)
  // the basis in s direction has the following structure
  // 0,1,2,3,4,5,6,7,8 ===> 3   2   1
  //                        4       0/8
  //                        5   6   7
  const int num_ins = (nFunc_x_0-1) / 4 - 2;

  // The matching pattern is patch 0's top, patch 1's top and patch 2's bottom
  // meets at the interface
  // The first group of master-slave relationship : 
  // patch 0 : 2,6 on top is master,
  // patch 1 : 2,6 on top is slave
  // patch 2 : 2,6 on bottom is slave
  const int patch0_top_start = (nFunc_z_0-1) * nFxFy;
  const int patch1_top_start = (nFunc_z_1-1) * nFxFy + f_start_1;
  
  per_master_nodes.reserve(nFxFy*2 + nFunc_y * (nFunc_z_0 + nFunc_z_1 + nFunc_z_2) );
  per_slave_nodes.reserve(nFxFy*2 + nFunc_y * (nFunc_z_0 + nFunc_z_1 + nFunc_z_2) );
 
  for(int ii=1; ii<nFunc_y-1; ++ii)
  {
    // 2
    per_master_nodes.push_back( patch0_top_start + ii * nFunc_x + 2 + num_ins );
    per_slave_nodes.push_back(  patch1_top_start + ii * nFunc_x + 2 + num_ins );
    per_master_nodes.push_back( patch0_top_start + ii * nFunc_x + 2 + num_ins );
    per_slave_nodes.push_back( f_start_2 + ii * nFunc_x + 2 + num_ins );
    // 6
    per_master_nodes.push_back( patch0_top_start + ii * nFunc_x + 6 + 3*num_ins );
    per_slave_nodes.push_back(  patch1_top_start + ii * nFunc_x + 6 + 3*num_ins );
    per_master_nodes.push_back( patch0_top_start + ii * nFunc_x + 6 + 3*num_ins );
    per_slave_nodes.push_back(  f_start_2        + ii * nFunc_x + 6 + 3*num_ins );
  } 
  
  for(int ii=1; ii<nFunc_y-1; ++ii)
  {
    const int level_t = ii * nFunc_x;
    // The second group of master-slave relationship :
    // patch 0 : 4 - 6 top
    // patch 1 : 4 - 6 top
    for( int jj= 3 + num_ins; jj<6+3*num_ins; ++jj )
    {
      per_master_nodes.push_back( patch0_top_start + level_t + jj );
      per_slave_nodes.push_back( patch1_top_start + level_t + jj );
    }
    // The third group of master-slave relationship
    // patch 0 : 1 top  --> patch 2 : 3 reverse order
    // patch 0 : 7 top  --> patch 2 : 5 reverse order
    for(int jj=1; jj<2+num_ins; ++jj)
    {
      per_master_nodes.push_back( patch0_top_start + level_t + jj);
      per_slave_nodes.push_back(f_start_2 + level_t + 4 + 2*num_ins - jj);    
      
      per_master_nodes.push_back( patch0_top_start + level_t + 6 + 3 * num_ins + jj);
      per_slave_nodes.push_back( f_start_2 + level_t + 6 + 3 * num_ins - jj);
    
      // The fourth group of master-slave relationship
      // patch 1 : 1 top --> patch 2 : 1 bottom
      // patch 1 : 7 top --> patch 2 : 7 bottom
      per_master_nodes.push_back( patch1_top_start + level_t + jj);
      per_slave_nodes.push_back( f_start_2 + level_t + jj );
      per_master_nodes.push_back( patch1_top_start + level_t + 6 + 3 * num_ins + jj);
      per_slave_nodes.push_back( f_start_2 + level_t + 6 + 3 * num_ins + jj );
    }
    // The fifth group of master-slave relationship
    // patch 0 : 0 top --> patch 0 : 8 top
    // patch 0 : 0 top --> patch 2 : 5 bottom
    per_master_nodes.push_back(patch0_top_start + level_t);
    per_slave_nodes.push_back( patch0_top_start + level_t + 8 + 4*num_ins);
    per_master_nodes.push_back(patch0_top_start + level_t);
    per_slave_nodes.push_back( f_start_2 + level_t + 4 + 2 * num_ins);
    
    // The sixth group of master-slave relationship
    // patch 1 : 0 top --> patch 1 : 8 top
    // patch 1 : 0 top --> patch 2 : 0 bottom
    // patch 1 : 0 top --> patch 2 : 8 bottom
    per_master_nodes.push_back(patch1_top_start + level_t);
    per_slave_nodes.push_back(patch1_top_start + level_t + 8 + 4*num_ins);
    per_master_nodes.push_back(patch1_top_start + level_t);
    per_slave_nodes.push_back(f_start_2 + level_t );
    per_master_nodes.push_back(patch1_top_start + level_t);
    per_slave_nodes.push_back(f_start_2 + level_t + 8 + 4 * num_ins);
  }
  

  // Periodic paris within each patch
  // patch 0 : front - back, e37 - e26, e13 - e02, c3-c2
  per_master_nodes.insert(per_master_nodes.end(), back_front_0 + nbf_0, back_front_0 + 2*nbf_0);
  per_slave_nodes.insert(per_slave_nodes.end(), back_front_0, back_front_0 + nbf_0);
  
  per_master_nodes.insert(per_master_nodes.end(), edge04_0+3*nez_0, edge04_0+4*nez_0);
  per_slave_nodes.insert(per_slave_nodes.end(), edge04_0+2*nez_0, edge04_0+3*nez_0);

  per_master_nodes.insert(per_master_nodes.end(), edge02_0+ney_0, edge02_0+2*ney_0);
  per_slave_nodes.insert(per_slave_nodes.end(), edge02_0, edge02_0+ney_0);
  
  per_master_nodes.push_back(corner_0[3]);
  per_slave_nodes.push_back(corner_0[2]);
  
  per_master_nodes.push_back(corner_0[7]);
  per_slave_nodes.push_back(corner_0[6]);

  // patch 1 : front - back, e37 - e26, e13 - e02, c3-c2
  per_master_nodes.insert(per_master_nodes.end(), back_front_1 + nbf_1, back_front_1 + 2*nbf_1);
  per_slave_nodes.insert(per_slave_nodes.end(), back_front_1, back_front_1 + nbf_1);
  
  per_master_nodes.insert(per_master_nodes.end(), edge04_1+3*nez_1, edge04_1+4*nez_1);
  per_slave_nodes.insert(per_slave_nodes.end(), edge04_1+2*nez_1, edge04_1+3*nez_1);

  per_master_nodes.insert(per_master_nodes.end(), edge02_1+ney_1, edge02_1+2*ney_1);
  per_slave_nodes.insert(per_slave_nodes.end(), edge02_1, edge02_1+ney_1);
  
  per_master_nodes.push_back(corner_1[3]);
  per_slave_nodes.push_back(corner_1[2]); 
  
  per_master_nodes.push_back(corner_1[7]);
  per_slave_nodes.push_back(corner_1[6]); 
  
  // patch 2 : front - back, e37 - e26, e57 - e46, c7-c6
  per_master_nodes.insert(per_master_nodes.end(), back_front_2 + nbf_2, back_front_2 + 2*nbf_2);
  per_slave_nodes.insert(per_slave_nodes.end(), back_front_2, back_front_2 + nbf_2);
  
  per_master_nodes.insert(per_master_nodes.end(), edge04_2+3*nez_2, edge04_2+4*nez_2);
  per_slave_nodes.insert(per_slave_nodes.end(), edge04_2+2*nez_2, edge04_2+3*nez_2);

  per_master_nodes.insert(per_master_nodes.end(), edge02_2+3*ney_2, edge02_2+4*ney_2);
  per_slave_nodes.insert( per_slave_nodes.end(),  edge02_2+2*ney_2, edge02_2+3*ney_2);
  
  per_master_nodes.push_back(corner_2[3]);
  per_slave_nodes.push_back(corner_2[2]); 
  
  per_master_nodes.push_back(corner_2[7]);
  per_slave_nodes.push_back(corner_2[6]); 

  VEC_T::shrink2fit(per_master_nodes);
  VEC_T::shrink2fit(per_slave_nodes);


  // Dirichlet nodes
  dir_nodes.reserve(nlr_0 + nlr_1 + nlr_2 + 2*nex_0 + 2*nex_1 + 2*nex_2 
      + 2*nez_0 + 2*nez_1 + 2*nez_2);

  // patch 0 : f_left, e01, e04, e15, c0, c1
  dir_nodes.insert(dir_nodes.end(), left_right_0, left_right_0 + nlr_0);
  dir_nodes.insert(dir_nodes.end(), edge01_0, edge01_0 + nex_0 );
  dir_nodes.insert(dir_nodes.end(), edge01_0 + 2*nex_0, edge01_0 + 3*nex_0 );
  dir_nodes.insert(dir_nodes.end(), edge04_0, edge04_0 + 2*nez_0 );
  dir_nodes.push_back(corner_0[0]);
  dir_nodes.push_back(corner_0[1]);
  dir_nodes.push_back(corner_0[4]);
  dir_nodes.push_back(corner_0[5]);

  // patch 1 : f_left, e01, e04, e15, c0, c1
  dir_nodes.insert(dir_nodes.end(), left_right_1, left_right_1 + nlr_1);
  dir_nodes.insert(dir_nodes.end(), edge01_1, edge01_1 + nex_1 );
  dir_nodes.insert(dir_nodes.end(), edge01_1 + 2*nex_1, edge01_1 + 3*nex_1 );
  dir_nodes.insert(dir_nodes.end(), edge04_1, edge04_1 + 2*nez_1 );
  dir_nodes.push_back(corner_1[0]);
  dir_nodes.push_back(corner_1[1]);
  dir_nodes.push_back(corner_1[4]);
  dir_nodes.push_back(corner_1[5]);

  // patch 2 : f_left, e04, e15, e45, c4, c5
  dir_nodes.insert(dir_nodes.end(), left_right_2, left_right_2 + nlr_2);
  dir_nodes.insert(dir_nodes.end(), edge01_2, edge01_2 + nex_2 );
  dir_nodes.insert(dir_nodes.end(), edge01_2 + 2*nex_2, edge01_2 + 3*nex_2 );
  dir_nodes.insert(dir_nodes.end(), edge04_2, edge04_2 + 2*nez_2 );
  dir_nodes.push_back(corner_2[0]);
  dir_nodes.push_back(corner_2[1]);
  dir_nodes.push_back(corner_2[4]);
  dir_nodes.push_back(corner_2[5]);
  
  VEC_T::shrink2fit(dir_nodes);

  delete [] back_front_0; delete [] left_right_0; delete [] bottom_top_0;
  delete [] edge01_0; delete [] edge02_0; delete [] edge04_0; delete [] corner_0;
  delete [] back_front_1; delete [] left_right_1; delete [] bottom_top_1;
  delete [] edge01_1; delete [] edge02_1; delete [] edge04_1; delete [] corner_1;
  delete [] back_front_2; delete [] left_right_2; delete [] bottom_top_2;
  delete [] edge01_2; delete [] edge02_2; delete [] edge04_2; delete [] corner_2;

  num_per_nodes = per_master_nodes.size();
  num_dir_nodes = dir_nodes.size();

  // Element face integral specification
  front_elem.clear(); back_elem.clear(); left_elem.clear(); right_elem.clear();
  top_elem.clear(); bottom_elem.clear();
  
  num_front_elem = 0; num_back_elem = 0;
  num_left_elem = 0; num_right_elem = 0;

  VEC_T::shrink2fit(front_elem); VEC_T::shrink2fit(back_elem);
  VEC_T::shrink2fit(left_elem); VEC_T::shrink2fit(right_elem);

  std::vector<int> front, back, left, right, top, bottom;
  Generate_BCElems_B(mesh->get_patch_mesh(0), front, back, left, right, top, bottom);

  bottom_elem.insert(bottom_elem.end(), bottom.begin(), bottom.end());

  Generate_BCElems_B(mesh->get_patch_mesh(1), front, back, left, right, top, bottom);

  bottom_elem.insert(bottom_elem.end(), bottom.begin(), bottom.end());
  
  Generate_BCElems_B(mesh->get_patch_mesh(2), front, back, left, right, top, bottom);

  top_elem.insert(top_elem.end(), top.begin(), top.end());

  VEC_T::shrink2fit(bottom_elem); VEC_T::shrink2fit(top_elem);

  num_top_elem = top_elem.size();
  num_bottom_elem = bottom_elem.size();

  std::cout<<"-----> 3-Patch Coronary geometry: NBC on patch 0 & 1 bottom, patch 2 top; ";
  std::cout<<" Front-Back C0 periodic; Left Dirichlet, Right nothing. \n";
}


void BoundaryCond::BC_type_102( const IMesh * const &mesh )
{
  const int numPat = mesh->get_num_patch();

  if(numPat != 3)
  {
    std::cerr<<"Error: BC_type_102 is designed for 3-patch coronary geometry only. \n";
    exit(EXIT_FAILURE);
  }

  // We extract the boundary nodes for each patch
  // Patch 0
  int nFunc_x_0, nFunc_y_0, nFunc_z_0, f_start_0;
  nFunc_x_0 = mesh->get_patch_mesh(0)->get_nFunc_x();
  nFunc_y_0 = mesh->get_patch_mesh(0)->get_nFunc_y();
  nFunc_z_0 = mesh->get_patch_mesh(0)->get_nFunc_z();
  f_start_0 = mesh->get_patch_mesh(0)->get_nFunc_start();
  
  int * back_front_0, * left_right_0, * bottom_top_0, * edge01_0, 
      * edge02_0, * edge04_0, * corner_0; 

  int nbf_0, nlr_0, nbt_0, nex_0, ney_0, nez_0;

  Generate_BCNodes_D(nFunc_x_0, nFunc_y_0, nFunc_z_0, f_start_0,
      back_front_0, nbf_0, left_right_0, nlr_0, bottom_top_0, nbt_0,
      edge01_0, nex_0, edge02_0, ney_0, edge04_0, nez_0, corner_0);
  
  // Patch 1
  int nFunc_x_1, nFunc_y_1, nFunc_z_1, f_start_1;
  nFunc_x_1 = mesh->get_patch_mesh(1)->get_nFunc_x();
  nFunc_y_1 = mesh->get_patch_mesh(1)->get_nFunc_y();
  nFunc_z_1 = mesh->get_patch_mesh(1)->get_nFunc_z();
  f_start_1 = mesh->get_patch_mesh(1)->get_nFunc_start();
 
  int * back_front_1, * left_right_1, * bottom_top_1, * edge01_1, 
      * edge02_1, * edge04_1, * corner_1; 

  int nbf_1, nlr_1, nbt_1, nex_1, ney_1, nez_1;

  Generate_BCNodes_D(nFunc_x_1, nFunc_y_1, nFunc_z_1, f_start_1,
      back_front_1, nbf_1, left_right_1, nlr_1, bottom_top_1, nbt_1,
      edge01_1, nex_1, edge02_1, ney_1, edge04_1, nez_1, corner_1);


  // Patch 2
  int nFunc_x_2, nFunc_y_2, nFunc_z_2, f_start_2;
  nFunc_x_2 = mesh->get_patch_mesh(2)->get_nFunc_x();
  nFunc_y_2 = mesh->get_patch_mesh(2)->get_nFunc_y();
  nFunc_z_2 = mesh->get_patch_mesh(2)->get_nFunc_z();
  f_start_2 = mesh->get_patch_mesh(2)->get_nFunc_start();
 
  int * back_front_2, * left_right_2, * bottom_top_2, * edge01_2, 
      * edge02_2, * edge04_2, * corner_2; 

  int nbf_2, nlr_2, nbt_2, nex_2, ney_2, nez_2;

  Generate_BCNodes_D(nFunc_x_2, nFunc_y_2, nFunc_z_2, f_start_2,
      back_front_2, nbf_2, left_right_2, nlr_2, bottom_top_2, nbt_2,
      edge01_2, nex_2, edge02_2, ney_2, edge04_2, nez_2, corner_2);


  // This is a 3-patch coronary geometry, with top/bottom faces of each patch as
  // the interface
  if( (nFunc_x_0 != nFunc_x_1) || (nFunc_y_0 != nFunc_y_1) )
  {
    std::cerr<<"Error: Patch 0 and Patch 1 interface function does not match. \n";
    exit(EXIT_FAILURE);
  }

  if( (nFunc_x_0 != nFunc_x_2) || (nFunc_y_0 != nFunc_y_2) )
  {
    std::cerr<<"Error: Patch 0 and Patch 2 interface function does not match. \n";
    exit(EXIT_FAILURE);
  }
  const int nFunc_x = nFunc_x_0;
  const int nFunc_y = nFunc_y_0;
  
  const int nFxFy = nFunc_x * nFunc_y;

  if( (nFunc_x-1) % 4 != 0 )
  {
    std::cerr<<"Error: nFunc_x = 4 (2+inserted_s) + 1. \n ";
    exit(EXIT_FAILURE);
  }

  // num_ins = number of basis functions between two C0 knots
  // The structure is
  // 0 , ... ,  2+num_ins, ..., 2(2+num_ins), ..., 3(2+num_ins), ... , 4(2+num_ins) 
  // e.g. num_ins = 0 (aka no knot insertion or degree elevation)
  // the basis in s direction has the following structure
  // 0,1,2,3,4,5,6,7,8 ===> 3   2   1
  //                        4       0/8
  //                        5   6   7
  const int num_ins = (nFunc_x-1) / 4 - 2;

  // ============= Interface Condition =============
  // The matching pattern is patch 0's top, patch 1's bottom and patch 2's top
  // meets at the interface
  // The first group of master-slave relationship : 
  // patch 0 : 2,6 on top is master,
  // patch 1 : 2,6 on bottom is slave
  // patch 2 : 2,6 on top is slave
  const int patch0_start = (nFunc_z_0-1) * nFxFy;
  const int patch1_start = f_start_1;
  const int patch2_start = (nFunc_z_2-1) * nFxFy + f_start_2;
  
  per_master_nodes.reserve(nFxFy*2 + nFunc_y * (nFunc_z_0 + nFunc_z_1 + nFunc_z_2) );
  per_slave_nodes.reserve(nFxFy*2 + nFunc_y * (nFunc_z_0 + nFunc_z_1 + nFunc_z_2) );
 
  for(int ii=1; ii<nFunc_y-1; ++ii)
  {
    // 2
    per_master_nodes.push_back( patch0_start + ii * nFunc_x + 2 + num_ins );
    per_slave_nodes.push_back(  patch1_start + ii * nFunc_x + 2 + num_ins );
    per_master_nodes.push_back( patch0_start + ii * nFunc_x + 2 + num_ins );
    per_slave_nodes.push_back(  patch2_start + ii * nFunc_x + 2 + num_ins );
    // 6
    per_master_nodes.push_back( patch0_start + ii * nFunc_x + 6 + 3*num_ins );
    per_slave_nodes.push_back(  patch1_start + ii * nFunc_x + 6 + 3*num_ins );
    per_master_nodes.push_back( patch0_start + ii * nFunc_x + 6 + 3*num_ins );
    per_slave_nodes.push_back(  patch2_start + ii * nFunc_x + 6 + 3*num_ins );
  } 
  
  for(int ii=1; ii<nFunc_y-1; ++ii)
  {
    const int level_t = ii * nFunc_x;
    // The second group of master-slave relationship :
    // patch 0 : 4 - 6 top
    // patch 1 : 4 - 6 top
    for( int jj= 3 + num_ins; jj<6+3*num_ins; ++jj )
    {
      per_master_nodes.push_back( patch0_start + level_t + jj );
      per_slave_nodes.push_back( patch1_start + level_t + jj );
    }
    // The third group of master-slave relationship
    // patch 0 : 1 top  --> patch 2 : 3 reverse order
    // patch 0 : 7 top  --> patch 2 : 5 reverse order
    for(int jj=1; jj<2+num_ins; ++jj)
    {
      per_master_nodes.push_back( patch0_start + level_t + jj);
      per_slave_nodes.push_back(  patch2_start + level_t + 4 + 2*num_ins - jj);    
      
      per_master_nodes.push_back( patch0_start + level_t + 6 + 3 * num_ins + jj);
      per_slave_nodes.push_back(  patch2_start + level_t + 6 + 3 * num_ins - jj);
    
      // The fourth group of master-slave relationship
      // patch 1 : 1 top --> patch 2 : 1 bottom
      // patch 1 : 7 top --> patch 2 : 7 bottom
      per_master_nodes.push_back( patch1_start + level_t + jj);
      per_slave_nodes.push_back(  patch2_start + level_t + jj );
      per_master_nodes.push_back( patch1_start + level_t + 6 + 3 * num_ins + jj);
      per_slave_nodes.push_back(  patch2_start + level_t + 6 + 3 * num_ins + jj );
    }
    // The fifth group of master-slave relationship
    // patch 0 : 0 top --> patch 0 : 8 top
    // patch 0 : 0 top --> patch 2 : 5 bottom
    per_master_nodes.push_back(patch0_start + level_t);
    per_slave_nodes.push_back( patch0_start + level_t + 8 + 4*num_ins);
    per_master_nodes.push_back(patch0_start + level_t);
    per_slave_nodes.push_back( patch2_start + level_t + 4 + 2 * num_ins);
    
    // The sixth group of master-slave relationship
    // patch 1 : 0 top --> patch 1 : 8 top
    // patch 1 : 0 top --> patch 2 : 0 bottom
    // patch 1 : 0 top --> patch 2 : 8 bottom
    per_master_nodes.push_back(patch1_start + level_t);
    per_slave_nodes.push_back(patch1_start + level_t + 8 + 4*num_ins);
    per_master_nodes.push_back(patch1_start + level_t);
    per_slave_nodes.push_back(patch2_start + level_t );
    per_master_nodes.push_back(patch1_start + level_t);
    per_slave_nodes.push_back(patch2_start + level_t + 8 + 4 * num_ins);
  }
  

  // Periodic paris within each patch
  // patch 0 : front - back, e15 - e04, e13 - e02, c1-c0, c5-c4
  per_master_nodes.insert(per_master_nodes.end(), back_front_0 + nbf_0, back_front_0 + 2*nbf_0);
  per_slave_nodes.insert(per_slave_nodes.end(), back_front_0, back_front_0 + nbf_0);
  
  per_master_nodes.insert(per_master_nodes.end(), edge04_0+1*nez_0, edge04_0+2*nez_0);
  per_slave_nodes.insert(per_slave_nodes.end(), edge04_0+0*nez_0, edge04_0+1*nez_0);

  per_master_nodes.insert(per_master_nodes.end(), edge02_0+ney_0, edge02_0+2*ney_0);
  per_slave_nodes.insert(per_slave_nodes.end(), edge02_0, edge02_0+ney_0);
  
  per_master_nodes.push_back(corner_0[1]);
  per_slave_nodes.push_back(corner_0[0]);
  
  per_master_nodes.push_back(corner_0[5]);
  per_slave_nodes.push_back(corner_0[4]);

  // patch 1 : front - back, e15 - e04, e57 - e46, c5-c4, c1-c0
  per_master_nodes.insert(per_master_nodes.end(), back_front_1 + nbf_1, back_front_1 + 2*nbf_1);
  per_slave_nodes.insert(per_slave_nodes.end(), back_front_1, back_front_1 + nbf_1);
  
  per_master_nodes.insert(per_master_nodes.end(), edge04_1+1*nez_1, edge04_1+2*nez_1);
  per_slave_nodes.insert(per_slave_nodes.end(), edge04_1+0*nez_1, edge04_1+1*nez_1);

  per_master_nodes.insert(per_master_nodes.end(), edge02_1+3*ney_1, edge02_1+4*ney_1);
  per_slave_nodes.insert(per_slave_nodes.end(), edge02_1+2*ney_1, edge02_1+3*ney_1);
  
  per_master_nodes.push_back(corner_1[5]);
  per_slave_nodes.push_back(corner_1[4]); 
  
  per_master_nodes.push_back(corner_1[1]);
  per_slave_nodes.push_back(corner_1[0]); 
  
  // patch 2 : front - back, e15 - e04, e13 - e02, c5-c4, c1-c0
  per_master_nodes.insert(per_master_nodes.end(), back_front_2 + nbf_2, back_front_2 + 2*nbf_2);
  per_slave_nodes.insert(per_slave_nodes.end(), back_front_2, back_front_2 + nbf_2);
  
  per_master_nodes.insert(per_master_nodes.end(), edge04_2+1*nez_2, edge04_2+2*nez_2);
  per_slave_nodes.insert(per_slave_nodes.end(), edge04_2+0*nez_2, edge04_2+1*nez_2);

  per_master_nodes.insert(per_master_nodes.end(), edge02_2+1*ney_2, edge02_2+2*ney_2);
  per_slave_nodes.insert( per_slave_nodes.end(),  edge02_2+0*ney_2, edge02_2+1*ney_2);
  
  per_master_nodes.push_back(corner_2[5]);
  per_slave_nodes.push_back(corner_2[4]); 
  
  per_master_nodes.push_back(corner_2[1]);
  per_slave_nodes.push_back(corner_2[0]); 

  VEC_T::shrink2fit(per_master_nodes);
  VEC_T::shrink2fit(per_slave_nodes);


  // Dirichlet nodes
  dir_nodes.reserve(nlr_0 + nlr_1 + nlr_2 + 2*nex_0 + 2*nex_1 + 2*nex_2 
      + 2*nez_0 + 2*nez_1 + 2*nez_2);

  // patch 0 : f_right, e23, e67, e26, e37, c2,3,6,7
  dir_nodes.insert(dir_nodes.end(), left_right_0 + nlr_0, left_right_0 + 2*nlr_0);
  dir_nodes.insert(dir_nodes.end(), edge01_0 + nex_0, edge01_0 + 2*nex_0 );
  dir_nodes.insert(dir_nodes.end(), edge01_0 + 3*nex_0, edge01_0 + 4*nex_0 );
  dir_nodes.insert(dir_nodes.end(), edge04_0 + 2*nez_0, edge04_0 + 4*nez_0 );
  dir_nodes.push_back(corner_0[2]);
  dir_nodes.push_back(corner_0[3]);
  dir_nodes.push_back(corner_0[6]);
  dir_nodes.push_back(corner_0[7]);

  // patch 1 : f_right, e23, e67, e26, e37, c2,3,6,7
  dir_nodes.insert(dir_nodes.end(), left_right_1 + nlr_1, left_right_1 + 2*nlr_1);
  dir_nodes.insert(dir_nodes.end(), edge01_1 + nex_1, edge01_1 + 2*nex_1 );
  dir_nodes.insert(dir_nodes.end(), edge01_1 + 3*nex_1, edge01_1 + 4*nex_1 );
  dir_nodes.insert(dir_nodes.end(), edge04_1 + 2*nez_1, edge04_1 + 4*nez_1 );
  dir_nodes.push_back(corner_1[2]);
  dir_nodes.push_back(corner_1[3]);
  dir_nodes.push_back(corner_1[6]);
  dir_nodes.push_back(corner_1[7]);

  // patch 2 : f_right, e23, e67, e26, e37, c2,3,6,7
  dir_nodes.insert(dir_nodes.end(), left_right_2+nlr_2, left_right_2 + 2*nlr_2);
  dir_nodes.insert(dir_nodes.end(), edge01_2 + nex_2, edge01_2 + 2*nex_2 );
  dir_nodes.insert(dir_nodes.end(), edge01_2 + 3*nex_2, edge01_2 + 4*nex_2 );
  dir_nodes.insert(dir_nodes.end(), edge04_2 + 2*nez_2, edge04_2 + 4*nez_2 );
  dir_nodes.push_back(corner_2[2]);
  dir_nodes.push_back(corner_2[3]);
  dir_nodes.push_back(corner_2[6]);
  dir_nodes.push_back(corner_2[7]);
  
  VEC_T::shrink2fit(dir_nodes);

  delete [] back_front_0; delete [] left_right_0; delete [] bottom_top_0;
  delete [] edge01_0; delete [] edge02_0; delete [] edge04_0; delete [] corner_0;
  delete [] back_front_1; delete [] left_right_1; delete [] bottom_top_1;
  delete [] edge01_1; delete [] edge02_1; delete [] edge04_1; delete [] corner_1;
  delete [] back_front_2; delete [] left_right_2; delete [] bottom_top_2;
  delete [] edge01_2; delete [] edge02_2; delete [] edge04_2; delete [] corner_2;

  num_per_nodes = per_master_nodes.size();
  num_dir_nodes = dir_nodes.size();

  // Element face integral specification
  front_elem.clear(); back_elem.clear(); left_elem.clear(); right_elem.clear();
  top_elem.clear(); bottom_elem.clear();
  
  num_front_elem = 0; num_back_elem = 0;
  num_left_elem = 0; num_right_elem = 0;

  VEC_T::shrink2fit(front_elem); VEC_T::shrink2fit(back_elem);
  VEC_T::shrink2fit(left_elem); VEC_T::shrink2fit(right_elem);

  std::vector<int> front, back, left, right, top, bottom;
  Generate_BCElems_B(mesh->get_patch_mesh(0), front, back, left, right, top, bottom);

  bottom_elem.insert(bottom_elem.end(), bottom.begin(), bottom.end());

  Generate_BCElems_B(mesh->get_patch_mesh(1), front, back, left, right, top, bottom);

  top_elem.insert(top_elem.end(), top.begin(), top.end());
  
  Generate_BCElems_B(mesh->get_patch_mesh(2), front, back, left, right, top, bottom);

  bottom_elem.insert(bottom_elem.end(), bottom.begin(), bottom.end());

  VEC_T::shrink2fit(bottom_elem); VEC_T::shrink2fit(top_elem);

  num_top_elem = top_elem.size();
  num_bottom_elem = bottom_elem.size();

  std::cout<<"-----> 3-Patch Coronary geometry: NBC on patch 0 & 2 bottom, patch 1 top; ";
  std::cout<<" Front-Back C0 periodic; Left Dirichlet, Right nothing. \n";
}


void BoundaryCond::BC_type_103( const IMesh * const &mesh )
{
  const int numPat = mesh->get_num_patch();

  if(numPat != 3)
  {
    std::cerr<<"Error: BC_type_102 is designed for 3-patch coronary geometry only. \n";
    exit(EXIT_FAILURE);
  }

  // We extract the boundary nodes for each patch
  // Patch 0
  int nFunc_x_0, nFunc_y_0, nFunc_z_0, f_start_0;
  nFunc_x_0 = mesh->get_patch_mesh(0)->get_nFunc_x();
  nFunc_y_0 = mesh->get_patch_mesh(0)->get_nFunc_y();
  nFunc_z_0 = mesh->get_patch_mesh(0)->get_nFunc_z();
  f_start_0 = mesh->get_patch_mesh(0)->get_nFunc_start();
  
  int * back_front_0, * left_right_0, * bottom_top_0, * edge01_0, 
      * edge02_0, * edge04_0, * corner_0; 

  int nbf_0, nlr_0, nbt_0, nex_0, ney_0, nez_0;

  Generate_BCNodes_D(nFunc_x_0, nFunc_y_0, nFunc_z_0, f_start_0,
      back_front_0, nbf_0, left_right_0, nlr_0, bottom_top_0, nbt_0,
      edge01_0, nex_0, edge02_0, ney_0, edge04_0, nez_0, corner_0);
  
  // Patch 1
  int nFunc_x_1, nFunc_y_1, nFunc_z_1, f_start_1;
  nFunc_x_1 = mesh->get_patch_mesh(1)->get_nFunc_x();
  nFunc_y_1 = mesh->get_patch_mesh(1)->get_nFunc_y();
  nFunc_z_1 = mesh->get_patch_mesh(1)->get_nFunc_z();
  f_start_1 = mesh->get_patch_mesh(1)->get_nFunc_start();
 
  int * back_front_1, * left_right_1, * bottom_top_1, * edge01_1, 
      * edge02_1, * edge04_1, * corner_1; 

  int nbf_1, nlr_1, nbt_1, nex_1, ney_1, nez_1;

  Generate_BCNodes_D(nFunc_x_1, nFunc_y_1, nFunc_z_1, f_start_1,
      back_front_1, nbf_1, left_right_1, nlr_1, bottom_top_1, nbt_1,
      edge01_1, nex_1, edge02_1, ney_1, edge04_1, nez_1, corner_1);


  // Patch 2
  int nFunc_x_2, nFunc_y_2, nFunc_z_2, f_start_2;
  nFunc_x_2 = mesh->get_patch_mesh(2)->get_nFunc_x();
  nFunc_y_2 = mesh->get_patch_mesh(2)->get_nFunc_y();
  nFunc_z_2 = mesh->get_patch_mesh(2)->get_nFunc_z();
  f_start_2 = mesh->get_patch_mesh(2)->get_nFunc_start();
 
  int * back_front_2, * left_right_2, * bottom_top_2, * edge01_2, 
      * edge02_2, * edge04_2, * corner_2; 

  int nbf_2, nlr_2, nbt_2, nex_2, ney_2, nez_2;

  Generate_BCNodes_D(nFunc_x_2, nFunc_y_2, nFunc_z_2, f_start_2,
      back_front_2, nbf_2, left_right_2, nlr_2, bottom_top_2, nbt_2,
      edge01_2, nex_2, edge02_2, ney_2, edge04_2, nez_2, corner_2);


  // This is a 3-patch coronary geometry, with top/bottom faces of each patch as
  // the interface
  if( (nFunc_x_0 != nFunc_x_1) || (nFunc_y_0 != nFunc_y_1) )
  {
    std::cerr<<"Error: Patch 0 and Patch 1 interface function does not match. \n";
    exit(EXIT_FAILURE);
  }

  if( (nFunc_x_0 != nFunc_x_2) || (nFunc_y_0 != nFunc_y_2) )
  {
    std::cerr<<"Error: Patch 0 and Patch 2 interface function does not match. \n";
    exit(EXIT_FAILURE);
  }
  const int nFunc_x = nFunc_x_0;
  const int nFunc_y = nFunc_y_0;
  
  const int nFxFy = nFunc_x * nFunc_y;

  if( (nFunc_x-1) % 4 != 0 )
  {
    std::cerr<<"Error: nFunc_x = 4 (2+inserted_s) + 1. \n ";
    exit(EXIT_FAILURE);
  }

  // num_ins = number of basis functions between two C0 knots
  // The structure is
  // 0 , ... ,  2+num_ins, ..., 2(2+num_ins), ..., 3(2+num_ins), ... , 4(2+num_ins) 
  // e.g. num_ins = 0 (aka no knot insertion or degree elevation)
  // the basis in s direction has the following structure
  // 0,1,2,3,4,5,6,7,8 ===> 3   2   1
  //                        4       0/8
  //                        5   6   7
  const int num_ins = (nFunc_x-1) / 4 - 2;

  // ============= Interface Condition =============
  // The matching pattern is patch 0's top, patch 1's bottom and patch 2's top
  // meets at the interface
  // The first group of master-slave relationship : 
  // patch 0 : 2,6 on top is master,
  // patch 1 : 2,6 on bottom is slave
  // patch 2 : 2,6 on top is slave
  const int patch0_start = (nFunc_z_0-1) * nFxFy;
  const int patch1_start = f_start_1;
  const int patch2_start = (nFunc_z_2-1) * nFxFy + f_start_2;
  
  per_master_nodes.reserve(nFxFy*2 + nFunc_y * (nFunc_z_0 + nFunc_z_1 + nFunc_z_2) );
  per_slave_nodes.reserve(nFxFy*2 + nFunc_y * (nFunc_z_0 + nFunc_z_1 + nFunc_z_2) );
 
  for(int ii=1; ii<nFunc_y-1; ++ii)
  {
    // 2
    per_master_nodes.push_back( patch0_start + ii * nFunc_x + 2 + num_ins );
    per_slave_nodes.push_back(  patch1_start + ii * nFunc_x + 2 + num_ins );
    per_master_nodes.push_back( patch0_start + ii * nFunc_x + 2 + num_ins );
    per_slave_nodes.push_back(  patch2_start + ii * nFunc_x + 2 + num_ins );
    // 6
    per_master_nodes.push_back( patch0_start + ii * nFunc_x + 6 + 3*num_ins );
    per_slave_nodes.push_back(  patch1_start + ii * nFunc_x + 6 + 3*num_ins );
    per_master_nodes.push_back( patch0_start + ii * nFunc_x + 6 + 3*num_ins );
    per_slave_nodes.push_back(  patch2_start + ii * nFunc_x + 6 + 3*num_ins );
  } 
  
  for(int ii=1; ii<nFunc_y-1; ++ii)
  {
    const int level_t = ii * nFunc_x;
    // The second group of master-slave relationship :
    // patch 0 : 4 - 6 top
    // patch 1 : 4 - 6 top
    for( int jj= 3 + num_ins; jj<6+3*num_ins; ++jj )
    {
      per_master_nodes.push_back( patch0_start + level_t + jj );
      per_slave_nodes.push_back( patch1_start + level_t + jj );
    }
    // The third group of master-slave relationship
    // patch 0 : 1 top  --> patch 2 : 3 reverse order
    // patch 0 : 7 top  --> patch 2 : 5 reverse order
    for(int jj=1; jj<2+num_ins; ++jj)
    {
      per_master_nodes.push_back( patch0_start + level_t + jj);
      per_slave_nodes.push_back(  patch2_start + level_t + 4 + 2*num_ins - jj);    
      
      per_master_nodes.push_back( patch0_start + level_t + 6 + 3 * num_ins + jj);
      per_slave_nodes.push_back(  patch2_start + level_t + 6 + 3 * num_ins - jj);
    
      // The fourth group of master-slave relationship
      // patch 1 : 1 top --> patch 2 : 1 bottom
      // patch 1 : 7 top --> patch 2 : 7 bottom
      per_master_nodes.push_back( patch1_start + level_t + jj);
      per_slave_nodes.push_back(  patch2_start + level_t + jj );
      per_master_nodes.push_back( patch1_start + level_t + 6 + 3 * num_ins + jj);
      per_slave_nodes.push_back(  patch2_start + level_t + 6 + 3 * num_ins + jj );
    }
    // The fifth group of master-slave relationship
    // patch 0 : 0 top --> patch 0 : 8 top
    // patch 0 : 0 top --> patch 2 : 5 bottom
    per_master_nodes.push_back(patch0_start + level_t);
    per_slave_nodes.push_back( patch0_start + level_t + 8 + 4*num_ins);
    per_master_nodes.push_back(patch0_start + level_t);
    per_slave_nodes.push_back( patch2_start + level_t + 4 + 2 * num_ins);
    
    // The sixth group of master-slave relationship
    // patch 1 : 0 top --> patch 1 : 8 top
    // patch 1 : 0 top --> patch 2 : 0 bottom
    // patch 1 : 0 top --> patch 2 : 8 bottom
    per_master_nodes.push_back(patch1_start + level_t);
    per_slave_nodes.push_back(patch1_start + level_t + 8 + 4*num_ins);
    per_master_nodes.push_back(patch1_start + level_t);
    per_slave_nodes.push_back(patch2_start + level_t );
    per_master_nodes.push_back(patch1_start + level_t);
    per_slave_nodes.push_back(patch2_start + level_t + 8 + 4 * num_ins);
  }
  

  // Periodic paris within each patch
  // patch 0 : front - back, e15 - e04, c5-c4
  per_master_nodes.insert(per_master_nodes.end(), back_front_0 + nbf_0, back_front_0 + 2*nbf_0);
  per_slave_nodes.insert(per_slave_nodes.end(), back_front_0, back_front_0 + nbf_0);
  
  per_master_nodes.insert(per_master_nodes.end(), edge04_0+1*nez_0, edge04_0+2*nez_0);
  per_slave_nodes.insert(per_slave_nodes.end(), edge04_0+0*nez_0, edge04_0+1*nez_0);

  per_master_nodes.push_back(corner_0[5]);
  per_slave_nodes.push_back(corner_0[4]);

  // patch 1 : front - back, e15 - e04, e57 - e46, c5-c4, c1-c0
  per_master_nodes.insert(per_master_nodes.end(), back_front_1 + nbf_1, back_front_1 + 2*nbf_1);
  per_slave_nodes.insert(per_slave_nodes.end(), back_front_1, back_front_1 + nbf_1);
  
  per_master_nodes.insert(per_master_nodes.end(), edge04_1+1*nez_1, edge04_1+2*nez_1);
  per_slave_nodes.insert(per_slave_nodes.end(), edge04_1+0*nez_1, edge04_1+1*nez_1);

  per_master_nodes.insert(per_master_nodes.end(), edge02_1+3*ney_1, edge02_1+4*ney_1);
  per_slave_nodes.insert(per_slave_nodes.end(), edge02_1+2*ney_1, edge02_1+3*ney_1);
  
  per_master_nodes.push_back(corner_1[5]);
  per_slave_nodes.push_back(corner_1[4]); 
  
  per_master_nodes.push_back(corner_1[1]);
  per_slave_nodes.push_back(corner_1[0]); 
  
  // patch 2 : front - back, e15 - e04, e13 - e02, c1-c0, c7-c6
  per_master_nodes.insert(per_master_nodes.end(), back_front_2 + nbf_2, back_front_2 + 2*nbf_2);
  per_slave_nodes.insert(per_slave_nodes.end(), back_front_2, back_front_2 + nbf_2);
  
  per_master_nodes.insert(per_master_nodes.end(), edge04_2+1*nez_2, edge04_2+2*nez_2);
  per_slave_nodes.insert(per_slave_nodes.end(), edge04_2+0*nez_2, edge04_2+1*nez_2);

  per_master_nodes.insert(per_master_nodes.end(), edge02_2+1*ney_2, edge02_2+2*ney_2);
  per_slave_nodes.insert( per_slave_nodes.end(),  edge02_2+0*ney_2, edge02_2+1*ney_2);
  
  per_master_nodes.push_back(corner_2[1]);
  per_slave_nodes.push_back(corner_2[0]); 
  
  per_master_nodes.push_back(corner_2[5]);
  per_slave_nodes.push_back(corner_2[4]); 

  VEC_T::shrink2fit(per_master_nodes);
  VEC_T::shrink2fit(per_slave_nodes);


  // Dirichlet nodes
  dir_nodes.reserve(nlr_0 + nlr_1 + nlr_2 + 2*nex_0 + 2*nex_1 + 2*nex_2 
      + 2*nez_0 + 2*nez_1 + 2*nez_2);

  // patch 0 : f_right, f_bottom, e01, e23, e67, e26, e37, e02, e13, c0,1,2,3,6,7
  dir_nodes.insert(dir_nodes.end(), left_right_0 + nlr_0, left_right_0 + 2*nlr_0);
  dir_nodes.insert(dir_nodes.end(), bottom_top_0, bottom_top_0 + nbt_0);
  dir_nodes.insert(dir_nodes.end(), edge01_0, edge01_0 + 2*nex_0 );
  dir_nodes.insert(dir_nodes.end(), edge01_0 + 3*nex_0, edge01_0 + 4*nex_0 );
  dir_nodes.insert(dir_nodes.end(), edge04_0 + 2*nez_0, edge04_0 + 4*nez_0 );
  dir_nodes.insert(dir_nodes.end(), edge02_0, edge02_0 + 2*ney_0 );
  dir_nodes.push_back(corner_0[0]);
  dir_nodes.push_back(corner_0[1]);
  dir_nodes.push_back(corner_0[2]);
  dir_nodes.push_back(corner_0[3]);
  dir_nodes.push_back(corner_0[6]);
  dir_nodes.push_back(corner_0[7]);

  // patch 1 : f_right, e23, e67, e26, e37, c2,3,6,7
  dir_nodes.insert(dir_nodes.end(), left_right_1 + nlr_1, left_right_1 + 2*nlr_1);
  dir_nodes.insert(dir_nodes.end(), edge01_1 + nex_1, edge01_1 + 2*nex_1 ); // e23
  dir_nodes.insert(dir_nodes.end(), edge01_1 + 3*nex_1, edge01_1 + 4*nex_1 ); // e67
  dir_nodes.insert(dir_nodes.end(), edge04_1 + 2*nez_1, edge04_1 + 4*nez_1 ); // e26 e37
  dir_nodes.push_back(corner_1[2]);
  dir_nodes.push_back(corner_1[3]);
  dir_nodes.push_back(corner_1[6]);
  dir_nodes.push_back(corner_1[7]);

  // patch 2 : f_right, e23, e67, e26, e37, c2,3,6,7
  dir_nodes.insert(dir_nodes.end(), left_right_2 + nlr_2, left_right_2 + 2*nlr_2);
  dir_nodes.insert(dir_nodes.end(), edge01_2 + nex_2, edge01_2 + 2*nex_2 );
  dir_nodes.insert(dir_nodes.end(), edge01_2 + 3*nex_2, edge01_2 + 4*nex_2 );
  dir_nodes.insert(dir_nodes.end(), edge04_2 + 2*nez_2, edge04_2 + 4*nez_2 );
  dir_nodes.push_back(corner_2[2]);
  dir_nodes.push_back(corner_2[3]);
  dir_nodes.push_back(corner_2[6]);
  dir_nodes.push_back(corner_2[7]);
  
  VEC_T::shrink2fit(dir_nodes);

  delete [] back_front_0; delete [] left_right_0; delete [] bottom_top_0;
  delete [] edge01_0; delete [] edge02_0; delete [] edge04_0; delete [] corner_0;
  delete [] back_front_1; delete [] left_right_1; delete [] bottom_top_1;
  delete [] edge01_1; delete [] edge02_1; delete [] edge04_1; delete [] corner_1;
  delete [] back_front_2; delete [] left_right_2; delete [] bottom_top_2;
  delete [] edge01_2; delete [] edge02_2; delete [] edge04_2; delete [] corner_2;

  num_per_nodes = per_master_nodes.size();
  num_dir_nodes = dir_nodes.size();

  // Element face integral specification
  front_elem.clear(); back_elem.clear(); left_elem.clear(); right_elem.clear();
  top_elem.clear(); bottom_elem.clear();
  
  num_front_elem = 0; num_back_elem = 0;
  num_left_elem = 0; num_right_elem = 0;

  VEC_T::shrink2fit(front_elem); VEC_T::shrink2fit(back_elem);
  VEC_T::shrink2fit(left_elem); VEC_T::shrink2fit(right_elem);

  std::vector<int> front, back, left, right, top, bottom;
  
  Generate_BCElems_B(mesh->get_patch_mesh(1), front, back, left, right, top, bottom);
  top_elem.insert(top_elem.end(), top.begin(), top.end());
  
  Generate_BCElems_B(mesh->get_patch_mesh(2), front, back, left, right, top, bottom);
  bottom_elem.insert(bottom_elem.end(), bottom.begin(), bottom.end());

  VEC_T::shrink2fit(bottom_elem); VEC_T::shrink2fit(top_elem);

  num_top_elem = top_elem.size();
  num_bottom_elem = bottom_elem.size();

  std::cout<<"-----> 3-Patch Coronary geometry: NBC on patch 2 bottom, patch 1 top; ";
  std::cout<<" Dirichlet on patch 0 bottom. ";
  std::cout<<" Front-Back C0 periodic; Left Nothing, Right Dirichlet. \n";
}

void BoundaryCond::Read_SH_interface_pair( const std::string &interface_file,
    std::vector<int> &patch1, std::vector<int> &patch2,
    std::vector<int> &face1, std::vector<int> &face2 ) const
{
  patch1.clear(); patch2.clear();
  face1.clear(); face2.clear();
  
  std::ifstream infile( interface_file.c_str(), std::ifstream::in );
  if( !infile.is_open() )
  {
    std::cerr<<"Error: Cannot open the file: "<<interface_file<<std::endl;
    exit(EXIT_FAILURE);
  }

  std::istringstream sstrm;
  std::string sline;

  getline(infile, sline);
  sstrm.str(sline);
  int num_intface;
  sstrm>>num_intface;
  sstrm.clear();
  
  int p1, p2, f1, f2;
  for(int ii=0; ii<num_intface; ++ii)
  {
    getline(infile, sline);
    sstrm.str(sline);
    sstrm>>p1;
    sstrm>>p2;
    sstrm>>f1;
    sstrm>>f2;
    sstrm.clear();
    p1 = p1 - 1;
    p2 = p2 - 1;

    patch1.push_back(p1);
    patch2.push_back(p2);
    face1.push_back(f1);
    face2.push_back(f2);
  }

}


void BoundaryCond::BC_type_104( const IMesh * const &mesh, const std::vector<double> &cPts )
{
  const int numPat = mesh->get_num_patch();
  
  const double tol = 1.0e-6;

  if(numPat != 3)
  {
    std::cerr<<"Error: BC_type_104 is designed for 3-patch geometry only. \n";
    exit(EXIT_FAILURE);
  }

  char * char_home_dir = getenv("HOME");
  std::string interface_file(char_home_dir);
  interface_file.append("/PERIGEE/input/coronary_3patch/global_data.dat");

  // Two patches intpat1[i] and intpat2[i] meet with intpat1's intfac1[i] and
  // intpat2's intfac2[i]. 
  std::vector<int> intpat1, intpat2, intfac1, intfac2;

  // Read in the values for the four int std::vectors from data file to specify all
  // interior interfaces.
  // Note: in shaolie's numbering, first patch starts with index 1. In the
  // reader, we correct it to start with 0.
  Read_SH_interface_pair( interface_file, intpat1, intpat2, intfac1, intfac2 );
  
  std::cout<<"Read interface interface pair from "<<interface_file<<std::endl;


  const int nFunc_x_0 = mesh->get_patch_mesh(0)->get_nFunc_x();
  const int nFunc_y_0 = mesh->get_patch_mesh(0)->get_nFunc_y();

  
  per_master_nodes.clear();
  per_slave_nodes.clear();
  
  dir_nodes.clear();

  // Now we loop over patches
  for(int ii=0; ii<numPat; ++ii)
  {
    int nFunc_x = mesh->get_patch_mesh(ii)->get_nFunc_x();
    int nFunc_y = mesh->get_patch_mesh(ii)->get_nFunc_y();
    int nFunc_z = mesh->get_patch_mesh(ii)->get_nFunc_z();
    int f_start = mesh->get_patch_mesh(ii)->get_nFunc_start();

    if( nFunc_x != nFunc_x_0 || nFunc_y != nFunc_y_0 )
    {
      std::cerr<<"Error: number of basis functions in x/y-direction does not match for patch "
        <<ii<<" !"<<std::endl;
      exit(EXIT_FAILURE);
    }

    int * back_front, * left_right, * bottom_top, * edge01, * edge02, * edge04, * corner;
    int nbf, nlr, nbt, nex, ney, nez;

    // Generate the indices on all boundary faces, edges, points.
    Generate_BCNodes_D(nFunc_x, nFunc_y, nFunc_z, f_start,
        back_front, nbf, left_right, nlr, bottom_top, nbt,
        edge01, nex, edge02, ney, edge04, nez, corner);

    // Impose periodic pairs for front - back, 
    // e15 - e04, e13 - e02, e37 - e26, e57 - e46
    // c1 - c0, c5 - c4, c3 - c2, c7 - c6
    per_master_nodes.insert( per_master_nodes.end(), back_front+nbf, back_front+2*nbf );
    per_slave_nodes.insert( per_slave_nodes.end(), back_front, back_front + nbf );

    per_master_nodes.insert( per_master_nodes.end(), edge04 + nez, edge04 + 2*nez );
    per_slave_nodes.insert( per_slave_nodes.end(), edge04, edge04 + nez );
    
    per_master_nodes.insert( per_master_nodes.end(), edge02+ney, edge02 + 2*ney );
    per_slave_nodes.insert( per_slave_nodes.end(), edge02, edge02 + ney );

    per_master_nodes.insert( per_master_nodes.end(), edge04+3*nez, edge04+4*nez );
    per_slave_nodes.insert( per_slave_nodes.end(), edge04+2*nez, edge04+3*nez );

    per_master_nodes.insert( per_master_nodes.end(), edge02+3*ney, edge02+4*ney );
    per_slave_nodes.insert( per_slave_nodes.end(), edge02+2*ney, edge02+3*ney);

    per_master_nodes.push_back(corner[5]);
    per_slave_nodes.push_back(corner[4]);
    
    per_master_nodes.push_back(corner[1]);
    per_slave_nodes.push_back(corner[0]);
    
    per_master_nodes.push_back(corner[7]);
    per_slave_nodes.push_back(corner[6]);
    
    per_master_nodes.push_back(corner[3]);
    per_slave_nodes.push_back(corner[2]);

    // Impose Dirichlet nodes on Right, e23, e67, e26, e37, c2,3,6,7
    // Note: There are nodes that are both per and dir. It is OK for now. It
    // will be handled in the creation of ID array. We will specify Dir nodes in
    // the last step in ID generation, so the repeated nodes will have Dirichlet
    // condition only.
    dir_nodes.insert( dir_nodes.end(), left_right + nlr, left_right + 2*nlr );
    dir_nodes.insert( dir_nodes.end(), edge01 + nex, edge01 + 2*nex );
    dir_nodes.insert( dir_nodes.end(), edge01 + 3*nex, edge01 + 4*nex );
    dir_nodes.insert( dir_nodes.end(), edge04 + 2*nez, edge04 + 4*nez );
    dir_nodes.push_back(corner[2]);
    dir_nodes.push_back(corner[3]);
    dir_nodes.push_back(corner[6]);
    dir_nodes.push_back(corner[7]);

    delete [] back_front; delete [] left_right; delete [] bottom_top;
    delete [] edge01; delete [] edge02; delete [] edge04; delete [] corner;
  }

  // Now we create a pesudoID array
  std::vector<int> pID;
  pID.clear();
  pID.resize(mesh->get_nFunc());
  for(unsigned int ii=0; ii<pID.size(); ++ii)
    pID[ii] = ii;
  for(unsigned int ii=0; ii<per_master_nodes.size(); ++ii)
    pID[per_slave_nodes[ii]] = per_master_nodes[ii];


  // Now we loop over each interior pair
  int patchm, patchs, facem, faces, patchm_start, patchs_start;
  int master_index, slave_index;
  double mx, my, mz, sx, sy, sz;
  double dist;
  for(unsigned int ii=0; ii<intpat1.size(); ++ii)
  {
    // we set the patch with smaller index as the slave; with bigger index as
    // the master. If two patch has identical index, report error.
    if( intpat1[ii] > intpat2[ii] )
    {
      patchm = intpat1[ii]; facem = intfac1[ii];
      patchs = intpat2[ii]; faces = intfac2[ii];
    }
    else if( intpat1[ii] < intpat2[ii] )
    {
      patchm = intpat2[ii]; facem = intfac2[ii];
      patchs = intpat1[ii]; faces = intfac1[ii];
    }
    else
    {
      std::cerr<<"Error: the "<<ii<<"th interior interface pair involves patches "<<
        intpat1[ii]<<" and "<<intpat2[ii]<<std::endl;
      exit(EXIT_FAILURE);
    }
   
    if(patchm >= numPat || patchs >= numPat)
    {
      std::cerr<<"Error: the patches have index "<<patchm<<'\t'<<patchs
        <<". This BC is for "<<numPat<<" patch geometry. \n";
      exit(EXIT_FAILURE);
    }

    // Now we determine the starting index for the top/bottom interface.
    int nFunc_z_m = mesh->get_patch_mesh(patchm)->get_nFunc_z();
    int nFunc_z_s = mesh->get_patch_mesh(patchs)->get_nFunc_z();
    int f_start_m = mesh->get_patch_mesh(patchm)->get_nFunc_start();
    int f_start_s = mesh->get_patch_mesh(patchs)->get_nFunc_start(); 

    if( facem == 1)
      patchm_start = f_start_m;
    else if( facem == 6 )
      patchm_start = f_start_m + (nFunc_z_m-1) * nFunc_x_0 * nFunc_y_0;
    else
    {
      std::cerr<<"Error: the "<<ii
        <<"th interior master interface has non-top/non-bottom face involved. \n";
      exit(EXIT_FAILURE);
    } 

    if( faces == 1 )
      patchs_start = f_start_s;
    else if( faces == 6 )
      patchs_start = f_start_s + (nFunc_z_s-1) * nFunc_x_0 * nFunc_y_0;
    else
    {
      std::cerr<<"Error: the "<<ii
        <<"th interior slave interface has non-top/non-bottom face involved. \n";
      exit(EXIT_FAILURE);
    } 

    int counter = 0;
    for(int jj=0; jj<nFunc_y_0; ++jj)
    {
      for(int kk=0; kk<nFunc_x_0; ++kk)
      {
        master_index = patchm_start + jj * nFunc_x_0 + kk;
        mx = cPts[master_index * 4 + 0];
        my = cPts[master_index * 4 + 1];
        mz = cPts[master_index * 4 + 2];
        
        for(int nn=0; nn<nFunc_y_0; ++nn)
        {
          for(int mm=0; mm<nFunc_x_0; ++mm)
          {
            slave_index = patchs_start + nn * nFunc_x_0 + mm;
            sx = cPts[slave_index * 4 + 0];
            sy = cPts[slave_index * 4 + 1];
            sz = cPts[slave_index * 4 + 2];
            
            dist = sqrt( (sx - mx) * (sx - mx) + (sy - my) * (sy - my)
               + (sz - mz) * (sz - mz) ); 
            
            if( dist < tol )
            {
              counter += 1;
              
              if( pID[slave_index] < master_index )
                pID[slave_index] = master_index;
            }
          }
        }
      }
    }

  } // Finish the loop for interior interfaces

  for(unsigned int ii=0; ii<dir_nodes.size(); ++ii)
    pID[dir_nodes[ii]] = -1;


  per_master_nodes.clear();
  per_slave_nodes.clear();
  for(int ii=0; ii<(int) pID.size(); ++ii)
  {
    if(pID[ii] > 0)
    {
      if(pID[ii] > ii)
      {
        per_master_nodes.push_back(pID[ii]);
        per_slave_nodes.push_back(ii);
      }
      else if(pID[ii] < ii)
      {
        std::cerr<<"Error: pseudoID cannot be smaller that its index. \n";
        exit(EXIT_FAILURE);
      }
    }
  }

  num_per_nodes = per_master_nodes.size();
  num_dir_nodes = dir_nodes.size();


  // Now element face integral specification for NBC
  front_elem.clear(); back_elem.clear(); left_elem.clear(); right_elem.clear();
  top_elem.clear(); bottom_elem.clear();

  num_front_elem = 0; num_back_elem = 0;
  num_left_elem = 0; num_right_elem = 0;

  VEC_T::shrink2fit(front_elem); VEC_T::shrink2fit(back_elem);
  VEC_T::shrink2fit(left_elem); VEC_T::shrink2fit(right_elem);

  std::vector<int> front, back, left, right, top, bottom;
  Generate_BCElems_B(mesh->get_patch_mesh(0), front, back, left, right, top, bottom);

  bottom_elem.insert(bottom_elem.end(), bottom.begin(), bottom.end());

  Generate_BCElems_B(mesh->get_patch_mesh(1), front, back, left, right, top, bottom);

  top_elem.insert(top_elem.end(), top.begin(), top.end());

  Generate_BCElems_B(mesh->get_patch_mesh(2), front, back, left, right, top, bottom);

  bottom_elem.insert(bottom_elem.end(), bottom.begin(), bottom.end());

  VEC_T::shrink2fit(bottom_elem); VEC_T::shrink2fit(top_elem);

  num_top_elem = top_elem.size();
  num_bottom_elem = bottom_elem.size();

  std::cout<<"-----> 3-Patch Coronary geometry: NBC on patch 0 & 2 bottom, patch 1 top; ";
  std::cout<<" Front-Back C0 periodic; Right Dirichlet, Left nothing. \n";


}


void BoundaryCond::BC_type_105( const IMesh * const &mesh, const std::vector<double> &cPts )
{
  const int numPat = mesh->get_num_patch();
  
  const double tol = 1.0e-6;

  if(numPat != 3)
  {
    std::cerr<<"Error: BC_type_105 is designed for 3-patch coronary geometry only. \n";
    exit(EXIT_FAILURE);
  }

  char * char_home_dir = getenv("HOME");
  std::string interface_file(char_home_dir);
  interface_file.append("/PERIGEE/input/coronary_3patch/global_data.dat");

  // Two patches intpat1[i] and intpat2[i] meet with intpat1's intfac1[i] and
  // intpat2's intfac2[i]. 
  std::vector<int> intpat1, intpat2, intfac1, intfac2;

  // Read in the values for the four int std::vectors from data file to specify all
  // interior interfaces.
  // Note: in shaolie's numbering, first patch starts with index 1. In the
  // reader, we correct it to start with 0.
  Read_SH_interface_pair( interface_file, intpat1, intpat2, intfac1, intfac2 );
  
  std::cout<<"Read interface interface pair from "<<interface_file<<std::endl;


  const int nFunc_x_0 = mesh->get_patch_mesh(0)->get_nFunc_x();
  const int nFunc_y_0 = mesh->get_patch_mesh(0)->get_nFunc_y();

  
  per_master_nodes.clear();
  per_slave_nodes.clear();
  
  dir_nodes.clear();

  // Now we loop over patches
  for(int ii=0; ii<numPat; ++ii)
  {
    int nFunc_x = mesh->get_patch_mesh(ii)->get_nFunc_x();
    int nFunc_y = mesh->get_patch_mesh(ii)->get_nFunc_y();
    int nFunc_z = mesh->get_patch_mesh(ii)->get_nFunc_z();
    int f_start = mesh->get_patch_mesh(ii)->get_nFunc_start();

    if( nFunc_x != nFunc_x_0 || nFunc_y != nFunc_y_0 )
    {
      std::cerr<<"Error: number of basis functions in x/y-direction does not match for patch "
        <<ii<<" !"<<std::endl;
      exit(EXIT_FAILURE);
    }

    int * back_front, * left_right, * bottom_top, * edge01, * edge02, * edge04, * corner;
    int nbf, nlr, nbt, nex, ney, nez;

    // Generate the indices on all boundary faces, edges, points.
    Generate_BCNodes_D(nFunc_x, nFunc_y, nFunc_z, f_start,
        back_front, nbf, left_right, nlr, bottom_top, nbt,
        edge01, nex, edge02, ney, edge04, nez, corner);

    // Impose periodic pairs for front - back, 
    // e15 - e04, e13 - e02, e37 - e26, e57 - e46
    // c1 - c0, c5 - c4, c3 - c2, c7 - c6
    per_master_nodes.insert( per_master_nodes.end(), back_front+nbf, back_front+2*nbf );
    per_slave_nodes.insert( per_slave_nodes.end(), back_front, back_front + nbf );

    per_master_nodes.insert( per_master_nodes.end(), edge04 + nez, edge04 + 2*nez );
    per_slave_nodes.insert( per_slave_nodes.end(), edge04, edge04 + nez );
    
    per_master_nodes.insert( per_master_nodes.end(), edge02+ney, edge02 + 2*ney );
    per_slave_nodes.insert( per_slave_nodes.end(), edge02, edge02 + ney );

    per_master_nodes.insert( per_master_nodes.end(), edge04+3*nez, edge04+4*nez );
    per_slave_nodes.insert( per_slave_nodes.end(), edge04+2*nez, edge04+3*nez );

    per_master_nodes.insert( per_master_nodes.end(), edge02+3*ney, edge02+4*ney );
    per_slave_nodes.insert( per_slave_nodes.end(), edge02+2*ney, edge02+3*ney);

    per_master_nodes.push_back(corner[5]);
    per_slave_nodes.push_back(corner[4]);
    
    per_master_nodes.push_back(corner[1]);
    per_slave_nodes.push_back(corner[0]);
    
    per_master_nodes.push_back(corner[7]);
    per_slave_nodes.push_back(corner[6]);
    
    per_master_nodes.push_back(corner[3]);
    per_slave_nodes.push_back(corner[2]);


    // Mater-slave nodes on LEFT face, e01, e45
    for( int zz = 0; zz < nFunc_z; ++zz )
    {
      int mnode = zz * nFunc_x * nFunc_y + nFunc_x - 1;
      for( int xx = 1; xx < nFunc_x - 1; ++xx )
      {
        int snode = zz * nFunc_x * nFunc_y + xx;

        per_master_nodes.push_back(mnode);
        per_slave_nodes.push_back(snode);
      }
    }

    // Impose Dirichlet nodes on Right, e23, e67, e26, e37, c2,3,6,7
    // Note: There are nodes that are both per and dir. It is OK for now. It
    // will be handled in the creation of ID array. We will specify Dir nodes in
    // the last step in ID generation, so the repeated nodes will have Dirichlet
    // condition only.
    dir_nodes.insert( dir_nodes.end(), left_right + nlr, left_right + 2*nlr );
    dir_nodes.insert( dir_nodes.end(), edge01 + nex, edge01 + 2*nex );
    dir_nodes.insert( dir_nodes.end(), edge01 + 3*nex, edge01 + 4*nex );
    dir_nodes.insert( dir_nodes.end(), edge04 + 2*nez, edge04 + 4*nez );
    dir_nodes.push_back(corner[2]);
    dir_nodes.push_back(corner[3]);
    dir_nodes.push_back(corner[6]);
    dir_nodes.push_back(corner[7]);

    delete [] back_front; delete [] left_right; delete [] bottom_top;
    delete [] edge01; delete [] edge02; delete [] edge04; delete [] corner;
  }

  // Now we create a pesudoID array
  std::vector<int> pID;
  pID.clear();
  pID.resize(mesh->get_nFunc());
  for(unsigned int ii=0; ii<pID.size(); ++ii)
    pID[ii] = ii;
  for(unsigned int ii=0; ii<per_master_nodes.size(); ++ii)
    pID[per_slave_nodes[ii]] = per_master_nodes[ii];


  // Now we loop over each interior pair
  int patchm, patchs, facem, faces, patchm_start, patchs_start;
  int master_index, slave_index;
  double mx, my, mz, sx, sy, sz;
  double dist;
  for(unsigned int ii=0; ii<intpat1.size(); ++ii)
  {
    // we set the patch with smaller index as the slave; with bigger index as
    // the master. If two patch has identical index, report error.
    if( intpat1[ii] > intpat2[ii] )
    {
      patchm = intpat1[ii]; facem = intfac1[ii];
      patchs = intpat2[ii]; faces = intfac2[ii];
    }
    else if( intpat1[ii] < intpat2[ii] )
    {
      patchm = intpat2[ii]; facem = intfac2[ii];
      patchs = intpat1[ii]; faces = intfac1[ii];
    }
    else
    {
      std::cerr<<"Error: the "<<ii<<"th interior interface pair involves patches "<<
        intpat1[ii]<<" and "<<intpat2[ii]<<std::endl;
      exit(EXIT_FAILURE);
    }
   
    if(patchm >= numPat || patchs >= numPat)
    {
      std::cerr<<"Error: the patches have index "<<patchm<<'\t'<<patchs
        <<". This BC is for "<<numPat<<" patch geometry. \n";
      exit(EXIT_FAILURE);
    }

    // Now we determine the starting index for the top/bottom interface.
    int nFunc_z_m = mesh->get_patch_mesh(patchm)->get_nFunc_z();
    int nFunc_z_s = mesh->get_patch_mesh(patchs)->get_nFunc_z();
    int f_start_m = mesh->get_patch_mesh(patchm)->get_nFunc_start();
    int f_start_s = mesh->get_patch_mesh(patchs)->get_nFunc_start(); 

    if( facem == 1)
      patchm_start = f_start_m;
    else if( facem == 6 )
      patchm_start = f_start_m + (nFunc_z_m-1) * nFunc_x_0 * nFunc_y_0;
    else
    {
      std::cerr<<"Error: the "<<ii
        <<"th interior master interface has non-top/non-bottom face involved. \n";
      exit(EXIT_FAILURE);
    } 

    if( faces == 1 )
      patchs_start = f_start_s;
    else if( faces == 6 )
      patchs_start = f_start_s + (nFunc_z_s-1) * nFunc_x_0 * nFunc_y_0;
    else
    {
      std::cerr<<"Error: the "<<ii
        <<"th interior slave interface has non-top/non-bottom face involved. \n";
      exit(EXIT_FAILURE);
    } 

    int counter = 0;
    for(int jj=0; jj<nFunc_y_0; ++jj)
    {
      for(int kk=0; kk<nFunc_x_0; ++kk)
      {
        master_index = patchm_start + jj * nFunc_x_0 + kk;
        mx = cPts[master_index * 4 + 0];
        my = cPts[master_index * 4 + 1];
        mz = cPts[master_index * 4 + 2];
        
        for(int nn=0; nn<nFunc_y_0; ++nn)
        {
          for(int mm=0; mm<nFunc_x_0; ++mm)
          {
            slave_index = patchs_start + nn * nFunc_x_0 + mm;
            sx = cPts[slave_index * 4 + 0];
            sy = cPts[slave_index * 4 + 1];
            sz = cPts[slave_index * 4 + 2];
            
            dist = sqrt( (sx - mx) * (sx - mx) + (sy - my) * (sy - my)
               + (sz - mz) * (sz - mz) ); 
            
            if( dist < tol )
            {
              counter += 1;
              
              if( pID[slave_index] < master_index )
                pID[slave_index] = master_index;
            }
          }
        }
      }
    }

  } // Finish the loop for interior interfaces

  for(unsigned int ii=0; ii<dir_nodes.size(); ++ii)
    pID[dir_nodes[ii]] = -1;


  per_master_nodes.clear();
  per_slave_nodes.clear();
  for(int ii=0; ii<(int) pID.size(); ++ii)
  {
    if(pID[ii] > 0)
    {
      if(pID[ii] > ii)
      {
        per_master_nodes.push_back(pID[ii]);
        per_slave_nodes.push_back(ii);
      }
      else if(pID[ii] < ii)
      {
        std::cerr<<"Error: pseudoID cannot be smaller that its index. \n";
        exit(EXIT_FAILURE);
      }
    }
  }

  num_per_nodes = per_master_nodes.size();
  num_dir_nodes = dir_nodes.size();


  // Now element face integral specification for NBC
  front_elem.clear(); back_elem.clear(); left_elem.clear(); right_elem.clear();
  top_elem.clear(); bottom_elem.clear();

  num_front_elem = 0; num_back_elem = 0;
  num_left_elem = 0; num_right_elem = 0;

  VEC_T::shrink2fit(front_elem); VEC_T::shrink2fit(back_elem);
  VEC_T::shrink2fit(left_elem); VEC_T::shrink2fit(right_elem);

  std::vector<int> front, back, left, right, top, bottom;
  Generate_BCElems_B(mesh->get_patch_mesh(0), front, back, left, right, top, bottom);

  bottom_elem.insert(bottom_elem.end(), bottom.begin(), bottom.end());

  Generate_BCElems_B(mesh->get_patch_mesh(1), front, back, left, right, top, bottom);

  top_elem.insert(top_elem.end(), top.begin(), top.end());

  Generate_BCElems_B(mesh->get_patch_mesh(2), front, back, left, right, top, bottom);

  bottom_elem.insert(bottom_elem.end(), bottom.begin(), bottom.end());

  VEC_T::shrink2fit(bottom_elem); VEC_T::shrink2fit(top_elem);

  num_top_elem = top_elem.size();
  num_bottom_elem = bottom_elem.size();

  std::cout<<"-----> 3-Patch Coronary geometry: NBC on patch 0 & 2 bottom, patch 1 top; ";
  std::cout<<" Front-Back C0 periodic; Right Dirichlet, Left face follows edge15. \n";

}


void BoundaryCond::BC_type_200( const IMesh * const &mesh, const std::vector<double> &cPts )
{
  const int numPat = mesh->get_num_patch();
  
  const double tol = 1.0e-6;

  if(numPat != 107)
  {
    std::cerr<<"Error: BC_type_200 is designed for 107-patch geometry only. \n";
    exit(EXIT_FAILURE);
  }

  char * char_home_dir = getenv("HOME");
  std::string interface_file(char_home_dir);
  interface_file.append("/PERIGEE/input/pig_coronary_107patch/global_data.dat");

  // Two patches intpat1[i] and intpat2[i] meet with intpat1's intfac1[i] and
  // intpat2's intfac2[i]. 
  std::vector<int> intpat1, intpat2, intfac1, intfac2;

  // Read in the values for the four int std::vectors from data file to specify all
  // interior interfaces.
  // Note: in shaolie's numbering, first patch starts with index 1. In the
  // reader, we correct it to start with 0.
  Read_SH_interface_pair( interface_file, intpat1, intpat2, intfac1, intfac2 );
  
  std::cout<<"Read interface interface pair from "<<interface_file<<std::endl;


  const int nFunc_x_0 = mesh->get_patch_mesh(0)->get_nFunc_x();
  const int nFunc_y_0 = mesh->get_patch_mesh(0)->get_nFunc_y();

  
  per_master_nodes.clear();
  per_slave_nodes.clear();
  
  dir_nodes.clear();

  // Now we loop over patches
  for(int ii=0; ii<numPat; ++ii)
  {
    int nFunc_x = mesh->get_patch_mesh(ii)->get_nFunc_x();
    int nFunc_y = mesh->get_patch_mesh(ii)->get_nFunc_y();
    int nFunc_z = mesh->get_patch_mesh(ii)->get_nFunc_z();
    int f_start = mesh->get_patch_mesh(ii)->get_nFunc_start();

    if( nFunc_x != nFunc_x_0 || nFunc_y != nFunc_y_0 )
    {
      std::cerr<<"Error: number of basis functions in x/y-direction does not match for patch "
        <<ii<<" !"<<std::endl;
      exit(EXIT_FAILURE);
    }

    int * back_front, * left_right, * bottom_top, * edge01, * edge02, * edge04, * corner;
    int nbf, nlr, nbt, nex, ney, nez;

    // Generate the indices on all boundary faces, edges, points.
    Generate_BCNodes_D(nFunc_x, nFunc_y, nFunc_z, f_start,
        back_front, nbf, left_right, nlr, bottom_top, nbt,
        edge01, nex, edge02, ney, edge04, nez, corner);

    // Impose periodic pairs for front - back, 
    // e15 - e04, e13 - e02, e37 - e26, e57 - e46
    // c1 - c0, c5 - c4, c3 - c2, c7 - c6
    per_master_nodes.insert( per_master_nodes.end(), back_front+nbf, back_front+2*nbf );
    per_slave_nodes.insert( per_slave_nodes.end(), back_front, back_front + nbf );

    per_master_nodes.insert( per_master_nodes.end(), edge04 + nez, edge04 + 2*nez );
    per_slave_nodes.insert( per_slave_nodes.end(), edge04, edge04 + nez );
    
    per_master_nodes.insert( per_master_nodes.end(), edge02+ney, edge02 + 2*ney );
    per_slave_nodes.insert( per_slave_nodes.end(), edge02, edge02 + ney );

    per_master_nodes.insert( per_master_nodes.end(), edge04+3*nez, edge04+4*nez );
    per_slave_nodes.insert( per_slave_nodes.end(), edge04+2*nez, edge04+3*nez );

    per_master_nodes.insert( per_master_nodes.end(), edge02+3*ney, edge02+4*ney );
    per_slave_nodes.insert( per_slave_nodes.end(), edge02+2*ney, edge02+3*ney);

    per_master_nodes.push_back(corner[5]);
    per_slave_nodes.push_back(corner[4]);
    
    per_master_nodes.push_back(corner[1]);
    per_slave_nodes.push_back(corner[0]);
    
    per_master_nodes.push_back(corner[7]);
    per_slave_nodes.push_back(corner[6]);
    
    per_master_nodes.push_back(corner[3]);
    per_slave_nodes.push_back(corner[2]);

    // Impose Dirichlet nodes on Right, e23, e67, e26, e37, c2,3,6,7
    // Note: There are nodes that are both per and dir. It is OK for now. It
    // will be handled in the creation of ID array. We will specify Dir nodes in
    // the last step in ID generation, so the repeated nodes will have Dirichlet
    // condition only.
    dir_nodes.insert( dir_nodes.end(), left_right + nlr, left_right + 2*nlr );
    dir_nodes.insert( dir_nodes.end(), edge01 + nex, edge01 + 2*nex );
    dir_nodes.insert( dir_nodes.end(), edge01 + 3*nex, edge01 + 4*nex );
    dir_nodes.insert( dir_nodes.end(), edge04 + 2*nez, edge04 + 4*nez );
    dir_nodes.push_back(corner[2]);
    dir_nodes.push_back(corner[3]);
    dir_nodes.push_back(corner[6]);
    dir_nodes.push_back(corner[7]);

    delete [] back_front; delete [] left_right; delete [] bottom_top;
    delete [] edge01; delete [] edge02; delete [] edge04; delete [] corner;
  }

  // Now we create a pesudoID array
  std::vector<int> pID;
  pID.clear();
  pID.resize(mesh->get_nFunc());
  for(unsigned int ii=0; ii<pID.size(); ++ii)
    pID[ii] = ii;
  for(unsigned int ii=0; ii<per_master_nodes.size(); ++ii)
    pID[per_slave_nodes[ii]] = per_master_nodes[ii];


  // Now we loop over each interior pair
  int patchm, patchs, facem, faces, patchm_start, patchs_start;
  int master_index, slave_index;
  double mx, my, mz, sx, sy, sz;
  double dist;
  for(unsigned int ii=0; ii<intpat1.size(); ++ii)
  {
    // we set the patch with smaller index as the slave; with bigger index as
    // the master. If two patch has identical index, report error.
    if( intpat1[ii] > intpat2[ii] )
    {
      patchm = intpat1[ii]; facem = intfac1[ii];
      patchs = intpat2[ii]; faces = intfac2[ii];
    }
    else if( intpat1[ii] < intpat2[ii] )
    {
      patchm = intpat2[ii]; facem = intfac2[ii];
      patchs = intpat1[ii]; faces = intfac1[ii];
    }
    else
    {
      std::cerr<<"Error: the "<<ii<<"th interior interface pair involves patches "<<
        intpat1[ii]<<" and "<<intpat2[ii]<<std::endl;
      exit(EXIT_FAILURE);
    }
   
    if(patchm >= numPat || patchs >= numPat)
    {
      std::cerr<<"Error: the patches have index "<<patchm<<'\t'<<patchs
        <<". This BC is for "<<numPat<<" patch geometry. \n";
      exit(EXIT_FAILURE);
    }

    // Now we determine the starting index for the top/bottom interface.
    int nFunc_z_m = mesh->get_patch_mesh(patchm)->get_nFunc_z();
    int nFunc_z_s = mesh->get_patch_mesh(patchs)->get_nFunc_z();
    int f_start_m = mesh->get_patch_mesh(patchm)->get_nFunc_start();
    int f_start_s = mesh->get_patch_mesh(patchs)->get_nFunc_start(); 

    if( facem == 1)
      patchm_start = f_start_m;
    else if( facem == 6 )
      patchm_start = f_start_m + (nFunc_z_m-1) * nFunc_x_0 * nFunc_y_0;
    else
    {
      std::cerr<<"Error: the "<<ii
        <<"th interior master interface has non-top/non-bottom face involved. \n";
      exit(EXIT_FAILURE);
    } 

    if( faces == 1 )
      patchs_start = f_start_s;
    else if( faces == 6 )
      patchs_start = f_start_s + (nFunc_z_s-1) * nFunc_x_0 * nFunc_y_0;
    else
    {
      std::cerr<<"Error: the "<<ii
        <<"th interior slave interface has non-top/non-bottom face involved. \n";
      exit(EXIT_FAILURE);
    } 

    int counter = 0;
    for(int jj=0; jj<nFunc_y_0; ++jj)
    {
      for(int kk=0; kk<nFunc_x_0; ++kk)
      {
        master_index = patchm_start + jj * nFunc_x_0 + kk;
        mx = cPts[master_index * 4 + 0];
        my = cPts[master_index * 4 + 1];
        mz = cPts[master_index * 4 + 2];
        
        for(int nn=0; nn<nFunc_y_0; ++nn)
        {
          for(int mm=0; mm<nFunc_x_0; ++mm)
          {
            slave_index = patchs_start + nn * nFunc_x_0 + mm;
            sx = cPts[slave_index * 4 + 0];
            sy = cPts[slave_index * 4 + 1];
            sz = cPts[slave_index * 4 + 2];
            
            dist = sqrt( (sx - mx) * (sx - mx) + (sy - my) * (sy - my)
               + (sz - mz) * (sz - mz) ); 
            
            if( dist < tol )
            {
              counter += 1;
              
              if( pID[slave_index] < master_index )
                pID[slave_index] = master_index;
            }
          }
        }
      }
    }

  } // Finish the loop for interior interfaces

  for(unsigned int ii=0; ii<dir_nodes.size(); ++ii)
    pID[dir_nodes[ii]] = -1;


  per_master_nodes.clear();
  per_slave_nodes.clear();
  for(int ii=0; ii<(int) pID.size(); ++ii)
  {
    if(pID[ii] > 0)
    {
      if(pID[ii] > ii)
      {
        per_master_nodes.push_back(pID[ii]);
        per_slave_nodes.push_back(ii);
      }
      else if(pID[ii] < ii)
      {
        std::cerr<<"Error: pseudoID cannot be smaller that its index. \n";
        exit(EXIT_FAILURE);
      }
    }
  }

  num_per_nodes = per_master_nodes.size();
  num_dir_nodes = dir_nodes.size();


  // Now element face integral specification for NBC
  front_elem.clear(); back_elem.clear(); left_elem.clear(); right_elem.clear();
  top_elem.clear(); bottom_elem.clear();

  num_front_elem = 0; num_back_elem = 0;
  num_left_elem = 0; num_right_elem = 0;

  VEC_T::shrink2fit(front_elem); VEC_T::shrink2fit(back_elem);
  VEC_T::shrink2fit(left_elem); VEC_T::shrink2fit(right_elem);

  std::vector<int> front, back, left, right, top, bottom;
  Generate_BCElems_B(mesh->get_patch_mesh(0), front, back, left, right, top, bottom);

  bottom_elem.insert(bottom_elem.end(), bottom.begin(), bottom.end());

  VEC_T::shrink2fit(bottom_elem); VEC_T::shrink2fit(top_elem);

  num_top_elem = top_elem.size();
  num_bottom_elem = bottom_elem.size();

  std::cout<<"-----> 107-Patch Coronary geometry: NBC on patch 0  bottom; ";
  std::cout<<" Front-Back C0 periodic; Right Dirichlet, Left nothing. \n";


}


void BoundaryCond::BC_type_201( const IMesh * const &mesh, const std::vector<double> &cPts )
{
  const int numPat = mesh->get_num_patch();
  
  const double tol = 1.0e-6;

  if(numPat != 107)
  {
    std::cerr<<"Error: BC_type_201 is designed for 107-patch geometry only. \n";
    exit(EXIT_FAILURE);
  }

  char * char_home_dir = getenv("HOME");
  std::string interface_file(char_home_dir);
  interface_file.append("/IsoPETSc3D/trunk/input/pig_coronary_107patch/global_data.dat");

  // Two patches intpat1[i] and intpat2[i] meet with intpat1's intfac1[i] and
  // intpat2's intfac2[i]. 
  std::vector<int> intpat1, intpat2, intfac1, intfac2;

  // Read in the values for the four int std::vectors from data file to specify all
  // interior interfaces.
  // Note: in shaolie's numbering, first patch starts with index 1. In the
  // reader, we correct it to start with 0.
  Read_SH_interface_pair( interface_file, intpat1, intpat2, intfac1, intfac2 );
  
  std::cout<<"Read interface interface pair from "<<interface_file<<std::endl;


  const int nFunc_x_0 = mesh->get_patch_mesh(0)->get_nFunc_x();
  const int nFunc_y_0 = mesh->get_patch_mesh(0)->get_nFunc_y();

  
  per_master_nodes.clear();
  per_slave_nodes.clear();
  
  dir_nodes.clear();

  // Now we loop over patches
  for(int ii=0; ii<numPat; ++ii)
  {
    int nFunc_x = mesh->get_patch_mesh(ii)->get_nFunc_x();
    int nFunc_y = mesh->get_patch_mesh(ii)->get_nFunc_y();
    int nFunc_z = mesh->get_patch_mesh(ii)->get_nFunc_z();
    int f_start = mesh->get_patch_mesh(ii)->get_nFunc_start();

    if( nFunc_x != nFunc_x_0 || nFunc_y != nFunc_y_0 )
    {
      std::cerr<<"Error: number of basis functions in x/y-direction does not match for patch "
        <<ii<<" !"<<std::endl;
      exit(EXIT_FAILURE);
    }

    int * back_front, * left_right, * bottom_top, * edge01, * edge02, * edge04, * corner;
    int nbf, nlr, nbt, nex, ney, nez;

    // Generate the indices on all boundary faces, edges, points.
    Generate_BCNodes_D(nFunc_x, nFunc_y, nFunc_z, f_start,
        back_front, nbf, left_right, nlr, bottom_top, nbt,
        edge01, nex, edge02, ney, edge04, nez, corner);

    // Impose periodic pairs for front - back, 
    // e15 - e04, e13 - e02, e37 - e26, e57 - e46
    // c1 - c0, c5 - c4, c3 - c2, c7 - c6
    per_master_nodes.insert( per_master_nodes.end(), back_front+nbf, back_front+2*nbf );
    per_slave_nodes.insert( per_slave_nodes.end(), back_front, back_front + nbf );

    per_master_nodes.insert( per_master_nodes.end(), edge04 + nez, edge04 + 2*nez );
    per_slave_nodes.insert( per_slave_nodes.end(), edge04, edge04 + nez );
    
    per_master_nodes.insert( per_master_nodes.end(), edge02+ney, edge02 + 2*ney );
    per_slave_nodes.insert( per_slave_nodes.end(), edge02, edge02 + ney );

    per_master_nodes.insert( per_master_nodes.end(), edge04+3*nez, edge04+4*nez );
    per_slave_nodes.insert( per_slave_nodes.end(), edge04+2*nez, edge04+3*nez );

    per_master_nodes.insert( per_master_nodes.end(), edge02+3*ney, edge02+4*ney );
    per_slave_nodes.insert( per_slave_nodes.end(), edge02+2*ney, edge02+3*ney);

    per_master_nodes.push_back(corner[5]);
    per_slave_nodes.push_back(corner[4]);
    
    per_master_nodes.push_back(corner[1]);
    per_slave_nodes.push_back(corner[0]);
    
    per_master_nodes.push_back(corner[7]);
    per_slave_nodes.push_back(corner[6]);
    
    per_master_nodes.push_back(corner[3]);
    per_slave_nodes.push_back(corner[2]);


    // Master-slave relation on LEFT face, e01, e45
    for( int zz = 0; zz < nFunc_z; ++zz )
    {
      int mnode = zz * nFunc_x * nFunc_y + nFunc_x - 1;
      for( int xx = 1; xx < nFunc_x - 1; ++xx )
      {
        int snode = zz * nFunc_x * nFunc_y + xx;

        per_master_nodes.push_back(mnode);
        per_slave_nodes.push_back(snode);
      }
    }


    // Impose Dirichlet nodes on Right, e23, e67, e26, e37, c2,3,6,7
    // Note: There are nodes that are both per and dir. It is OK for now. It
    // will be handled in the creation of ID array. We will specify Dir nodes in
    // the last step in ID generation, so the repeated nodes will have Dirichlet
    // condition only.
    dir_nodes.insert( dir_nodes.end(), left_right + nlr, left_right + 2*nlr );
    dir_nodes.insert( dir_nodes.end(), edge01 + nex, edge01 + 2*nex );
    dir_nodes.insert( dir_nodes.end(), edge01 + 3*nex, edge01 + 4*nex );
    dir_nodes.insert( dir_nodes.end(), edge04 + 2*nez, edge04 + 4*nez );
    dir_nodes.push_back(corner[2]);
    dir_nodes.push_back(corner[3]);
    dir_nodes.push_back(corner[6]);
    dir_nodes.push_back(corner[7]);

    delete [] back_front; delete [] left_right; delete [] bottom_top;
    delete [] edge01; delete [] edge02; delete [] edge04; delete [] corner;
  }

  // Now we create a pesudoID array
  std::vector<int> pID;
  pID.clear();
  pID.resize(mesh->get_nFunc());
  for(unsigned int ii=0; ii<pID.size(); ++ii)
    pID[ii] = ii;
  for(unsigned int ii=0; ii<per_master_nodes.size(); ++ii)
    pID[per_slave_nodes[ii]] = per_master_nodes[ii];


  // Now we loop over each interior pair
  int patchm, patchs, facem, faces, patchm_start, patchs_start;
  int master_index, slave_index;
  double mx, my, mz, sx, sy, sz;
  double dist;
  for(unsigned int ii=0; ii<intpat1.size(); ++ii)
  {
    // we set the patch with smaller index as the slave; with bigger index as
    // the master. If two patch has identical index, report error.
    if( intpat1[ii] > intpat2[ii] )
    {
      patchm = intpat1[ii]; facem = intfac1[ii];
      patchs = intpat2[ii]; faces = intfac2[ii];
    }
    else if( intpat1[ii] < intpat2[ii] )
    {
      patchm = intpat2[ii]; facem = intfac2[ii];
      patchs = intpat1[ii]; faces = intfac1[ii];
    }
    else
    {
      std::cerr<<"Error: the "<<ii<<"th interior interface pair involves patches "<<
        intpat1[ii]<<" and "<<intpat2[ii]<<std::endl;
      exit(EXIT_FAILURE);
    }

    if(patchm >= numPat || patchs >= numPat)
    {
      std::cerr<<"Error: the patches have index "<<patchm<<'\t'<<patchs
        <<". This BC is for "<<numPat<<" patch geometry. \n";
      exit(EXIT_FAILURE);
    }

    // Now we determine the starting index for the top/bottom interface.
    int nFunc_z_m = mesh->get_patch_mesh(patchm)->get_nFunc_z();
    int nFunc_z_s = mesh->get_patch_mesh(patchs)->get_nFunc_z();
    int f_start_m = mesh->get_patch_mesh(patchm)->get_nFunc_start();
    int f_start_s = mesh->get_patch_mesh(patchs)->get_nFunc_start(); 

    if( facem == 1)
      patchm_start = f_start_m;
    else if( facem == 6 )
      patchm_start = f_start_m + (nFunc_z_m-1) * nFunc_x_0 * nFunc_y_0;
    else
    {
      std::cerr<<"Error: the "<<ii
        <<"th interior master interface has non-top/non-bottom face involved. \n";
      exit(EXIT_FAILURE);
    } 

    if( faces == 1 )
      patchs_start = f_start_s;
    else if( faces == 6 )
      patchs_start = f_start_s + (nFunc_z_s-1) * nFunc_x_0 * nFunc_y_0;
    else
    {
      std::cerr<<"Error: the "<<ii
        <<"th interior slave interface has non-top/non-bottom face involved. \n";
      exit(EXIT_FAILURE);
    } 

    int counter = 0;
    for(int jj=0; jj<nFunc_y_0; ++jj)
    {
      for(int kk=0; kk<nFunc_x_0; ++kk)
      {
        master_index = patchm_start + jj * nFunc_x_0 + kk;
        mx = cPts[master_index * 4 + 0];
        my = cPts[master_index * 4 + 1];
        mz = cPts[master_index * 4 + 2];

        for(int nn=0; nn<nFunc_y_0; ++nn)
        {
          for(int mm=0; mm<nFunc_x_0; ++mm)
          {
            slave_index = patchs_start + nn * nFunc_x_0 + mm;
            sx = cPts[slave_index * 4 + 0];
            sy = cPts[slave_index * 4 + 1];
            sz = cPts[slave_index * 4 + 2];

            dist = sqrt( (sx - mx) * (sx - mx) + (sy - my) * (sy - my)
                + (sz - mz) * (sz - mz) ); 

            if( dist < tol )
            {
              counter += 1;

              if( pID[slave_index] < master_index )
                pID[slave_index] = master_index;
            }
          }
        }
      }
    }

  } // Finish the loop for interior interfaces

  for(unsigned int ii=0; ii<dir_nodes.size(); ++ii)
    pID[dir_nodes[ii]] = -1;


  per_master_nodes.clear();
  per_slave_nodes.clear();
  for(int ii=0; ii<(int) pID.size(); ++ii)
  {
    if(pID[ii] > 0)
    {
      if(pID[ii] > ii)
      {
        per_master_nodes.push_back(pID[ii]);
        per_slave_nodes.push_back(ii);
      }
      else if(pID[ii] < ii)
      {
        std::cerr<<"Error: pseudoID cannot be smaller that its index. \n";
        exit(EXIT_FAILURE);
      }
    }
  }

  num_per_nodes = per_master_nodes.size();
  num_dir_nodes = dir_nodes.size();


  // Now element face integral specification for NBC
  front_elem.clear(); back_elem.clear(); left_elem.clear(); right_elem.clear();
  top_elem.clear(); bottom_elem.clear();

  num_front_elem = 0; num_back_elem = 0;
  num_left_elem = 0; num_right_elem = 0;

  VEC_T::shrink2fit(front_elem); VEC_T::shrink2fit(back_elem);
  VEC_T::shrink2fit(left_elem); VEC_T::shrink2fit(right_elem);

  std::vector<int> front, back, left, right, top, bottom;
  Generate_BCElems_B(mesh->get_patch_mesh(0), front, back, left, right, top, bottom);

  bottom_elem.insert(bottom_elem.end(), bottom.begin(), bottom.end());

  VEC_T::shrink2fit(bottom_elem); VEC_T::shrink2fit(top_elem);

  num_top_elem = top_elem.size();
  num_bottom_elem = bottom_elem.size();

  std::cout<<"-----> 107-Patch Coronary geometry: NBC on patch 0  bottom; ";
  std::cout<<" Front-Back C0 periodic; Right Dirichlet, Left master-slave pairing with e15. \n";

}

// EOF
