#include "INodalBC.hpp"

INodalBC::INodalBC()
{}


INodalBC::~INodalBC()
{
  VEC_T::clean(dir_nodes);
  VEC_T::clean(per_slave_nodes);
  VEC_T::clean(per_master_nodes);
  VEC_T::clean(ID);

  VEC_T::clean(num_dir_nodes);
  VEC_T::clean(num_per_nodes);
}


void INodalBC::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"======== BC info ======="<<std::endl;
  if(num_dir_nodes.size() > 0)
  {
    std::cout<<"Dirichlet nodes: "<<std::endl;
    for(unsigned int ii=0; ii<num_dir_nodes.size(); ++ii)
    {
      std::cout<<"nbc_id " << ii << ": ";
      for(unsigned int jj=0; jj<num_dir_nodes[ii]; ++jj)
        std::cout<<dir_nodes[ii][jj]<<'\t';
      std::cout<<std::endl;
    }
  }

  if(num_per_nodes.size() > 0)
  {
    std::cout<<"Periodic master - slave nodes: "<<std::endl;
    for(unsigned int ii=0; ii<num_per_nodes.size(); ++ii)
    {
      std::cout<<"nbc_id " << ii << ": ";
      for(unsigned int jj=0; jj<num_per_nodes[ii]; ++jj)
        std::cout<<per_master_nodes[ii][jj]<<'\t'<<per_slave_nodes[ii][jj]<<std::endl;
    }
  }

  std::cout<<std::endl<<"ID array: "<<std::endl;
  for(unsigned int ii=0; ii<ID.size(); ++ii)
    std::cout<<ID[ii]<<'\t';
  std::cout<<'\n';
  std::cout<<std::endl<<"========================"<<std::endl;
}


void INodalBC::Create_ID(const unsigned int &num_node)
{
  ID.clear();
  ID.resize(num_node);
  VEC_T::shrink2fit(ID);

  for(unsigned int ii = 0; ii<ID.size(); ++ii)
    ID[ii] = ii;

  for(unsigned int ii = 0; ii<dir_nodes.size(); ++ii)
  {
    for(unsigned int jj = 0; jj<dir_nodes[ii].size(); ++jj)
      ID[ dir_nodes[ii][jj] ] = -1;
  }

  for(unsigned int ii = 0; ii<per_slave_nodes.size(); ++ii)
  {
    for(unsigned int jj = 0; jj<per_slave_nodes[ii].size(); ++jj)
      ID[ per_slave_nodes[ii][jj] ] = per_master_nodes[ii][jj];
  }
}


void INodalBC::Generate_BCNode_2D_A( const int &nFunc_x, const int &nFunc_y,
    std::vector<unsigned int> &left, std::vector<unsigned int> &right,
    std::vector<unsigned int> &top, std::vector<unsigned int> &bottom,
    std::vector<unsigned int> &corner ) const
{
  left.clear(); right.clear(); top.clear(); bottom.clear(); corner.clear();

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

  corner.push_back(0); 
  corner.push_back(nFunc_x - 1);
  corner.push_back(nFunc_x * (nFunc_y-1));
  corner.push_back(nFunc_x * nFunc_y - 1);

  VEC_T::shrink2fit(left); VEC_T::shrink2fit(right);
  VEC_T::shrink2fit(top);  VEC_T::shrink2fit(bottom);
  VEC_T::shrink2fit(corner);
}


void INodalBC::Generate_BCNode_2D_B( const int &nFunc_x, const int &nFunc_y,
    std::vector<unsigned int> &lef1, std::vector<unsigned int> &rig1,
    std::vector<unsigned int> &top1, std::vector<unsigned int> &bot1,
    std::vector<unsigned int> &lef2, std::vector<unsigned int> &rig2,
    std::vector<unsigned int> &top2, std::vector<unsigned int> &bot2,
    std::vector<unsigned int> &corner1, std::vector<unsigned int> &corner2,
    std::vector<unsigned int> &corner3, std::vector<unsigned int> &corner4 ) const
{
  lef1.clear(); rig1.clear(); top1.clear(); bot1.clear();
  lef2.clear(); rig2.clear(); top2.clear(); bot2.clear();
  corner1.clear(); corner2.clear(); corner3.clear(); corner4.clear();

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
  corner3.push_back(2*nFunc_x - 1);
  corner3.push_back(nFunc_x * (nFunc_y-2));
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


void INodalBC::Generate_BCNode_3D_A( const int &nFunc_x,
    const int &nFunc_y, const int &nFunc_z,
    std::vector<int> &front, std::vector<int> &back, std::vector<int> &left,
    std::vector<int> &right, std::vector<int> &top, std::vector<int> &bottom ) const
{
  // This function reads in basic mesh object, based on the nfunc,
  // nFunc_x, ... , nElem, etc. to determine the boundary nodes/elements.
  // This approach is restricted to NURBS single patch mesh.

  front.clear(); back.clear(); left.clear(); right.clear();
  top.clear(); bottom.clear();

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


void INodalBC::Generate_BCNode_3D_B( const int &nFunc_x, const int &nFunc_y,
    const int &nFunc_z,
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


void INodalBC::Generate_BCNode_3D_C( const int &nFunc_x,
    const int &nFunc_y, const int &nFunc_z,
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


void INodalBC::Generate_BCNode_3D_D( const int &nFunc_x, const int &nFunc_y,
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
    std::cerr<<"Error: in Generate_BCNode_3D_D, nFunc_z has to be greater than 2. \n";
    exit(EXIT_FAILURE);
  }
  if(nFunc_y <= 2)
  {
    std::cerr<<"Error: in Generate_BCNode_3D_D, nFunc_y has to be greater than 2. \n";
    exit(EXIT_FAILURE);
  }
  if(nFunc_x <= 2)
  {
    std::cerr<<"Error: in Generate_BCNode_3D_D, nFunc_x has to be greater than 2. \n";
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

// EOF
