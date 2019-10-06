#include "IElemBC.hpp"

IElemBC::IElemBC()
{
  left_elem.clear();
  right_elem.clear();
  top_elem.clear();
  bottom_elem.clear();
  front_elem.clear();
  back_elem.clear();
}


IElemBC::~IElemBC()
{
  VEC_T::clean(left_elem);
  VEC_T::clean(right_elem);
  VEC_T::clean(top_elem);
  VEC_T::clean(bottom_elem);
  VEC_T::clean(front_elem);
  VEC_T::clean(back_elem);
}



void IElemBC::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"======== Elem BC info ========"<<std::endl;
  std::cout<<"bottom elem: \t";
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
  std::cout<<std::endl<<"front elem: \t";
  for(unsigned int ii=0; ii<get_num_front_elem(); ++ii)
    std::cout<<get_front_elem(ii)<<'\t';
  std::cout<<std::endl<<"back elem: \t";
  for(unsigned int ii=0; ii<get_num_back_elem(); ++ii)
    std::cout<<get_back_elem(ii)<<'\t';
  std::cout<<std::endl<<"=============================="<<std::endl;
}



void IElemBC::Generate_BCElem_2D( const int &nElem_x, const int &nElem_y,
    std::vector<unsigned int> &fro, std::vector<unsigned int> &bac,
    std::vector<unsigned int> &lef, std::vector<unsigned int> &rig )
{
  lef.clear(); rig.clear(); fro.clear(); bac.clear();

  for(int ii=0; ii<nElem_y; ++ii)
  {
    bac.push_back( ii * nElem_x );
    fro.push_back( (ii + 1)*nElem_x - 1  );
  }

  for(int ii=0; ii<nElem_x; ++ii)
  {
    lef.push_back(ii);
    rig.push_back(ii + nElem_x * (nElem_y-1));
  }
}



void IElemBC::Generate_BCElem_3D_A( const int &nElem_x, const int &nElem_y,
    const int &nElem_z,
    std::vector<unsigned int> &front, std::vector<unsigned int> &back, 
    std::vector<unsigned int> &left, std::vector<unsigned int> &right, 
    std::vector<unsigned int> &top, std::vector<unsigned int> &bottom ) const
{
  // Reads single patch NURBS information from mesh, then determine the 
  // elements that have faces on the front, back, ..., bottom boundary faces.
  front.clear(); back.clear(); left.clear(); right.clear();
  top.clear(); bottom.clear();

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



void IElemBC::Generate_BCElem_3D_B( const int &nElem_x,
    const int &nElem_y, const int &nElem_z, const int &es,
    std::vector<unsigned int> &front, std::vector<unsigned int> &back, 
    std::vector<unsigned int> &left, std::vector<unsigned int> &right, 
    std::vector<unsigned int> &top, std::vector<unsigned int> &bottom ) const
{
  front.clear(); back.clear(); left.clear(); right.clear();
  top.clear(); bottom.clear();

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


// EOF
