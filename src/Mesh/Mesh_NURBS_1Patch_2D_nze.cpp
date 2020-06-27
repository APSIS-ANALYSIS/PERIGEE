#include "Mesh_NURBS_1Patch_2D_nze.hpp"

Mesh_NURBS_1Patch_2D_nze::Mesh_NURBS_1Patch_2D_nze( IMesh const * const &inmesh )
{
  s_degree = inmesh->get_s_degree();
  t_degree = inmesh->get_t_degree();
  
  hx_max   = inmesh->get_hx_max();
  hy_max   = inmesh->get_hy_max();
  hx_min   = inmesh->get_hx_min();
  hy_min   = inmesh->get_hy_min();

  nFunc_x  = inmesh->get_nFunc_x();
  nFunc_y  = inmesh->get_nFunc_y();
  nFunc    = inmesh->get_nFunc();

  nElem_x  = inmesh->get_nElem_x_nz();
  nElem_y  = inmesh->get_nElem_y_nz();
  nElem    = inmesh->get_nElem_nz();

  nLocBas  = inmesh->get_nLocBas();

  hx.clear();
  hy.clear();

  double temp = 0;
  const int nex = inmesh->get_nElem_x();
  for(int ii=0; ii<nex; ++ii)
  {
    temp = inmesh->get_hx(ii);
    if(temp > 0.0)
      hx.push_back(temp);
  }

  for(int ii=0; ii<inmesh->get_nElem_y(); ++ii)
  {
    temp = inmesh->get_hy(ii*nex);
    if(temp > 0.0)
      hy.push_back(temp);
  }

  VEC_T::shrink2fit(hx);
  VEC_T::shrink2fit(hy);

  if(int(hx.size()) != nElem_x)
  {
    std::cerr<<"Error: hx.size = "<<hx.size()<<" nElem_x = "<<nElem_x<<std::endl;
    exit(EXIT_FAILURE); 
  }

  if(int(hy.size()) != nElem_y)
  {
    std::cerr<<"Error: hy.size = "<<hy.size()<<" nElem_y = "<<nElem_y<<std::endl;
    exit(EXIT_FAILURE); 
  }
}


Mesh_NURBS_1Patch_2D_nze::~Mesh_NURBS_1Patch_2D_nze()
{
  std::cout<<"-- Mesh_NURBS_1Patch_2D_nze deleted. \n";
}


void Mesh_NURBS_1Patch_2D_nze::get_elem_index( const int &ee,
    int &ex, int &ey ) const
{
  ex = ee % nElem_x;
  ey = (ee - ex) / nElem_x;
}


double Mesh_NURBS_1Patch_2D_nze::get_hx( int ee ) const
{
  int ex, ey;
  get_elem_index(ee, ex, ey);
  return hx[ex];
}


double Mesh_NURBS_1Patch_2D_nze::get_hy( int ee ) const
{
  int ex, ey;
  get_elem_index(ee, ex, ey);
  return hy[ey];
}


double Mesh_NURBS_1Patch_2D_nze::get_hz( int ee ) const
{
  std::cerr<<"Error: this is a 2D mesh class, get_hz is invalid! \n";

  return -1.0;
}


void Mesh_NURBS_1Patch_2D_nze::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"====== Mesh_NURBS_1Patch_2D ======"<<std::endl;
  std::cout<<"Degree: "<<get_s_degree()<<'\t'<<get_t_degree()<<std::endl;
  std::cout<<"Max mesh size: "<<get_hx_max()<<'\t'<<get_hy_max()<<std::endl;
  std::cout<<"Min mesh size: "<<get_hx_min()<<'\t'<<get_hy_min()<<std::endl;
  std::cout<<"Elem #: "<<get_nElem_x()<<'\t'<<get_nElem_y()<<std::endl;
  std::cout<<"Total Elem: "<<get_nElem()<<std::endl;
  std::cout<<"Func #: "<<get_nFunc_x()<<'\t'<<get_nFunc_y()<<std::endl;
  std::cout<<"Total Func: "<<get_nFunc()<<std::endl;
  std::cout<<"Local Basis #: "<<get_nLocBas()<<std::endl;
  std::cout<<"Patch index: "<<get_patch_index()<<std::endl;
  std::cout<<"Start element: "<<get_nElem_start()<<std::endl;
  std::cout<<"Start function: "<<get_nFunc_start()<<std::endl;
  std::cout<<"=================================="<<std::endl;
}

// EOF
