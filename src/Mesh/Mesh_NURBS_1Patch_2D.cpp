#include "Mesh_NURBS_1Patch_2D.hpp"

Mesh_NURBS_1Patch_2D::Mesh_NURBS_1Patch_2D( int in_s_degree, int in_t_degree,
    double input_hx_max, double input_hy_max,
    double input_hx_min, double input_hy_min,
    const std::vector<double> &in_sknot,
    const std::vector<double> &in_tknot )
{
  s_degree = in_s_degree;
  t_degree = in_t_degree;

  hx_max = input_hx_max;
  hy_max = input_hy_max;
  hx_min = input_hx_min;
  hy_min = input_hy_min;

  nFunc_x = in_sknot.size() - s_degree - 1;
  nFunc_y = in_tknot.size() - t_degree - 1;
  nFunc   = nFunc_x * nFunc_y;

  nElem_x = nFunc_x - s_degree;
  nElem_y = nFunc_y - t_degree;
  nElem   = nElem_x * nElem_y;

  nLocBas = ( s_degree + 1 ) * ( t_degree + 1 );

  hx.clear(); hy.clear();

  for(int ii=0; ii<nElem_x; ++ii)
    hx.push_back(in_sknot[s_degree+ii+1] - in_sknot[s_degree+ii]);

  for(int ii=0; ii<nElem_y; ++ii)
    hy.push_back(in_tknot[t_degree+ii+1] - in_tknot[t_degree+ii]);

  VEC_T::shrink2fit(hx);
  VEC_T::shrink2fit(hy);

  int num_zero_elem_x = 0;
  for(int ii=0; ii<nElem_x; ++ii)
  {
    if(hx[ii] == 0.0)
      num_zero_elem_x += 1;
  }

  int num_zero_elem_y = 0;
  for(int ii=0; ii<nElem_y; ++ii)
  {
    if(hy[ii] == 0.0)
      num_zero_elem_y += 1;
  }

  // NURBS has tensor-product structure. Hence, the nonzero element has the same
  // structure and we can multiply the number of nonzero element in each
  // direction to get the total number of nonzero element
  nElem_x_nz = nElem_x - num_zero_elem_x;
  nElem_y_nz = nElem_y - num_zero_elem_y;
  nElem_nz = nElem_x_nz * nElem_y_nz;
}


Mesh_NURBS_1Patch_2D::~Mesh_NURBS_1Patch_2D()
{
  std::cout<<"-- Mesh_NURBS_1Patch_2D deleted. \n";
}


void Mesh_NURBS_1Patch_2D::get_elem_index( const int &ee,
    int &ex, int &ey ) const
{
  ex = ee % nElem_x;
  ey = (ee - ex) / nElem_x;
}


double Mesh_NURBS_1Patch_2D::get_hx( int ee ) const
{
  int ex, ey;

  get_elem_index(ee, ex, ey);

  return hx[ex];
}


double Mesh_NURBS_1Patch_2D::get_hy( int ee ) const 
{
  int ex, ey;

  get_elem_index(ee, ex, ey);

  return hy[ey];
}


double Mesh_NURBS_1Patch_2D::get_hz( int ee ) const
{
  std::cerr<<"Error: this is a 2D mesh class, get_hz is invalid! \n";

  return -1.0;
}


void Mesh_NURBS_1Patch_2D::print_mesh_info() const
{
  std::cout<<std::endl;
  std::cout<<"====== Mesh_NURBS_1Patch_2D ======"<<std::endl;
  std::cout<<"Degree: "<<get_s_degree()<<'\t'<<get_t_degree()<<std::endl;
  std::cout<<"Max mesh size: "<<get_hx_max()<<'\t'<<get_hy_max()<<std::endl;
  std::cout<<"Min mesh size: "<<get_hx_min()<<'\t'<<get_hy_min()<<std::endl;
  std::cout<<"Elem #: "<<get_nElem_x()<<'\t'<<get_nElem_y()<<std::endl;
  std::cout<<"Elem_nz #: "<<get_nElem_x_nz()<<'\t'<<get_nElem_y_nz()<<std::endl;
  std::cout<<"Total Elem: "<<get_nElem()<<std::endl;
  std::cout<<"Total nonzero Elem: "<<get_nElem_nz()<<std::endl;
  std::cout<<"Func #: "<<get_nFunc_x()<<'\t'<<get_nFunc_y()<<std::endl;
  std::cout<<"Total Func: "<<get_nFunc()<<std::endl;
  std::cout<<"Local Basis #: "<<get_nLocBas()<<std::endl;
  std::cout<<"Patch index: "<<get_patch_index()<<std::endl;
  std::cout<<"Start element: "<<get_nElem_start()<<std::endl;
  std::cout<<"Start function: "<<get_nFunc_start()<<std::endl;
  std::cout<<"=================================="<<std::endl;
}


// EOF
