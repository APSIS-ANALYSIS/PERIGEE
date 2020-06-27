#include "Mesh_NURBS_1Patch_3D.hpp"

Mesh_NURBS_1Patch_3D::Mesh_NURBS_1Patch_3D( 
    int in_s_degree, int in_t_degree, int in_u_degree,
    double input_hx_max, double input_hy_max, double input_hz_max,
    double input_hx_min, double input_hy_min, double input_hz_min,
    const std::vector<double> &in_sknot, const std::vector<double> &in_tknot,
    const std::vector<double> &in_uknot 
    )
{
  s_degree = in_s_degree;
  t_degree = in_t_degree;
  u_degree = in_u_degree;

  hx_max = input_hx_max;
  hy_max = input_hy_max;
  hz_max = input_hz_max;

  hx_min = input_hx_min;
  hy_min = input_hy_min;
  hz_min = input_hz_min;

  nFunc_x = in_sknot.size() - s_degree - 1;
  nFunc_y = in_tknot.size() - t_degree - 1;
  nFunc_z = in_uknot.size() - u_degree - 1;
  nFunc   = nFunc_x * nFunc_y * nFunc_z;

  nElem_x = nFunc_x - s_degree;
  nElem_y = nFunc_y - t_degree;
  nElem_z = nFunc_z - u_degree;
  nElem   = nElem_x * nElem_y * nElem_z;

  nLocBas = (s_degree + 1) * (t_degree + 1) * (u_degree + 1);

  hx.clear(); hy.clear(); hz.clear();

  for(int ii=0; ii<nElem_x; ++ii)
    hx.push_back(in_sknot[s_degree+ii+1] - in_sknot[s_degree+ii]);

  for(int ii=0; ii<nElem_y; ++ii)
    hy.push_back(in_tknot[t_degree+ii+1] - in_tknot[t_degree+ii]);

  for(int ii=0; ii<nElem_z; ++ii)
    hz.push_back(in_uknot[u_degree+ii+1] - in_uknot[u_degree+ii]);

  VEC_T::shrink2fit(hx);
  VEC_T::shrink2fit(hy);
  VEC_T::shrink2fit(hz);

  int num_zero_elem_x = 0;
  int num_zero_elem_y = 0;
  int num_zero_elem_z = 0;
  for(int ii=0; ii<nElem_x; ++ii)
  {
    if(hx[ii] == 0.0) num_zero_elem_x += 1;
  }
  for(int ii=0; ii<nElem_y; ++ii)
  {
    if(hy[ii] == 0.0) num_zero_elem_y += 1;
  }
  for(int ii=0; ii<nElem_z; ++ii)
  {
    if(hz[ii] == 0.0) num_zero_elem_z += 1;
  }

  nElem_x_nz = nElem_x - num_zero_elem_x;
  nElem_y_nz = nElem_y - num_zero_elem_y;
  nElem_z_nz = nElem_z - num_zero_elem_z;
  nElem_nz = nElem_x_nz * nElem_y_nz * nElem_z_nz;
}


Mesh_NURBS_1Patch_3D::Mesh_NURBS_1Patch_3D( 
    const int &in_s_degree, const int &in_t_degree, const int &in_u_degree,
    const double &input_hx_max, const double &input_hy_max, const double &input_hz_max,
    const double &input_hx_min, const double &input_hy_min, const double &input_hz_min,
    const int &inelemx, const int &inelemy, const int &inelemz )
{
  s_degree = in_s_degree;
  t_degree = in_t_degree;
  u_degree = in_u_degree;

  hx_max = input_hx_max;
  hy_max = input_hy_max;
  hz_max = input_hz_max;

  hx_min = input_hx_min;
  hy_min = input_hy_min;
  hz_min = input_hz_min;

  nElem_x = inelemx;
  nElem_y = inelemy;
  nElem_z = inelemz;
  nElem   = nElem_x * nElem_y * nElem_z;

  nFunc_x = nElem_x + s_degree;
  nFunc_y = nElem_y + t_degree;
  nFunc_z = nElem_z + u_degree;
  nFunc   = nFunc_x * nFunc_y * nFunc_z;

  nLocBas = (s_degree + 1) * (t_degree +1) * (u_degree+1);

  hx.clear(); hy.clear(); hz.clear();
  const double hhx = 1.0 / double(nElem_x);
  const double hhy = 1.0 / double(nElem_y);
  const double hhz = 1.0 / double(nElem_z);
  for(int ii=0; ii<nElem_x; ++ii)
    hx.push_back(hhx);
  for(int ii=0; ii<nElem_y; ++ii)
    hy.push_back(hhy);
  for(int ii=0; ii<nElem_z; ++ii)
    hz.push_back(hhz);

  VEC_T::shrink2fit(hx);
  VEC_T::shrink2fit(hy);
  VEC_T::shrink2fit(hz);


  int num_zero_elem_x = 0;
  int num_zero_elem_y = 0;
  int num_zero_elem_z = 0;
  for(int ii=0; ii<nElem_x; ++ii)
  {
    if(hx[ii] == 0.0) num_zero_elem_x += 1;
  }
  for(int ii=0; ii<nElem_y; ++ii)
  {
    if(hy[ii] == 0.0) num_zero_elem_y += 1;
  }
  for(int ii=0; ii<nElem_z; ++ii)
  {
    if(hz[ii] == 0.0) num_zero_elem_z += 1;
  }

  nElem_x_nz = nElem_x - num_zero_elem_x;
  nElem_y_nz = nElem_y - num_zero_elem_y;
  nElem_z_nz = nElem_z - num_zero_elem_z;
  nElem_nz = nElem_x_nz * nElem_y_nz * nElem_z_nz;


  std::cout<<"--- Mesh_NURBS_1Patch_3D generated.";
  std::cout<<"This should be used for test code ONLY. \n";
}



Mesh_NURBS_1Patch_3D::~Mesh_NURBS_1Patch_3D()
{
  VEC_T::clean(hx);
  VEC_T::clean(hy);
  VEC_T::clean(hz);
  std::cout<<"-- Mesh_NURBS_1Patch_3D deleted. \n";
}

void Mesh_NURBS_1Patch_3D::get_elem_index( const int &ee, 
    int &ex, int &ey, int &ez ) const 
{
  int exy = ee % (nElem_x * nElem_y);
  ez = (ee - exy) / (nElem_x * nElem_y);
  ex = exy % nElem_x;
  ey = (exy - ex)/nElem_x;
}

double Mesh_NURBS_1Patch_3D::get_hx( int ee ) const
{
  int ex, ey, ez;
  get_elem_index(ee, ex, ey, ez);

  return hx[ex];
}

double Mesh_NURBS_1Patch_3D::get_hy( int ee) const
{
  int ex, ey, ez;
  get_elem_index(ee, ex, ey, ez);

  return hy[ey];
}

double Mesh_NURBS_1Patch_3D::get_hz( int ee) const
{
  int ex, ey, ez;
  get_elem_index(ee, ex, ey, ez);

  return hz[ez];
}

void Mesh_NURBS_1Patch_3D::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"====== Mesh_NURBS_1Patch_3D ======"<<std::endl;
  std::cout<<"Degree: "<<get_s_degree()<<'\t'<<get_t_degree()<<'\t'<<get_u_degree()<<std::endl;
  std::cout<<"Max mesh size: "<<get_hx_max()<<'\t'<<get_hy_max()<<'\t'<<get_hz_max()<<std::endl;
  std::cout<<"Min mesh size: "<<get_hx_min()<<'\t'<<get_hy_min()<<'\t'<<get_hz_min()<<std::endl;
  std::cout<<"Elem #: "<<get_nElem_x()<<'\t'<<get_nElem_y()<<'\t'<<get_nElem_z()<<std::endl;
  std::cout<<"Total Elem: "<<get_nElem()<<std::endl;
  std::cout<<"Func #: "<<get_nFunc_x()<<'\t'<<get_nFunc_y()<<'\t'<<get_nFunc_z()<<std::endl;
  std::cout<<"Total Func: "<<get_nFunc()<<std::endl;
  std::cout<<"Local Basis #: "<<get_nLocBas()<<std::endl;
  std::cout<<"Patch index: "<<get_patch_index()<<std::endl;
  std::cout<<"Start element: "<<get_nElem_start()<<std::endl;
  std::cout<<"Start function: "<<get_nFunc_start()<<std::endl;
  std::cout<<"=================================="<<std::endl;
}

// EOF
