#include "Mesh_NURBS_1Patch_3D_nze.hpp"

Mesh_NURBS_1Patch_3D_nze::Mesh_NURBS_1Patch_3D_nze( IMesh const * const &inmesh )
{
  s_degree = inmesh->get_s_degree();
  t_degree = inmesh->get_t_degree();
  u_degree = inmesh->get_u_degree();

  hx_max   = inmesh->get_hx_max();
  hy_max   = inmesh->get_hy_max();
  hz_max   = inmesh->get_hz_max();

  hx_min   = inmesh->get_hx_min();
  hy_min   = inmesh->get_hy_min();
  hz_min   = inmesh->get_hz_min();

  nFunc_x  = inmesh->get_nFunc_x();
  nFunc_y  = inmesh->get_nFunc_y();
  nFunc_z  = inmesh->get_nFunc_z();
  nFunc    = inmesh->get_nFunc();

  nElem_x  = inmesh->get_nElem_x_nz();
  nElem_y  = inmesh->get_nElem_y_nz();
  nElem_z  = inmesh->get_nElem_z_nz();
  nElem    = inmesh->get_nElem_nz();

  nLocBas = inmesh->get_nLocBas();

  hx.clear(); hy.clear(); hz.clear();

  double temp = 0.0;
  
  for(int ii=0; ii<inmesh->get_nElem_x(); ++ii)
  {
    temp = inmesh->get_hx(ii);
    if(temp > 0.0) hx.push_back(temp);
  }

  for(int ii=0; ii<inmesh->get_nElem_y(); ++ii)
  {
    temp = inmesh->get_hy( ii * inmesh->get_nElem_x() ); 
    if(temp > 0.0) hy.push_back(temp);
  }

  for(int ii=0; ii<inmesh->get_nElem_z(); ++ii)
  {
    temp = inmesh->get_hz( ii * inmesh->get_nElem_x() * inmesh->get_nElem_y() );
    if(temp > 0.0) hz.push_back(temp);
  }

  VEC_T::shrink2fit(hx);
  VEC_T::shrink2fit(hy);
  VEC_T::shrink2fit(hz);

  if( int(hx.size()) != nElem_x )
  {
    std::cerr<<"ERROR: hx.size ="<<hx.size()<<" nElem_x = "<<nElem_x<<std::endl;
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  if( int(hy.size()) != nElem_y )
  {
    std::cerr<<"ERROR: hy.size ="<<hy.size()<<" nElem_y = "<<nElem_y<<std::endl;
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  if( int(hz.size()) != nElem_z )
  {
    std::cerr<<"ERROR: hz.size ="<<hz.size()<<" nElem_z = "<<nElem_z<<std::endl;
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
}



Mesh_NURBS_1Patch_3D_nze::~Mesh_NURBS_1Patch_3D_nze()
{
  std::cout<<"-- Mesh_NURBS_1Patch_3D_nze deleted. \n";
}



void Mesh_NURBS_1Patch_3D_nze::print_info() const
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


double Mesh_NURBS_1Patch_3D_nze::get_hx(int ee) const
{
  int ex, ey, ez;
  get_elem_index(ee, ex, ey, ez);
  return hx[ex];
}


double Mesh_NURBS_1Patch_3D_nze::get_hy(int ee) const
{
  int ex, ey, ez;
  get_elem_index(ee, ex, ey, ez);
  return hy[ey];
}


double Mesh_NURBS_1Patch_3D_nze::get_hz(int ee) const
{
  int ex, ey, ez;
  get_elem_index(ee, ex, ey, ez);
  return hz[ez];
}


void Mesh_NURBS_1Patch_3D_nze::get_elem_index( const int &ee, int &ex, 
    int &ey, int &ez ) const
{
  int exy = ee % (nElem_x * nElem_y);
  ez = (ee - exy) / (nElem_x * nElem_y);
  ex = exy % nElem_x;
  ey = (exy - ex) / nElem_x;
}

// EOF
