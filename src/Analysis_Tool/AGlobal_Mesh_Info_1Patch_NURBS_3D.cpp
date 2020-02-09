#include "AGlobal_Mesh_Info_1Patch_NURBS_3D.hpp"

AGlobal_Mesh_Info_1Patch_NURBS_3D::AGlobal_Mesh_Info_1Patch_NURBS_3D(
    const HDF5_PartReader * const &h5reader )
{
  h5reader->get_GMI_degree(xdegree, ydegree, zdegree);
  h5reader->get_GMI_h_max(hx_max, hy_max, hz_max);
  h5reader->get_GMI_h_min(hx_min, hy_min, hz_min);
  h5reader->get_GMI_nElem(nElem, nElem_x, nElem_y, nElem_z);
  h5reader->get_GMI_nFunc(nFunc, nFunc_x, nFunc_y, nFunc_z);
  h5reader->get_GMI_nLocBas(nLocBas);
  h5reader->get_GMI_probDim(probDim);
  h5reader->get_EXT_elemType(elemType);

  if(probDim != 3)
  {
    SYS_T::commPrint("Error: The problem dimension is not 3. ");
    SYS_T::commPrint(" The Global Mesh Info class is for 3D NURBS single patch problem. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
}

AGlobal_Mesh_Info_1Patch_NURBS_3D::~AGlobal_Mesh_Info_1Patch_NURBS_3D()
{}

void AGlobal_Mesh_Info_1Patch_NURBS_3D::print_info() const
{
  std::cout<<"AGlobal_Mesh_Info_1Patch_NURBS_3D:"<<std::endl;
  std::cout<<"degree: "<<xdegree<<'\t'<<ydegree<<'\t'<<zdegree<<'\n';
  std::cout<<"h_max: "<<hx_max<<'\t'<<hy_max<<'\t'<<hz_max<<'\n';
  std::cout<<"h_min: "<<hx_min<<'\t'<<hy_min<<'\t'<<hz_min<<'\n';
  std::cout<<"nElem: "<<nElem<<'\t'<<nElem_x<<'\t'<<nElem_y<<'\t'<<nElem_z<<'\n';
  std::cout<<"nFunc: "<<nFunc<<'\t'<<nFunc_x<<'\t'<<nFunc_y<<'\t'<<nFunc_z<<'\n';
  std::cout<<"nLocBas: "<<nLocBas<<std::endl;
  std::cout<<"probDim: "<<probDim<<std::endl;
  std::cout<<"elemType: "<<get_elemType()<<std::endl;
}

// EOF
