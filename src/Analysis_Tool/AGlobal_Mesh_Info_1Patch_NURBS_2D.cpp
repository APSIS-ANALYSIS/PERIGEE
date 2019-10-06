#include "AGlobal_Mesh_Info_1Patch_NURBS_2D.hpp"

AGlobal_Mesh_Info_1Patch_NURBS_2D::AGlobal_Mesh_Info_1Patch_NURBS_2D(
    const HDF5_PartReader * const &h5reader )
{
  h5reader->get_GMI_degree(xdegree, ydegree);
  h5reader->get_GMI_h_max(hx_max, hy_max);
  h5reader->get_GMI_h_min(hx_min, hy_min);
  h5reader->get_GMI_nElem(nElem, nElem_x, nElem_y);
  h5reader->get_GMI_nFunc(nFunc, nFunc_x, nFunc_y);
  h5reader->get_GMI_nLocBas(nLocBas);
  h5reader->get_GMI_probDim(probDim);
  h5reader->get_EXT_elemType(elemType);

  if(probDim != 2)
  {
    SYS_T::commPrint("Error: The problem dimension is not 2. ");
    SYS_T::commPrint(" The Global Mesh Info class is for 2D single patch NURBS geometry. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
}

AGlobal_Mesh_Info_1Patch_NURBS_2D::~AGlobal_Mesh_Info_1Patch_NURBS_2D()
{}

void AGlobal_Mesh_Info_1Patch_NURBS_2D::print() const
{
  std::cout<<"AGlobal_Mesh_Info_1Patch_NURBS_2D:"<<std::endl;
  std::cout<<"degree: "<<xdegree<<'\t'<<ydegree<<'\n';
  std::cout<<"h_max: "<<hx_max<<'\t'<<hy_max<<'\n';
  std::cout<<"h_min: "<<hx_min<<'\t'<<hy_min<<'\n';
  std::cout<<"nElem: "<<nElem<<'\t'<<nElem_x<<'\t'<<nElem_y<<'\n';
  std::cout<<"nFunc: "<<nFunc<<'\t'<<nFunc_x<<'\t'<<nFunc_y<<'\n';
  std::cout<<"nLocBas: "<<nLocBas<<std::endl;
  std::cout<<"probDim: "<<probDim<<std::endl;
  std::cout<<"elemType: "<<elemType<<std::endl;
}

// EOF
