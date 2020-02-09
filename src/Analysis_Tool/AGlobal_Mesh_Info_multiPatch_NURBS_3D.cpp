#include "AGlobal_Mesh_Info_multiPatch_NURBS_3D.hpp"

AGlobal_Mesh_Info_multiPatch_NURBS_3D::AGlobal_Mesh_Info_multiPatch_NURBS_3D( 
    const HDF5_PartReader * const &h5reader )
{
  h5reader->get_GMI_degree(xdegree, ydegree, zdegree);
  h5reader->get_GMI_h_max(hx_max, hy_max, hz_max);
  h5reader->get_GMI_h_min(hx_min, hy_min, hz_min);
  h5reader->get_GMI_nElem(nElem);
  h5reader->get_GMI_nFunc(nFunc);
  h5reader->get_GMI_nLocBas(nLocBas);
  h5reader->get_GMI_probDim(probDim);
  h5reader->get_EXT_elemType(elemType);

  if(probDim != 3)
  {
    SYS_T::commPrint("Error: The problem dimension is not 3. ");
    SYS_T::commPrint(" This global mesh info class is for 3D multipatch geometries. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

}


AGlobal_Mesh_Info_multiPatch_NURBS_3D::~AGlobal_Mesh_Info_multiPatch_NURBS_3D()
{}


void AGlobal_Mesh_Info_multiPatch_NURBS_3D::print_info() const
{
  std::cout<<"AGlobal_Mesh_Info_multiPatch_NURBS_3D:"<<std::endl;
  std::cout<<"degree: "<<xdegree<<'\t'<<ydegree<<'\t'<<zdegree<<'\n';
  std::cout<<"h_max: "<<hx_max<<'\t'<<hy_max<<'\t'<<hz_max<<'\n';
  std::cout<<"h_min: "<<hx_min<<'\t'<<hy_min<<'\t'<<hz_min<<'\n';
  std::cout<<"nElem: "<<nElem<<'\n';
  std::cout<<"nFunc: "<<nFunc<<'\n';
  std::cout<<"nLocBas: "<<nLocBas<<std::endl;
  std::cout<<"probDim: "<<probDim<<std::endl;
  std::cout<<"elemType: "<<get_elemType()<<std::endl;
}


// EOF
