#include "AGlobal_Mesh_Info_TS_2D.hpp"


AGlobal_Mesh_Info_TS_2D::AGlobal_Mesh_Info_TS_2D( const HDF5_PartReader * const &h5reader )
{
  h5reader->get_GMI_nElem(nElem);
  h5reader->get_GMI_nFunc(nFunc);
  h5reader->get_GMI_probDim(probDim);
  h5reader->get_EXT_elemType(elemType);
}


AGlobal_Mesh_Info_TS_2D::~AGlobal_Mesh_Info_TS_2D()
{
}


void AGlobal_Mesh_Info_TS_2D::print() const
{
  SYS_T::cmdPrint("nElem", nElem);
  SYS_T::cmdPrint("nFunc", nFunc);
  SYS_T::cmdPrint("probDim", probDim);
  SYS_T::cmdPrint("elemType", elemType);
}




// EOF
