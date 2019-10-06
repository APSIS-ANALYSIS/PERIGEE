#include "ALocal_meshSize_3D_NURBS.hpp"

ALocal_meshSize_3D_NURBS::ALocal_meshSize_3D_NURBS(const HDF5_PartReader * const &h5reader)
{
  h5reader->get_hxyz(hx,hy,hz);
}

ALocal_meshSize_3D_NURBS::~ALocal_meshSize_3D_NURBS()
{}

void ALocal_meshSize_3D_NURBS::print() const
{
  std::cout<<"ALocal_meshSize_3D_NURBS: \n";
  std::cout<<"hx: "<<hx.size()<<'\n';
  VEC_T::print(hx);
  std::cout<<"\n hy "<<hy.size()<<'\n';
  VEC_T::print(hy);
  std::cout<<"\n hz "<<hz.size()<<'\n';
  VEC_T::print(hz);
}
