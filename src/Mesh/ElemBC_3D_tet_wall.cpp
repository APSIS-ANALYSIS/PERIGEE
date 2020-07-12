#include "ElemBC_3D_tet_wall.hpp"

ElemBC_3D_tet_wall::ElemBC_3D_tet_wall(
    const std::vector<std::string> &vtkfileList,
    const int &elemtype )
: ElemBC_3D_tet( vtkfileList, elemtype )
{}


ElemBC_3D_tet_wall::~ElemBC_3D_tet_wall()
{
  VEC_T::clean( radius ); VEC_T::clean( radial );
  VEC_T::clean( tangential ); VEC_T::clean( circumferential );
}

// EOF
