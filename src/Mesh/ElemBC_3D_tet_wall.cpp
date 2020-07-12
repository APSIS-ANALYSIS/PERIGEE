#include "ElemBC_3D_tet_wall.hpp"

ElemBC_3D_tet_wall::ElemBC_3D_tet_wall(
    const std::vector<std::string> &vtkfileList,
    const int &elemtype )
: ElemBC_3D_tet( vtkfileList, elemtype )
{
  int total_num_nodes = 0;
  for(int ii=0; ii<num_ebc; ++ii) total_num_nodes += num_node[ii];

  radius.resize( total_num_nodes );
  radial.resize( total_num_nodes );
  tangential.resize( total_num_nodes );
  circumferential.resize( total_num_nodes );

}


ElemBC_3D_tet_wall::~ElemBC_3D_tet_wall()
{
  VEC_T::clean( radius ); VEC_T::clean( radial );
  VEC_T::clean( tangential ); VEC_T::clean( circumferential );
}

// EOF
