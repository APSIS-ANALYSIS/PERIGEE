#include "ElemBC_3D_tet_wall.hpp"

ElemBC_3D_tet_wall::ElemBC_3D_tet_wall(
    const std::vector<std::string> &vtkfileList,
    const std::vector<double> thickness_to_radius,
    const int &elemtype )
: ElemBC_3D_tet( vtkfileList, elemtype )
{
  radius.resize(num_ebc);
  for(int ii=0; ii<num_ebc; ++ii) radius.reisze( num_node[ii] );

  const std::string centerlineFile("centerlines.vtp");

  vtkXMLPolyDataReader * reader = vtkXMLPolyDataReader::New();
  reader -> SetFileName( centerlineFile.c_str() );
  reader -> Update();

  line_polydata = reader -> GetOutput();
  reader -> Delete();

}


ElemBC_3D_tet_wall::~ElemBC_3D_tet_wall()
{
  VEC_T::clean( radius );
}

// EOF
