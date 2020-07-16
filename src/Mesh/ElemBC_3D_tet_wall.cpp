#include "ElemBC_3D_tet_wall.hpp"

ElemBC_3D_tet_wall::ElemBC_3D_tet_wall(
    const std::vector<std::string> &vtkfileList,
    const std::vector<double> &thickness_to_radius,
    const std::string &centerlineFile,
    const int &elemtype )
: ElemBC_3D_tet( vtkfileList, elemtype )
{
  radius.resize(num_ebc);
  for(int ii=0; ii<num_ebc; ++ii) radius.resize( num_node[ii] );

  vtkXMLPolyDataReader * reader = vtkXMLPolyDataReader::New();
  reader -> SetFileName( centerlineFile.c_str() );
  reader -> Update();

  vtkPolyData * centerlineData = reader -> GetOutput();
  
  vtkPointLocator * locator = vtkPointLocator::New();
  locator -> Initialize();
  locator -> SetDataSet( centerlineData );
  locator -> BuildLocator();

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    for(int ii=0; ii<num_node[ebc_id]; ++ii)
    {
      const double coor_x = pt_xyz[ebc_id][3*ii];
      const double coor_y = pt_xyz[ebc_id][3*ii+1];
      const double coor_z = pt_xyz[ebc_id][3*ii+2];

      const double pt[3] = {coor_x, coor_y, coor_z};

      const int closest_id = locator -> FindClosestPoint(&pt[0]);

      const double * closest_cl_pt = centerlineData -> GetPoints() -> GetPoint(closest_id);
      const double line_pt_x = closest_cl_pt[0];
      const double line_pt_y = closest_cl_pt[1];
      const double line_pt_z = closest_cl_pt[2];

      radius[ebc_id][ii] = MATH_T::norm2(line_pt_x - coor_x, line_pt_y - coor_y, line_pt_z - coor_z);
    }
  }

  // clean memory
  locator -> Delete();
  reader -> Delete();
}


ElemBC_3D_tet_wall::~ElemBC_3D_tet_wall()
{
  VEC_T::clean( radius );
}

// EOF
