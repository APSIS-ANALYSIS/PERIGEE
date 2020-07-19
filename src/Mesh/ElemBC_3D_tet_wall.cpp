#include "ElemBC_3D_tet_wall.hpp"

ElemBC_3D_tet_wall::ElemBC_3D_tet_wall(
    const std::vector<std::string> &vtkfileList,
    const std::vector<double> &thickness_to_radius,
    const std::string &centerlineFile,
    const int &elemtype )
: ElemBC_3D_tet( vtkfileList, elemtype )
{
  // Check input size
  SYS_T::print_fatal_if( thickness_to_radius.size() != vtkfileList.size(),
      "Error: thickness_to_radius length does not match that of the vtkfileList.\n");

  radius.resize(num_ebc);
  thickness.resize(num_ebc);
  for(int ii=0; ii<num_ebc; ++ii) 
  {
    radius[ii].resize( num_node[ii] );
    thickness[ii].resize( num_node[ii] );
  }

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
      const double pt[3] = {pt_xyz[ebc_id][3*ii], pt_xyz[ebc_id][3*ii+1], pt_xyz[ebc_id][3*ii+2]};

      const int closest_id = locator -> FindClosestPoint(&pt[0]);

      const double * cl_pt = centerlineData -> GetPoints() -> GetPoint(closest_id);

      radius[ebc_id][ii] = MATH_T::norm2(cl_pt[0] - pt[0], cl_pt[1] - pt[1], cl_pt[2] - pt[2]);
   
      thickness[ebc_id][ii] = radius[ebc_id][ii] * thickness_to_radius[ebc_id]; 
    }
  }

  // clean memory
  locator -> Delete();
  reader -> Delete();
}


ElemBC_3D_tet_wall::~ElemBC_3D_tet_wall()
{
  for(int ii=0; ii<num_ebc; ++ii) VEC_T::clean( radius[ii] );

  VEC_T::clean( radius );
}


void ElemBC_3D_tet_wall::print_info() const
{
  ElemBC_3D_tet::print_info();

  for(int face=0; face<num_ebc; ++face)
  {
    VEC_T::print( radius[face], "wall_id_" + SYS_T::to_string(face) + "_radius.txt", '\n');
  }
}

// EOF
