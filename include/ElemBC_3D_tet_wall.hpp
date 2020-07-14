#ifndef ELEMBC_3D_TET_WALL_HPP
#define ELEMBC_3D_TET_WALL_HPP
// ==================================================================
// ElemBC_3D_tet_wall.hpp
//
// A derived class from the ElemBC_3D_tet.hpp
//
// This class has additional information of the wall mesh.
//
// Date: July 12 2020
// ==================================================================
#include "ElemBC_3D_tet.hpp"
#include "Vector_3.hpp"

class ElemBC_3D_tet_wall : public ElemBC_3D_tet
{
  public:
    ElemBC_3D_tet_wall( const std::vector<std::string> &vtkfileList,
        const std::vector<double> thickness_to_radius,
        const std::string &centerlineFile = "centerlines.vtp",
        const int &elemtype = 501 );

    virtual ~ElemBC_3D_tet_wall();

  private:
    // num_ebc times num_node[ii] in size
    std::vector< std::vector<double> > radius;
};

#endif
