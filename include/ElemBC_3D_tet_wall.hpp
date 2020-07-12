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

class ElemBC_3D_tet_wall : public ElemBC_3D_tet
{
  public:
    ElemBC_3D_tet_outflow( const std::vector<std::string> &vtkfileList,
        const int &elemtype = 501 );

    virtual ~ElemBC_3D_tet_wall();

  private:
    std::vector<double> raidus;

    std::vector<Vector_3> radial, tangential, circumferential; 
};

#endif
