#ifndef ELEMBC_3D_TET_WALL_HPP
#define ELEMBC_3D_TET_WALL_HPP
// ==================================================================
// ElemBC_3D_tet_wall.hpp
//
// A derived class from the ElemBC_3D_tet.hpp
//
// This class has additional information of the wall mesh.
//
// Author: Ju Liu
// Date: July 12 2020
// ==================================================================
#include "ElemBC_3D_tet.hpp"
#include "Vector_3.hpp"

class ElemBC_3D_tet_wall : public ElemBC_3D_tet
{
  public:
    // Constructing the wall boundary conditions.
    // \para: vtkfileList stores the set of vtp files of the wall, the
    //        union of these vtp files is the whole wall surface. Its
    //        length is num_ebc.
    // \para: thickness_to_radius has the length num_ebc, which gives
    //        the thickness ratio to the radius at the corresponding
    //        wall surface vtp file. For most arteries, the ratio is
    //        ten percent to the diameter, that is twenty percent to
    //        the radius.
    ElemBC_3D_tet_wall( const std::vector<std::string> &vtkfileList,
        const std::vector<double> &thickness_to_radius,
        const std::vector<double> &youngsmod_alpha,
        const std::vector<double> &youngsmod_beta,
        const std::string &centerlineFile = "centerlines.vtp",
        const int &fluid_density = 1.065,
        const int &elemtype = 501 );

    virtual ~ElemBC_3D_tet_wall();

    virtual void print_info() const;

    virtual void get_wall_thickness( const int &ebc_id, std::vector<double> &th ) const
    {th = thickness[ebc_id];}

    virtual void get_wall_youngsmod( const int &ebc_id, std::vector<double> &E ) const
    {E = youngsmod[ebc_id];}

  private:
    // num_ebc times num_node[ii] in size, 0 <= ii <num_ebc
    // here num_ebc means the number of different wall regions, which
    // potentially have different properties (thickness, modulus, etc).
    std::vector< std::vector<double> > radius;
    
    // num_ebc times num_node[ii] in size, 0 <= ii <num_ebc
    std::vector< std::vector<double> > thickness;
    std::vector< std::vector<double> > youngsmod;

};

#endif
