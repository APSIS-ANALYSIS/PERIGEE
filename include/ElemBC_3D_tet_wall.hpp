#ifndef ELEMBC_3D_TET_WALL_HPP
#define ELEMBC_3D_TET_WALL_HPP
// ============================================================================
// ElemBC_3D_tet_wall.hpp
//
// A derived class from the ElemBC_3D_tet.hpp
//
// This class has additional information of the wall mesh used for
// thin-wall approximation in FSI simulations.
//
// Author: Ju Liu, Ingrid S. Lan
// Date: July 12 2020
// ============================================================================
#include "ElemBC_3D_tet.hpp"

class ElemBC_3D_tet_wall : public ElemBC_3D_tet
{
  public:
    // ------------------------------------------------------------------------
    // Constructing an empty wall
    // ------------------------------------------------------------------------
    ElemBC_3D_tet_wall( const int &elemtype = 501, const double &in_fluid_density = 1.065 );
    
    // ------------------------------------------------------------------------
    // Constructing wall properties with uniform thickness
    // \para: walls_combined contains a single vtp with the complete wall surface
    // \para: uniform_thickness is the wall thickness
    // \para: uniform_youngsmod is the wall Young's modulus
    // \para: uniform_ks is the spring constant for external tissue support
    //        (Kelvin-Voigt viscoelastic model)
    // \para: uniform_cs is the damping constant for external tissue support
    //        (Kelvin-Voigt viscoelastic model)
    // Note: This function is typically used for constructing a simple wall
    //       model without feeding a centerline vtp file. The input of fluid
    //       density is logically not needed, but to keep output file
    //       consistency.
    // ------------------------------------------------------------------------
    ElemBC_3D_tet_wall(const std::string &walls_combined,
        const double &uniform_thickness,
        const double &uniform_youngsmod,
        const double &uniform_ks = 0.0,
        const double &uniform_cs = 0.0,
        const int &elemtype = 501,
        const double &in_fluid_density = 1.065 );

    // ------------------------------------------------------------------------
    //  Constructing wall properties with a single spatial distribution.
    //  \para: walls_combined contains a single vtp with the complete wall surface
    //  \para: centerlines_combined is a vtp with the complete set of centerlines 
    //  \para: thickness2radius_combined is the thickness-to-radius ratio
    //         to be prescribed for the complete wall. For most arteries, 
    //         we can assume the thickness is ten percent of the diameter,
    //         or twenty percent of the radius.
    // ------------------------------------------------------------------------
    ElemBC_3D_tet_wall(const std::string &walls_combined,
        const std::string &centerlines_combined,
        const double &thickness2radius_combined,
        const int &elemtype = 501,
        const double &in_fluid_density = 1.065 );

    // ------------------------------------------------------------------------
    //  Constructing wall properties with multiple spatial distributions.
    //  The background wall properties will first be prescribed using the
    //  constructor above. Wall properties in wallsList will then be overwritten
    //  with the corresponding centerlinesList and thickness2radiusList, so
    //  these three lists must have the same length.
    //  \para: wallsList is a vector of wall surface vtp's, each a subset of
    //         the entire wall 
    //  \para: centerlinesList is a vector of corresponding centerline vtp's
    //  \para: thickness2radiusList is a vector of corresponding ratios.
    // ------------------------------------------------------------------------
    ElemBC_3D_tet_wall(const std::string &walls_combined,
        const std::string &centerlines_combined,
        const double &thickness2radius_combined,
        const std::vector<std::string> &wallsList,
        const std::vector<std::string> &centerlinesList,
        const std::vector<double> &thickness2radiusList,
        const int &elemtype = 501,
        const double &in_fluid_density = 1.065 );

    virtual ~ElemBC_3D_tet_wall();

    virtual void print_info() const;

    virtual std::vector<double> get_wall_thickness() const {return thickness;}

    virtual std::vector<double> get_wall_youngsmod() const {return youngsmod;}

    virtual void write_vtk( const int &ebc_id,
        const std::string &filename="elembc_surface" ) const;

    virtual double get_fluid_density() const {return fluid_density;}

  private:
    // fluid density used to compute the Young's modulus.
    // Its value will be passed to the analysis code and will be checked to
    // ensure that a consistent value is used in analysis code.
    const double fluid_density;

    // num_ebc = 1 for the wall, so these properties all have length num_node[0]
    std::vector<double> radius, thickness, youngsmod, ks, cs;
    
    // compute young's modulus for a wall node with the given radius & thickness
    // One may modify this function for different ways of prescribing the 
    // Young's modulus. The default one is adopted from N. Xiao, JCP, 244:22-40.
    void compute_youngsmod( const double &r, const double &th, double &E );

    // helper function for write_vtk to add wall properties to vtkPointSet
    void add_wall_data( vtkPointSet * const &grid_w, const int &ebc_id ) const;
};

#endif
