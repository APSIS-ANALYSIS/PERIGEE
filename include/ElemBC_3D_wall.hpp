#ifndef ELEMBC_3D_WALL_HPP
#define ELEMBC_3D_WALL_HPP
// ============================================================================
// ElemBC_3D_wall.hpp
//
// A derived class from the ElemBC_3D.hpp
//
// This class has additional information of the wall mesh used for the thin-wall
// approximation in FSI simulations.
//
// Author: Ju Liu, Ingrid S. Lan
// Date: July 12 2020
// ============================================================================
#include "ElemBC_3D.hpp"
#include "vtkCellLocator.h"
#include "vtkGenericCell.h"

class ElemBC_3D_wall : public ElemBC_3D
{
  public:
    // ------------------------------------------------------------------------
    // Constructing an empty wall
    // ------------------------------------------------------------------------
    ElemBC_3D_wall( const int &elemtype );
    
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
    //       model without feeding a centerline vtp file.
    // ------------------------------------------------------------------------
    ElemBC_3D_wall(const std::string &walls_combined,
        const double &uniform_thickness,
        const double &uniform_youngsmod,
        const double &uniform_springconst,
        const double &uniform_dampingconst,
        const int &elemtype);

    // ------------------------------------------------------------------------
    //  Constructing wall properties with a single spatial distribution.
    //  \para: walls_combined contains a single vtp with the complete wall surface
    //  \para: centerlines_combined is a vtp with the complete set of centerlines 
    //  \para: thickness2radius_combined is the thickness-to-radius ratio
    //         to be prescribed for the complete wall. For most arteries, 
    //         we can assume the thickness is ten percent of the diameter,
    //         or twenty percent of the radius.
    //  \para: ks_combined is the  spring constant over the entire wall.
    //  \para: cs_combined is the damping constant over the entire wall.
    // ------------------------------------------------------------------------
    ElemBC_3D_wall(const std::string &walls_combined,
        const std::string &centerlines_combined,
        const double &thickness2radius_combined,
        const double &springconst_combined,
        const double &dampingconst_combined,
        const int &elemtype);

    // ------------------------------------------------------------------------
    //  Constructing wall properties with multiple spatial distributions.
    //  The background wall properties will first be prescribed using the
    //  constructor above. Wall properties in wallsList will then be overwritten
    //  with the corresponding centerlinesList and lists of properties, so
    //  all five lists must have the same length.
    //  \para: wallsList is a vector of wall surface vtp's, each a subset of
    //         the entire wall 
    //  \para: centerlinesList is a vector of corresponding centerline vtp's
    //  \para: thickness2radiusList is a vector of corresponding ratios.
    //  \para: ksList is a vector of corresponding spring constants.
    //  \para: csList is a vector of corresponding damping constants.
    // ------------------------------------------------------------------------
    ElemBC_3D_wall(const std::string &walls_combined,
        const std::string &centerlines_combined,
        const double &thickness2radius_combined,
        const double &springconst_combined,
        const double &dampingconst_combined,
        const std::vector<std::string> &wallsList,
        const std::vector<std::string> &centerlinesList,
        const std::vector<double> &thickness2radiusList,
        const std::vector<double> &springconstList,
        const std::vector<double> &dampingconstList,
        const int &elemtype);

    virtual ~ElemBC_3D_wall() = default;

    // Overwrite wall properties from a vtp/vtu file by GlobalNodeID. All wall nodes
    // with GlobalNodeIDs that also exist in wallprop_vtk will be overwritten. 
    // This can be used to apply wall properties from svPre's Laplacian solver, which
    // generates a wallprop vtp consisting of all caps in addition to the wall. 
    // NOTE: The GlobalNodeIDs are assumed to be consistent (zero-indexed).
    // \para: wallprop_vtk is the vtp / vtu file with the property field 
    // \para: type - 0 (thickness), 1 (youngsmod), 2 (springconst), 3 (dampingconst)
    // \para: vtk_fieldname is the corresponding field name in wallprop_vtk 
    virtual void overwrite_from_vtk(const std::string &wallprop_vtk,
        const int &type, const std::string &vtk_fieldname );

    virtual void print_info() const;

    virtual std::vector<double> get_wall_thickness() const {return thickness;}

    virtual std::vector<double> get_wall_youngsmod() const {return youngsmod;}

    virtual std::vector<double> get_wall_springconst() const {return springconst;}

    virtual std::vector<double> get_wall_dampingconst() const {return dampingconst;}

    virtual void write_vtk( const int &ebc_id,
        const std::string &filename="elembc_surface" ) const;

  private:
    // num_ebc = 1 for the wall, so these properties all have length num_node[0]
    std::vector<double> radius, thickness, youngsmod, springconst, dampingconst;
    
    // compute young's modulus for a wall node with the given radius & thickness
    // One may modify this function for different ways of prescribing the 
    // Young's modulus. The default one is adopted from N. Xiao, JCP, 244:22-40.
    void compute_youngsmod( const double &r, const double &th, double &E );

    // helper function for write_vtk to add wall properties to vtkPointSet
    void add_wall_data( vtkPointSet * const &grid_w, const int &ebc_id ) const;
};

#endif
