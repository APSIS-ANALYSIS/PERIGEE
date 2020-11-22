#ifndef FEAELEMENT_TRIANGLE6_3D_DER0_HPP
#define FEAELEMENT_TRIANGLE6_3D_DER0_HPP
// ==================================================================
// FEAElement_Triangle6_3D_der0.hpp
// Element routine for the quadratic triangular element in 
// three-dimensional space, with evaluation of the basis functions
// only (no derivatives).
//
//     s
//     |
//     2
//     | -
//     |   -
//     |     -
//     |       -      t = 1 - r - s
//     5         4
//     |           -
//     |             -
//     |               -
//     0---------3------- 1 -- r
//
//
// Triangle6 means 6-node triangle; _3D means the element has
// coordinates in three-dimensions; _der0 means that only the function
// itself is evaluated.
//
// This class is designed for boundary integrations for elemental/
// natural boundary conditions in 3D problems.
//
// Author: Ju Liu
// Date Created: Feb. 17 2018.
// ==================================================================
#include "FEAElement.hpp"

class FEAElement_Triangle6_3D_der0 : public FEAElement
{
  public:
    FEAElement_Triangle6_3D_der0( const int &in_nqua );

    virtual ~FEAElement_Triangle6_3D_der0();

    virtual int get_elemDim() const {return 2;}

    // element type : 522
    // 5: simplicial element
    // 2: 2D element
    virtual int get_Type() const {return 522;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_numType() const {return 1;}

    virtual int get_nLocBas() const {return 6;}

    virtual void print_info() const;

    virtual double get_memory_usage() const;

    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );

    virtual void get_R( const int &quaindex, double * const &basis ) const;

    // Assuming the triangle nodes are arranged such that the outward
    // direction is given by dx_dr x dx_ds
    virtual void get_2d_normal_out( const int &quaindex,
        double &nx, double &ny, double &nz, double &len ) const;

    // If the triangle nodes are NOT arranged in any particular order,
    // use an interior node to define the outward direction. 
    virtual void get_normal_out( const int &quaindex,
        const double &sur_pt_x, const double &sur_pt_y, const double &sur_pt_z,
        const double &intpt_x, const double &intpt_y, const double &intpt_z,
        double &nx, double &ny, double &nz, double &len ) const;

    virtual double get_detJac(const int &quaindex) const
    {return detJac[quaindex];}

  private:
    const int numQuapts;

    // Container for R0, R1, R2, R3, R4, R5
    // 0 <= ii < 6 x numQuapts
    double * R;

    // Container for dx_dr, etc., each of length numQuapts
    double * dx_dr, * dx_ds;
    double * dy_dr, * dy_ds;
    double * dz_dr, * dz_ds;

    // unit normal vector components, each of length numQuapts
    double * unx, * uny, * unz;

    // Jacobian determinant, length numQuapts
    double * detJac; 
};

#endif
