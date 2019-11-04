#ifndef FEAELEMENT_TET10_V2_HPP
#define FEAELEMENT_TET10_V2_HPP
// ==================================================================
// FEAElement_Tet10_v2.hpp
// This is an implementation of the element routine for quadratic 10
// node tetrahedral element. Version 2 means this is a second version
// implementation of the Tet10 element to make the node numbering
// compatible with the vtk format, see the graph below.
// 
// Notice that only the node number 8 and 9 are switched. The rest
// numbering are identical to that of FEAElement_Tet10 class.
//
// Tet10 : 10-node tet element, aka, quadratic tet.
//
//                     u
//                     ^
//                     |
//                     3
//                    /| -
//                     |   -
//                  /  |     -
//                     7        9
//                /    |           -
//               8     |             -
//                     |                - 
//               /    /0--------6---------2-------> s
//              /    /                -
//              /   /              -
//                 4           -   
//             /  /        5    
//             / /     -
//              /  - 
//             1
//            /
//           *
//           r
//
// This class is designed for the volumetric integration in model
// assembly.
//
// Date created: Nov. 3 2019
// ==================================================================
#include "FEAElement.hpp"
#include "Matrix_double_6by6_Array.hpp"

class FEAElement_Tet10_v2 : public FEAElement
{
  public:
    FEAElement_Tet10_v2( const int &in_nqua );

    virtual ~FEAElement_Tet10_v2();

    virtual int get_elemDim() const {return 3;}

    // A unique number for this element. 
    // Recall that 532 is the original Tet10 element
    virtual int get_Type() const {return 533;}

    virtual int get_numType() const {return 1;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 10;}

    virtual void print_info() const;

    virtual double get_memory_usage() const;

    // Given the quadrature points and nodes' coordinates, evaluate
    // the basis functions and their derivatives up to second order.
    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );

    // Return the element size.
    // Here we adopt the algorithm for Tet4, and use the four vertex
    // nodes to calculate the element size
    virtual double get_h( const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const;

    // get_xxx functions give access to the basis functions' value
    // at the corresponding quadrature points
    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y, double * const &basis_z ) const;

    virtual void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y,
        double * const &basis_z ) const;

    virtual void get_3D_R_dR_d2R( const int &quaindex,
        double * const &basis, double * const &basis_x,
        double * const &basis_y, double * const &basis_z,
        double * const &basis_xx, double * const &basis_yy,
        double * const &basis_zz, double * const &basis_xy,
        double * const &basis_xz, double * const &basis_yz ) const;

    virtual void get_Jacobian(const int &quaindex,
        double * const &jac_value) const;

    virtual void get_invJacobian(const int &quaindex,
        double * const &jac_value) const;

    virtual double get_detJac(const int &quaindex) const
    {return detJac[quaindex];}

  private:
    const int numQuapts;

    // Auxiliary container for calculating second derivatives
    double RHS [6];
    double sol [6];

    // R: 0 <= ii < 10 numQuapts
    double * R, * dR_dx, * dR_dy, * dR_dz;
    double * d2R_dxx, * d2R_dyy, * d2R_dzz;
    double * d2R_dxy, * d2R_dxz, * d2R_dyz;

    // Container for dR_ds/t/u and second derivative at quad pts.
    // Length is 10
    double dR_dr [10], dR_ds [10], dR_dt [10];
    double d2R_drr [10], d2R_dss [10], d2R_dtt [10];
    double d2R_drs [10], d2R_drt [10], d2R_dst [10]; 

    // Container for
    // dx_dr : 0 <= ii < 9 numQuapts
    double * dx_dr;
    
    // dr_dx : 0 <= ii < 9 numQuapts
    double * dr_dx;

    // detJac : 0 <= ii < numQuapts
    double * detJac;
};

#endif
