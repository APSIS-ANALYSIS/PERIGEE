#ifndef FEAELEMENT_TET10_HPP
#define FEAELEMENT_TET10_HPP
// ==================================================================
// FEAElement_Tet10.hpp
// This is an implementation of the element routine for quadratic 10
// node tetrahedral element.
//
// Tet10 : 10-node tet element, aka, quadratic tet.
//
//                     t
//                     ^
//                     |
//                     3
//                    /| -
//                     |   -
//                  /  |     -
//                     7        8
//               /    |           -
//                9    |             -
//                     |                - 
//              /    /0--------6---------2-------> s
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
// Date created: Feb. 16 2018
// ==================================================================
#include "FEAElement.hpp"
#include "FE_Tools.hpp"

class FEAElement_Tet10 : public FEAElement
{
  public:
    FEAElement_Tet10( const int &in_nqua );

    virtual ~FEAElement_Tet10();

    virtual int get_elemDim() const {return 3;}

    virtual FEType get_Type() const {return FEType::Tet10;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 10;}

    virtual void print_info() const;

    virtual double get_memory_usage() const;

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

    // Get functions that give access to basis function values
    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual std::vector<double> get_R( const int &quaindex ) const;

    virtual void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y, double * const &basis_z ) const;

    virtual void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y,
        double * const &basis_z ) const;

    virtual std::vector<double> get_dR_dx( const int &quaindex ) const;

    virtual std::vector<double> get_dR_dy( const int &quaindex ) const;

    virtual std::vector<double> get_dR_dz( const int &quaindex ) const;

    virtual void get_3D_R_dR_d2R( const int &quaindex,
        double * const &basis, double * const &basis_x,
        double * const &basis_y, double * const &basis_z,
        double * const &basis_xx, double * const &basis_yy,
        double * const &basis_zz, double * const &basis_xy,
        double * const &basis_xz, double * const &basis_yz ) const;

    virtual void get_3D_R_gradR_LaplacianR( const int &quaindex,
        double * const &basis, double * const &basis_x,
        double * const &basis_y, double * const &basis_z,
        double * const &basis_xx, double * const &basis_yy,
        double * const &basis_zz ) const;

    virtual std::vector<double> get_d2R_dxx( const int &quaindex ) const;

    virtual std::vector<double> get_d2R_dyy( const int &quaindex ) const;

    virtual std::vector<double> get_d2R_dzz( const int &quaindex ) const;

    virtual std::vector<double> get_d2R_dxy( const int &quaindex ) const;

    virtual std::vector<double> get_d2R_dxz( const int &quaindex ) const;

    virtual std::vector<double> get_d2R_dyz( const int &quaindex ) const;

    virtual void get_Jacobian(const int &quaindex, double * const &jac_value) const;

    virtual std::array<double,9> get_Jacobian( const int &quaindex ) const;

    virtual void get_invJacobian(const int &quaindex, double * const &jac_value) const;

    virtual std::array<double,9> get_invJacobian( const int &quaindex ) const;
    
    virtual double get_detJac(const int &quaindex) const {return detJac[quaindex];}

  private:
    const int numQuapts;

    // R: 0 <= ii < 10 numQuapts
    double * R, * dR_dx, * dR_dy, * dR_dz;
    double * d2R_dxx, * d2R_dyy, * d2R_dzz;
    double * d2R_dxy, * d2R_dxz, * d2R_dyz;

    // Container for
    // dx_dr : 0 <= ii < 9 numQuapts
    double * dx_dr;

    // dr_dx : 0 <= ii < 9 numQuapts
    double * dr_dx;

    // detJac : 0 <= ii < numQuapts
    double * detJac;
};

#endif
