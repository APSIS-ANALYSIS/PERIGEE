#ifndef FEAELEMENT_HEX8_HPP
#define FEAELEMENT_HEX8_HPP
// ==================================================================
// FEAElement_Hex8.hpp
// Element routine for the linear hexagon element, with evaluation
// of shape functions and their derivatives.
// 
// Hex8 means 8-node hex, aka trilinear hex.
//
//                    t
//                    ^
//                    |
//                    4------------------7
//                   /.                 /|
//                  / .                / |
//                 /  .               /  |
//                /   .              /   |
//               /    .             /    |
//              /     .            /     |
//             5------------------6      |
//             |      0...........|......3-------> s
//             |     .            |     /
//             |    .             |    /
//             |   .              |   /
//             |  .               |  /
//             | .                | /
//             |.                 |/
//             1------------------2
//            /
//           *
//           r
//
// Date Created: Sep 6 2023
// ==================================================================
#include "FEAElement.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FE_Tools.hpp"

class FEAElement_Hex8 : public FEAElement
{
  public :
    FEAElement_Hex8( const int &in_nqua );

    virtual ~FEAElement_Hex8();

    virtual int get_elemDim() const {return 3;}

    // A unique number for this element.
    virtual int get_Type() const {return 601;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 8;}

    virtual void print_info() const;

    virtual double get_memory_usage() const;

    // Given the quadrature points and nodal coordinates, evaluate the basis 
    // functions and their derivatives up to second order
    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );
    
    virtual double get_h( const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const;
    
    // Get functions give access to function evaluations at the quadrature point 
    // corresponding to quaindex
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

    // Get the Jacobian matrix dx/dr
    virtual void get_Jacobian(const int &quaindex, double * const &jac_value) const;

    virtual std::array<double,9> get_Jacobian( const int &quaindex ) const;

    // Get the inverse Jacobian matrix dr/dx
    virtual void get_invJacobian(const int &quaindex, double * const &jac_value) const;

    virtual std::array<double,9> get_invJacobian( const int &quaindex ) const;

    // Get the determinant of the Jacobian matrix
    virtual double get_detJac(const int &quaindex) const {return detJac[quaindex];}

    // Build basis and build the boundary element
    //   Hex-Face-0 : Node 0 3 2 1
    //   Hex-Face-1 : Node 4 5 6 7
    //   Hex-Face-2 : Node 0 1 5 4
    //   Hex-Face-3 : Node 1 2 6 5
    //   Hex-Face-4 : Node 3 7 6 2
    //   Hex-Face-5 : Node 0 4 7 3
    virtual void buildBasisBoundary( const IQuadPts * const &quad_rule_s, const int &face_id,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );

    // Get the outwardnormal on faces after calling buildBasisBoundary
    virtual Vector_3 get_2d_normal_out( const int &quaindex, double &area ) const
    {return quadrilateral_face->get_2d_normal_out( quaindex, area );}

  private:
    // Number of quadrature points
    const int numQuapts;

    // R : 0 <= ii < 8 x numQuapts
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

    FEAElement * quadrilateral_face;
    bool face_built_flag;
};

#endif
