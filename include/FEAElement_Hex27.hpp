#ifndef FEAELEMENT_HEX27_HPP
#define FEAELEMENT_HEX27_HPP
// ==================================================================
// FEAElement_Hex27.hpp
// Element routine for the quadratic hexagon element, with evaluation
// of shape functions and their derivatives.
// 
// Hex27 means 27-node hex, aka triquadratic hex.
//
//                    t
//                    ^
//                    |
//                    4--------15--------7
//                   /.                 /|
//                  / .                / |
//                 12 .               14 |
//                /   16             /   19
//               /    .             /    |
//              /     .            /     |
//             5--------13--------6      |
//             |      0........11.|......3-------> s
//             |     .            |     /
//             |    .             |    /
//             17  8              18  10
//             |  .               |  /
//             | .                | /
//             |.                 |/
//             1--------9---------2
//            /
//           *
//           r
//
//                    t
//                    ^
//                    |
//                    4------------------7
//                   /.                 /|
//                  / .                / |
//                 /  .        20     /  |
//                /   .        .     /   |
//               /    .      back   /    |
//              /     .            /  23 |
//             5------------------6   *  |
//             |   .  0...........|......3-------> s
//             |  22 .            |     /
//             |    .    21       |    /
//             |   .     *        |   /
//             |  .      front    |  /
//             | .                | /
//             |.                 |/
//             1------------------2
//            /
//           *
//           r
//
//                    t
//                    ^
//                    |
//                    4------------------7
//                   /.                 /|
//                  / .      25        / |
//                 /  .      *        /  |
//                /   .      up      /   |
//               /    .             /    |
//              /     .  26 # center     |
//             5------------------6      |
//             |      0...........|......3-------> s
//             |     .            |     /
//             |    .      24     |    /
//             |   .       .      |   /
//             |  .       down    |  /
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
#include "FEAElement_Quad9_3D_der0.hpp"
#include "FE_Tools.hpp"

class FEAElement_Hex27 final : public FEAElement
{
  public :
    FEAElement_Hex27( const int &in_nqua );

    ~FEAElement_Hex27() override = default;

    int get_elemDim() const override {return 3;}

    // A unique number for this element.
    FEType get_Type() const override {return FEType::Hex27;}

    int get_numQuapts() const override {return numQuapts;}

    int get_nLocBas() const override {return 27;}

    void print_info() const override;

    // Given the quadrature points and nodal coordinates, evaluate the basis 
    // functions and their derivatives up to second order
    void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) override;
    
    double get_h( const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const override;
    
    // Get functions give access to function evaluations at the quadrature point 
    // corresponding to quaindex
    void get_R( const int &quaindex, double * const &basis ) const override;

    std::vector<double> get_R( const int &quaindex ) const override;

    void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y, double * const &basis_z ) const override;

    void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y,
        double * const &basis_z ) const override;

    void get_3D_R_dR_d2R( const int &quaindex, 
        double * const &basis, double * const &basis_x, 
        double * const &basis_y, double * const &basis_z,
        double * const &basis_xx, double * const &basis_yy, 
        double * const &basis_zz, double * const &basis_xy, 
        double * const &basis_xz, double * const &basis_yz ) const override;

    void get_3D_R_gradR_LaplacianR( const int &quaindex,
        double * const &basis, double * const &basis_x, 
        double * const &basis_y, double * const &basis_z, 
        double * const &basis_xx, double * const &basis_yy, 
        double * const &basis_zz ) const override;
    
    // Get the Jacobian matrix dx/dr
    void get_Jacobian(const int &quaindex, double * const &jac_value) const override;

    std::array<double,9> get_Jacobian( const int &quaindex ) const override;

    // Get the inverse Jacobian matrix dr/dx
    void get_invJacobian(const int &quaindex, double * const &jac_value) const override;

    std::array<double,9> get_invJacobian( const int &quaindex ) const override;

    // Get the determinant of the Jacobian matrix
    double get_detJac(const int &quaindex) const override {return detJac[quaindex];}

    // Build basis and build the boundary element
    //   Hex-Face-0 : Node 0 3 2 1 11 10 9 8 24
    //   Hex-Face-1 : Node 4 5 6 7 12 13 14 15 25
    //   Hex-Face-2 : Node 0 1 5 4 8 17 12 16 22
    //   Hex-Face-3 : Node 1 2 6 5 9 18 13 17 21
    //   Hex-Face-4 : Node 3 7 6 2 19 14 18 10 23
    //   Hex-Face-5 : Node 0 4 7 3 16 15 19 11 20
    void buildBasis( const int &face_id, const IQuadPts * const &quad_rule_s,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) override;

    // Get the outwardnormal on faces
    Vector_3 get_2d_normal_out( const int &quaindex, double &area ) const override
    {return quadrilateral_face->get_2d_normal_out( quaindex, area );}

  private:
    // Number of quadrature points
    const int numQuapts;

    // R : 0 <= ii < 27 x numQuapts
    std::vector<double> R{}, dR_dx{}, dR_dy{}, dR_dz {},
        d2R_dxx{}, d2R_dyy{}, d2R_dzz{}, d2R_dxy{}, d2R_dxz{}, d2R_dyz{};

    // Container for
    // dx_dr : 0 <= ii < 9 numQuapts
    std::vector<double> dx_dr {};

    // dr_dx : 0 <= ii < 9 numQuapts
    std::vector<double> dr_dx {};

    // detJac : 0 <= ii < numQuapts
    std::vector<double> detJac {};

    std::unique_ptr<FEAElement> quadrilateral_face;
};

#endif
