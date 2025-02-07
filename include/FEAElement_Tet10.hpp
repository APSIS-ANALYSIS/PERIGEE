#ifndef FEAELEMENT_TET10_HPP
#define FEAELEMENT_TET10_HPP
// ==================================================================
// FEAElement_Tet10.hpp
// Element routine for quadratic 10-node tetrahedral element.
// The node numbering is made compatible with the vtk format -- see 
// the graph below.
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
//                     7        9
//                 /   |           -
//                     |             -
//                8    |                - 
//                    /0--------6---------2-------> s
//               /   /                -
//                  /              -
//              /  4           -   
//                /        5    
//             / /     -
//              /  - 
//             1
//            /
//           *
//           r
//
// This class is designed for volumetric integration in model assembly.
//
// Date created: Nov. 3 2019
// ==================================================================
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FE_Tools.hpp"

class FEAElement_Tet10 final : public FEAElement
{
  public:
    FEAElement_Tet10( const int &in_nqua );

    ~FEAElement_Tet10() override = default;

    int get_elemDim() const override {return 3;}

    FEType get_Type() const override {return FEType::Tet10;}

    int get_numQuapts() const override {return numQuapts;}

    int get_nLocBas() const override {return nLocBas;}

    void print_info() const override;

    // Given the quadrature points and nodal coordinates, evaluate
    // the basis functions and their derivatives up to second order
    void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) override;

    // Return the element size.
    // Here we adopt the algorithm for Tet4 and use the four vertex
    // nodes to calculate the element size
    double get_h( const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const override;

    // get_xxx functions give access to function evaluations at the
    // quadrature point corresponding to quaindex
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

    std::array<double,9> get_Jacobian( const int &quaindex ) const override;

    std::array<double,9> get_invJacobian( const int &quaindex ) const override;

    double get_detJac(const int &quaindex) const override {return detJac[quaindex];}

    // Build basis and build the boundary element
    //   Tet-Face-0 : Node 1 2 3 5 9 8
    //   Tet-Face-1 : Node 0 3 2 7 9 6
    //   Tet-Face-2 : Node 0 1 3 4 8 7
    //   Tet-Face-3 : Node 0 2 1 6 5 4
    void buildBasis( const int &face_id, const IQuadPts * const &quad_rule_s,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) override;

    // Get the outwardnormal on faces
    Vector_3 get_2d_normal_out( const int &quaindex, double &area ) const override
    {return triangle_face->get_2d_normal_out( quaindex, area );}

  private:
    static constexpr int nLocBas = 10;

    const int numQuapts;

    // R: 0 <= ii < 10 numQuapts
    std::vector<double> R {}, dR_dx {}, dR_dy {}, dR_dz {},
        d2R_dxx {}, d2R_dyy {}, d2R_dzz {}, d2R_dxy {}, d2R_dxz {}, d2R_dyz {};

    // Container for
    // dx_dr : 0 <= ii < 9 numQuapts
    std::vector<double> dx_dr {};

    // dr_dx : 0 <= ii < 9 numQuapts
    std::vector<double> dr_dx {};

    // detJac : 0 <= ii < numQuapts
    std::vector<double> detJac {};

    std::unique_ptr<FEAElement> triangle_face;
};

#endif
