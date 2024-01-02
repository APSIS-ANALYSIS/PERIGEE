#ifndef FEAELEMENT_TET10_V2_HPP
#define FEAELEMENT_TET10_V2_HPP
// ==================================================================
// FEAElement_Tet10_v2.hpp
// Element routine for quadratic 10-node tetrahedral element.
// In this version v2, the node numbering is made compatible
// compatible with the vtk format -- see the graph below.
// 
// Notice that only nodes 8 and 9 are switched. The remaining nodes
// are numbered identically to those in FEAElement_Tet10 class.
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
#include "FEAElement.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FE_Tools.hpp"

class FEAElement_Tet10_v2 : public FEAElement
{
  public:
    FEAElement_Tet10_v2( const int &in_nqua );

    virtual ~FEAElement_Tet10_v2();

    virtual int get_elemDim() const {return 3;}

    // A unique number for this element. 
    virtual int get_Type() const {return 502;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 10;}

    virtual void print_info() const;

    virtual double get_memory_usage() const;

    // Given the quadrature points and nodal coordinates, evaluate
    // the basis functions and their derivatives up to second order
    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );

    // Return the element size.
    // Here we adopt the algorithm for Tet4 and use the four vertex
    // nodes to calculate the element size
    virtual double get_h( const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const;

    // get_xxx functions give access to function evaluations at the
    // quadrature point corresponding to quaindex
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

    // Build basis and build the boundary element
    //   Tet-Face-0 : Node 1 2 3 5 9 8
    //   Tet-Face-1 : Node 0 3 2 7 9 6
    //   Tet-Face-2 : Node 0 1 3 4 8 7
    //   Tet-Face-3 : Node 0 2 1 6 5 4
    virtual void buildBasis( const IQuadPts * const &quad_rule_s, const int &face_id,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );

    // Get the outwardnormal on faces
    virtual Vector_3 get_2d_normal_out( const int &quaindex, double &area ) const
    {return triangle_face->get_2d_normal_out( quaindex, area );}

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

    FEAElement * triangle_face;
};

#endif
