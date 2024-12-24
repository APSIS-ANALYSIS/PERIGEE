#ifndef FEAELEMENT_TET4_HPP
#define FEAELEMENT_TET4_HPP
// ============================================================================
// FEAElement_Tet4.hpp
// Element routine for the linear tetrahedral element, with evaluation of shape
// functions and their derivatives.
// 
// Tet4 means 4-node tet, aka linear tets.
//
// Date Created: Jan 19 2017
// ============================================================================
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FE_Tools.hpp"

class FEAElement_Tet4 final : public FEAElement
{
  public:
    FEAElement_Tet4( const int &in_nqua );

    ~FEAElement_Tet4() override = default;

    int get_elemDim() const override {return 3;}

    FEType get_Type() const override {return FEType::Tet4;}

    int get_numQuapts() const override {return numQuapts;}

    int get_nLocBas() const override {return 4;}

    void print_info() const override;

    // Given the quadrature points and nodal coordinates, evaluate the basis 
    // functions and their derivatives up to second order
    void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) override;

    // Return the element size.
    // For the linear tet element, we calculate the DIAMETER of the
    // circumscribing sphere
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

    std::array<double,9> get_Jacobian( const int &quaindex ) const override
    {
      ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_Jacobian function error.\n"  );
      return {{Jac[0], Jac[1], Jac[2], Jac[3], Jac[4], Jac[5], Jac[6], Jac[7], Jac[8]}};
    }

    // Get the inverse Jacobian matrix dr/dx
    void get_invJacobian(const int &quaindex, double * const &jac_value) const override;

    std::array<double,9> get_invJacobian( const int &quaindex ) const override
    {
      ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_invJacobian function error.\n"  );
      return {{Jac[9], Jac[10], Jac[11], Jac[12], Jac[13], Jac[14], Jac[15], Jac[16], Jac[17]}};
    }

    // Get the determinant of the Jacobian matrix
    double get_detJac(const int &quaindex) const override {return detJac;}

    // Build basis and build the boundary element
    //   Tet-Face-0 : Node 1 2 3
    //   Tet-Face-1 : Node 0 3 2
    //   Tet-Face-2 : Node 0 1 3
    //   Tet-Face-3 : Node 0 2 1
    void buildBasis( const int &face_id, const IQuadPts * const &quad_rule_s,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) override;

    // Get the outwardnormal on faces
    // This function requires the buildBasis with face_id provided, so that the
    // proper triangle face element is constructed
    // The node numbering of the face element guarantees the get_2d_normal_out
    // returns the outward normal vector
    // See FE_T::QuadPts_on_face function for more details.
    Vector_3 get_2d_normal_out( const int &quaindex, double &area ) const override
    {return triangle_face->get_2d_normal_out( quaindex, area );}

    std::array<std::vector<double>, 3> get_face_ctrlPts( const int &face_id,
        const double * const &volctrl_x,
        const double * const &volctrl_y,
        const double * const &volctrl_z ) const override;

  private:
    // Number of quadrature points
    const int numQuapts;

    // R : 0 <= ii < 4 x numQuapts
    std::vector<double> R {};

    // tet4 is linear, thus the first-order derivatives are constant
    std::array<double,4> dR_dx, dR_dy, dR_dz;

    // Container for
    // dx_dr : 0 <= ii < 9
    // dr_dx : 9 <= ii < 18
    std::array<double,18> Jac; 

    double detJac;

    std::unique_ptr<FEAElement> triangle_face;
};

#endif
