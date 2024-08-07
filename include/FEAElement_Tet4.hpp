#ifndef FEAELEMENT_TET4_HPP
#define FEAELEMENT_TET4_HPP
// ==================================================================
// FEAElement_Tet4.hpp
// Element routine for the linear tetrahedral element, with evaluation
// of shape functions and their derivatives.
// 
// Tet4 means 4-node tet, aka linear tets.
//
// Date Created: Jan 19 2017
// ==================================================================
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FE_Tools.hpp"

class FEAElement_Tet4 : public FEAElement
{
  public:
    FEAElement_Tet4( const int &in_nqua );

    virtual ~FEAElement_Tet4();

    virtual int get_elemDim() const {return 3;}

    // A unique number for this element.
    virtual int get_Type() const {return 501;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 4;}

    virtual void print_info() const;

    virtual double get_memory_usage() const;

    // Given the quadrature points and nodal coordinates, evaluate the basis 
    // functions and their derivatives up to second order
    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );

    // Return the element size.
    // For the linear tet element, we calculate the DIAMETER of the
    // circumscribing sphere
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

    virtual std::vector<double> get_d2R_dxx( const int &quaindex ) const
    {
      ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_d2R_dxx function error.\n" );
      return { 0.0, 0.0, 0.0, 0.0 };
    }

    virtual std::vector<double> get_d2R_dyy( const int &quaindex ) const
    {
      ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_d2R_dyy function error.\n" );
      return { 0.0, 0.0, 0.0, 0.0 };
    }

    virtual std::vector<double> get_d2R_dzz( const int &quaindex ) const
    {
      ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_d2R_dzz function error.\n"  );
      return { 0.0, 0.0, 0.0, 0.0 };
    }

    virtual std::vector<double> get_d2R_dxy( const int &quaindex ) const
    {
      ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_d2R_dxy function error.\n"  );
      return { 0.0, 0.0, 0.0, 0.0 };
    }

    virtual std::vector<double> get_d2R_dxz( const int &quaindex ) const
    {
      ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_d2R_dxz function error.\n"  );
      return { 0.0, 0.0, 0.0, 0.0 };
    }

    virtual std::vector<double> get_d2R_dyz( const int &quaindex ) const
    {
      ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_d2R_dyz function error.\n"  );
      return { 0.0, 0.0, 0.0, 0.0 };
    }

    // Get the Jacobian matrix dx/dr
    virtual void get_Jacobian(const int &quaindex, double * const &jac_value) const;

    virtual std::array<double,9> get_Jacobian( const int &quaindex ) const
    {
      ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_Jacobian function error.\n"  );
      return {{Jac[0], Jac[1], Jac[2], Jac[3], Jac[4], Jac[5], Jac[6], Jac[7], Jac[8]}};
    }

    // Get the inverse Jacobian matrix dr/dx
    virtual void get_invJacobian(const int &quaindex, double * const &jac_value) const;

    virtual std::array<double,9> get_invJacobian( const int &quaindex ) const
    {
      ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_invJacobian function error.\n"  );
      return {{Jac[9], Jac[10], Jac[11], Jac[12], Jac[13], Jac[14], Jac[15], Jac[16], Jac[17]}};
    }

    // Get the determinant of the Jacobian matrix
    virtual double get_detJac(const int &quaindex) const {return detJac;}

    // Build basis and build the boundary element
    //   Tet-Face-0 : Node 1 2 3
    //   Tet-Face-1 : Node 0 3 2
    //   Tet-Face-2 : Node 0 1 3
    //   Tet-Face-3 : Node 0 2 1
    virtual void buildBasis( const int &face_id, const IQuadPts * const &quad_rule_s,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );

    // Get the outwardnormal on faces
    // This function requires the buildBasis with face_id provided, so that the
    // proper triangle face element is constructed
    // The node numbering of the face element guarantees the get_2d_normal_out
    // returns the outward normal vector
    // See FE_T::QuadPts_on_face function for more details.
    virtual Vector_3 get_2d_normal_out( const int &quaindex, double &area ) const
    {return triangle_face->get_2d_normal_out( quaindex, area );}

    virtual void get_face_ctrlPts( const int &face_id,
        const double * const &volctrl_x,
        const double * const &volctrl_y,
        const double * const &volctrl_z,
        std::vector<double> &facectrl_x,
        std::vector<double> &facectrl_y,
        std::vector<double> &facectrl_z );

    virtual std::array<std::vector<double>, 3> get_face_ctrlPts( const int &face_id,
        const double * const &volctrl_x,
        const double * const &volctrl_y,
        const double * const &volctrl_z );

  private:
    // Number of quadrature points
    const int numQuapts;

    // R : 0 <= ii < 4 x numQuapts
    double * R;

    // tet4 is linear, thus the first-order derivatives are constant
    double dR_dx[4], dR_dy[4], dR_dz[4];

    // Container for
    // dx_dr : 0 <= ii < 9
    // dr_dx : 9 <= ii < 18
    double Jac[18]; 

    double detJac;

    FEAElement * triangle_face;
};

#endif
