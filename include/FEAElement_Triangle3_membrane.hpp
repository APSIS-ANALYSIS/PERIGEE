#ifndef FEAELEMENT_TRIANGLE3_MEMBRANE_HPP
#define FEAELEMENT_TRIANGLE3_MEMBRANE_HPP
// ==================================================================
// FEAElement_Triangle3_membrane.hpp
// Element routine for the linear triangular element in
// three-dimensional space for the coupled momentum method. 
// Evaluates the element Jacobian, basis functions, and their gradients
// with respect to the local lamina coordinates. Also computes the
// global-to-local rotation matrix Q.
// 
// Triangle3 means 3-node triangle; _membrane indicates a membrane
// assumption such that the displacement is only a function of the
// in-plane parametric coordinates.
//
// Date Created: Jan. 8 2021
// ==================================================================
#include "FEAElement.hpp"

class FEAElement_Triangle3_membrane : public FEAElement
{
  public:
    FEAElement_Triangle3_membrane( const int &in_nqua );

    virtual ~FEAElement_Triangle3_membrane();

    virtual int get_elemDim() const {return 2;}

    // element type : 523
    // 5: simplicial element
    // 2: 2D element
    virtual int get_Type() const {return 523;}

    virtual int get_numType() const {return 1;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual void print() const;

    virtual double get_memory_usage() const;

    // --------------------------------------------------------------
    // Input: 
    // \para : quad_rule quadrature points
    // \para ctrl_x/y/z : the control points' coordinates in the global
    //                    system.
    // This function will generate the global-to_lamina rotation matrix Q,
    // the basis functions, and the basis functions' gradient with respect
    // to the lamina coorindates at each quadratue point.
    // --------------------------------------------------------------
    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );

    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y, double * const &basis_z ) const;

    virtual void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y,
        double * const &basis_z ) const;

    virtual void get_rotationMatrix( const int &quaindex, Matrix_3x3 &rot_mat ) const;

    // Assumes the triangle nodes are arranged such that the outward
    // direction is given by dx_dr x dx_ds
    virtual void get_2d_normal_out( const int &quaindex,
        double &nx, double &ny, double &nz, double &area ) const;

    // If the triangle nodes are NOT arranged in any particular order,
    // use an interior node to define the outward direction.
    virtual void get_normal_out( const int &quaindex,
        const double &sur_pt_x, const double &sur_pt_y, const double &sur_pt_z,
        const double &intpt_x, const double &intpt_y, const double &intpt_z,
        double &nx, double &ny, double &nz, double &len ) const;

    virtual double get_detJac(const int &quaindex) const
    {return detJac;}

  private:
    const int nLocBas;

    int numQuapts;

    // Containers for rotated *local* coordinates
    double * ctrl_xl, * ctrl_yl, * ctrl_zl;

    // container for R0 = 1 - r - s, R1 = r, R2 = s :
    // 0 <= ii < 3 x numQuapts
    double * R;

    // containers for gradients of basis functions with respect to
    // rotated *local* coordinates. Derivatives are constant here.
    double * dR_dx, * dR_dy;

    // containers for unit vectors used to construct rotation matrix Q,
    // each of length 3 
    double e_r[3], e_s[3], e_l1[3], e_l2[3];

    // global-to-local 3x3 rotation matrix
    Matrix_3x3 Q;

    // unit outward normal vector
    double unx, uny, unz;

    // Container for rotated *local* 2D Jacobian and its inverse
    // dx_dr : 0 <= ii < 4
    // dr_dx : 4 <= ii < 8
    double Jac[8];

    // Jacobian determinant 
    double detJac;

    // deallocate the memory
    virtual void clearBasisCache();
};

#endif
