#ifndef FEAELEMENT_BS_2D_DER2_HPP
#define FEAELEMENT_BS_2D_DER2_HPP
// ============================================================================
// FEAElement_BS_2D_der2.hpp
// This is an implementation of the element routine for 2D B-splines. The basis
// functions are evaluated for the function and its derivatives up to the second
// order.
//
// The Jacobian and the inverse of the Jacobian with respect to s/t/u are
// evaluated. (Note: If one needs dx_dsprime, hx/hy/hz need to be applied for
// the change of variable from s to s_prime. This procedure is needed in the
// evaluation of the VMS stabilization parameter, for example.)
//
// The design is intended for no-cache-style programming. The objective is to
// reduce the heap usage and reduce the call of malloc and new operators. This
// element routine is also intended for fast evaluation.
// 
// A typical usage of this element routine is
// 1. Setup the memory allocation by calling the constructors;
// 2. Call buildBasis function to evaluate the basis functions;
// 3. Call get_xxx functions to obtain the basis function values as well as the
//    values of Jacobians, normal vectors, etc.
// 4. In case of the user needs to change the polynomial degree or the
//    quadrature points, call reset_xxx functions to re-allocate the container's
//    size.
//
// 
// Date: Sept. 28 2015
// ============================================================================

#include "FEAElement.hpp"

#include "Matrix_double_3by3_Array.hpp"

class FEAElement_BS_2D_der2 : public FEAElement
{
  public:
    FEAElement_BS_2D_der2( const int &in_sdeg, const int &in_tdeg,
        const int &in_nquas, const int &in_nquat );

    virtual ~FEAElement_BS_2D_der2();

    virtual int get_elemDim() const {return 2;}

    virtual int get_Type() const {return 622;}

    virtual int get_numType() const {return 1;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual void reset_degree(const int &new_sdeg, const int &new_tdeg);

    virtual void reset_numQua( const int &new_squa, const int &new_tqua );

    virtual void reset_numQua( const int &new_squa, const int &new_tqua,
        const int &new_uqua );

    virtual void print() const;

    virtual double get_memory_usage() const;

    virtual void buildBasis( const double &hx, const double &hy,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ext_x,
        const double * const &ext_y );

    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual double get_detJac(const int &quaindex) const;

    virtual void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y ) const;

    virtual void get_2D_R_dR_d2R( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y,
        double * const &basis_xx, double * const &basis_yy, double * const &basis_xy ) const; 

    virtual void get_2D_R_gradR_LaplacianR( const int &quaindex,
        double * const &basis, double * const &basis_x, double * const &basis_y,
        double * const &basis_xx, double * const &basis_yy ) const;

    virtual void get_Jacobian(const int &quaindex, double * const &jac_value) const;

    virtual void get_invJacobian(const int &quaindex, double * const &jac_value) const;

    virtual void get_2d_normal_back( const int &quaindex,
        double &nx, double &ny, double &line ) const;

    virtual void get_2d_normal_front( const int &quaindex,
        double &nx, double &ny, double &line ) const;

    virtual void get_2d_normal_left( const int &quaindex,
        double &nx, double &ny, double &line ) const;

    virtual void get_2d_normal_right( const int &quaindex,
        double &nx, double &ny, double &line ) const;

  private:
    int num_qua_s, num_qua_t, numQuapts;
    int sdeg, tdeg, sdp1, tdp1;
    int nLocBas;
    int rlength;

    // container for R, dR_dx, dR_dy, d2R_dxx, d2R_dyy, d2R_dxy
    // R       :           0 <= ii < rlength
    // dR_dx   :     rlength <= ii < 2 * rlength
    // dR_dy   : 2 * rlength <= ii < 3 * rlength
    // d2R_dxx : 3 * rlength <= ii < 4* rlength
    // d2R_dyy : 4 * rlength <= ii < 5 * rlength
    // d2R_dxy : 5 * rlength <= ii < 6 * rlength
    double * R;

    // container for dx_ds, ds_dx, and detJac
    // dx_ds  :             0 <= ii < 4 * numQuapts
    // ds_dx  : 4 * numQuapts <= ii < 8 * numQuapts
    // detJac : 8 * numQuapts <= ii < 9 * numQuapts
    double * Jac;

    // container for the Bezier element
    // dR_ds :           0 <= ii < nLocBas
    // dR_dt :     nLocBas <= ii < 2*nLocBas
    // d2R_dss : 2*nLocBas <= ii < 3*nLocBas
    // d2R_dtt : 3*nLocBas <= ii < 4*nLocBas
    // d2R_dst : 4*nLocBas <= ii < 5*nLocBas
    double * dRr;

    // coutainer for univariate component
    // Nns       :      0 <= ii < sdp1
    // dNns_ds   :   sdp1 <= ii < 2*sdp1
    // d2Nns_dss : 2*sdp1 <= ii < 3*sdp1
    // Nnt       :      0 <= ii < tdp1
    // dNnt_dt   :   tdp1 <= ii < 2*tdp1
    // d2Nnt_dtt : 2*tdp1 <= ii < 3*tdp1
    double * Nns, * Nnt;

    virtual void clearBasisCache();

    void BuildShape_atQua( const int &quaindex,
        const double &hx, const double &hy,
        const double * const &ctrl_x,
        const double * const &ctrl_y );

    void resize_container();

    void normalize_2d_vector(double &x, double &y) const
    {
      const double temp = 1.0 / sqrt(x*x+y*y);
      x = x * temp;
      y = y * temp;
    }
};


#endif
