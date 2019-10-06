#ifndef FEAELEMENT_NURBS_2D_HPP
#define FEAELEMENT_NURBS_2D_HPP
// ==================================================================
// FEAElement_NURBS_2D.hpp
// This is an implementation of the element routine for 2D NURBS. The
// basis functions are evaluated for the functions, their 1st order 
// and 2nd order derivatives. The Jacobian and its inverse w.r.t. s/
// t/u-coordinates are also evaluated.
//
// Author: Ju Liu
// Date: Sept. 2 2016
// ==================================================================

#include "FEAElement.hpp"
#include "Matrix_double_3by3_Array.hpp"

class FEAElement_NURBS_2D : public FEAElement
{
  public:
    FEAElement_NURBS_2D( const int &in_sdeg, const int &in_tdeg,
        const int &in_nquas, const int &in_nquat, const bool &in_is2ndder );

    virtual ~FEAElement_NURBS_2D();

    virtual int get_elemDim() const {return 2;}

    virtual int get_Type() const {return 122;}

    virtual int get_numType() const {return 1;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual void reset_degree( const int &new_sdeg, const int &new_tdeg );

    virtual void reset_numQua( const int &new_squa, const int &new_tqua );

    virtual void reset_numQua( const int &new_squa, const int &new_tqua,
        const int &new_uqua ) {reset_numQua(new_squa, new_tqua);}

    virtual void print() const;

    virtual double get_memory_usage() const;

    virtual void buildBasis( const double &hx, const double &hy,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_w,
        const double * const &ext_x,
        const double * const &ext_y );

    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual double get_detJac(const int &quaindex) const;

    virtual void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y ) const;

    virtual void get_2D_R_dR_d2R( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y,
        double * const &basis_xx, double * const &basis_yy, 
        double * const &basis_xy ) const;

    virtual void get_2D_R_gradR_LaplacianR( const int &quaindex,
        double * const &basis, double * const &basis_x, double * const &basis_y,
        double * const &basis_xx, double * const &basis_yy ) const;

    virtual void get_Jacobian(const int &quaindex, 
        double * const &jac_value) const;

    virtual void get_invJacobian(const int &quaindex, 
        double * const &jac_value) const;

    virtual void get_2d_normal_left( const int &quaindex,
        double &nx, double &ny, double &line ) const;

    virtual void get_2d_normal_right( const int &quaindex,
        double &nx, double &ny, double &line ) const;

    virtual void get_2d_normal_back( const int &quaindex,
        double &nx, double &ny, double &line ) const;

    virtual void get_2d_normal_front( const int &quaindex,
        double &nx, double &ny, double &line ) const;

  private:
    int num_qua_s, num_qua_t, numQuapts;
    int sdeg, tdeg, sdp1, tdp1;
    int nLocBas, rlength;
    bool is2ndder;

    // Container for
    // R     :             0 <= ii < rlength
    // dR_dx :       rlength <= ii < 2 * rlength
    // dR_dy :   2 * rlength <= ii < 3 * rlength
    // if is2ndder == true
    // d2R_dxx : 3 * rlength <= ii < 4 * rlength
    // d2R_dyy : 4 * rlength <= ii < 5 * rlength
    // d2R_dxy : 5 * rlength <= ii < 6 * rlength
    double * R;

    // Container for
    // dR_ds : 0 <= ii < nLocBas
    // dR_dt : nLocBas <= ii < 2 * nLocBas
    // if is2ndder == true
    // d2R_dss : 2 * nLocBas <= ii < 3 * nLocBas
    // d2R_dtt : 3 * nLocBas <= ii < 4 * nLocBas
    // d2R_dst : 4 * nLocBas <= ii < 5 * nLocBas
    double * dRr;

    // Container for
    // dx_ds  :             0 <= ii < 4 * numQuapts
    // ds_dx  : 4 * numQuapts <= ii < 8 * numQuapts
    // detJac : 8 * numQuapts <= ii < 9 * numQuapts
    double * Jac;

    // Container for univariate component
    // Nns :              0 <= ii < sdp1
    // dNns_ds :       sdp1 <= ii < 2 * sdp1
    // if is2ndder == true
    // d2Nns_dss : 2 * sdp1 <= ii < 3 * sdp1
    double * Nns;

    // Container for univariate component
    // Nnt :              0 <= ii < tdp1
    // dNnt_dt :       tdp1 <= ii < 2 * tdp1
    // if is2ndder == true
    // d2Nnt_dtt : 2 * tdp1 <= ii < 3 * tdp1
    double * Nnt;

    virtual void clearBasisCache();

    // Build only up to 1st order derivatives
    void BuildShape_atQua1( const int &quaindex, const double &hx,
        const double &hy, const double * const &ctrl_x,
        const double * const &ctrl_y, 
        const double * const &ctrl_w );

    // Build up to 2nd order derivatives
    void BuildShape_atQua2( const int &quaindex, const double &hx,
        const double &hy, const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_w );

    void resize_container();

    void normalize_2d_vector(double &x, double &y) const
    {
      const double temp = 1.0 / sqrt(x*x+y*y);
      x = x * temp;
      y = y * temp;
    }

};


#endif
