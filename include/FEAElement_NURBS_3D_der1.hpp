#ifndef FEAELEMENT_NURBS_3D_DER1_HPP
#define FEAELEMENT_NURBS_3D_DER1_HPP
// ==================================================================
// FEAElement_NURBS_3D_der1.hpp
// This is an implementation of the element routine for 3D NURBS.
// The basis function are evaluated for the function and its 1st order
// derivatives.
//
// The Jacobian and the inverse of the Jacobian w.r.t. s/t/u-coordinates
// are evaluated. In VMS, to evaluate the stabilization parameter, one 
// needs the inverse of the Jacobian w.r.t. xi/eta/zeta in [0,1] domain. 
// Hence, one need to multiply the ds_dxi = 1/hx, dt_deta = 1/hy, 
// du_dzeta = 1/hz:
//
//                 ds_d* by 1.0 / hx, 
//                 dt_d* by 1.0 / hy, 
//                 du_d* by 1.0 / hz. 
// This design is different from my previous NURBS element evaluations.
//
// The design is intended for a no-cache style programming. The objective
// is to reduce the heap usage and reduce the calls of malloc and new
// operators to optimize the use of the memory pool.
//
// A typical usage of this element routine is
// 1. Setup the memory allocation by calling the constructors;
// 2. Call buildBasis function to evaluate the basis functions;
// 3. Call get_xxx functions to obtain the basis function values as 
//    well as the values of Jacobians, normal vectors, etc.
// 4. In case of the user needs to change the polynomial degree or the
//    quadrature points, call reset_xxx functions to re-allocate the
//    container's size.
// 
// Author: Ju Liu
// Date: June 25 2016
// ==================================================================

#include "FEAElement.hpp"
#include "Matrix_double_6by6_Array.hpp"

class FEAElement_NURBS_3D_der1 : public FEAElement
{
  public:
    FEAElement_NURBS_3D_der1( const int &in_sdeg, const int &in_tdeg, 
        const int &in_udeg, const int &in_nquas, const int &in_nquat, 
        const int &in_nquau );

    virtual ~FEAElement_NURBS_3D_der1();

    virtual int get_elemDim() const {return 3;}

    virtual int get_Type() const {return 131;}

    virtual int get_numType() const {return 1;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual void reset_degree( const int &new_sdeg, const int &new_tdeg,
        const int &new_udeg );

    virtual void reset_numQua( const int &new_squa, const int &new_tqua,
        const int &new_uqua );

    virtual void print() const;

    virtual double get_memory_usage() const;

    virtual void buildBasis( const double &hx, 
        const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z,
        const double * const &ctrl_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z );

    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual double get_detJac(const int &quaindex) const;

    virtual void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y,
        double * const &basis_z ) const;

    virtual void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y, double * const &basis_z ) const;

    virtual void get_Jacobian(const int &quaindex, double * const &jac_value) const;

    virtual void get_invJacobian(const int &quaindex, double * const &jac_value) const;

    virtual void get_3d_normal_bottom( const int &quaindex,
        double &nx, double &ny, double &nz, double &surfaceArea ) const;

    virtual void get_3d_normal_top( const int &quaindex,
        double &nx, double &ny, double &nz, double &surfaceArea ) const;

    virtual void get_3d_normal_left( const int &quaindex,
        double &nx, double &ny, double &nz, double &surfaceArea ) const;

    virtual void get_3d_normal_right( const int &quaindex,
        double &nx, double &ny, double &nz, double &surfaceArea ) const;

    virtual void get_3d_normal_front( const int &quaindex,
        double &nx, double &ny, double &nz, double &surfaceArea ) const;

    virtual void get_3d_normal_back( const int &quaindex,
        double &nx, double &ny, double &nz, double &surfaceArea ) const;

  private:
    int num_qua_s, num_qua_t, num_qua_u, numQuapts;
    int sdeg, tdeg, udeg, sdp1, tdp1, udp1;
    int nLocBas, rlength;

    // Container for
    // R     :           0 <= ii < rlength
    // dR_dx :     rlength <= ii < 2 * rlength
    // dR_dy : 2 * rlength <= ii < 3* rlength
    // dR_dz : 3 * rlength <= ii < 4 * rlength
    double * R;

    // Container for
    // dx_ds  :              0 <= ii < 9 * numQuapts
    // ds_dx  :  9 * numQuapts <= ii < 18 * numQuapts
    // detJac : 18 * numQuapts <= ii < 19 * numQuapts
    double * Jac;


    // Bezier element
    // dR_ds   :           0 <= ii < nLocBas
    // dR_dt   :     nLocBas <= ii < 2 * nLocBas
    // dR_du   : 2 * nLocBas <= ii < 3 * nLocBas
    double * dRr;


    // univariate functions
    // Nns       :        0 <= ii < sdp1
    // dNns_ds   :     sdp1 <= ii < 2 * sdp1
    // Nnt       :        0 <= ii < tdp1
    // dNnt_dt   :     tdp1 <= ii < 2 * tdp1
    // Nnu       :        0 <= ii < udp1
    // dNnu_du   :     udp1 <= ii < 2 * udp1
    double * Nns, * Nnt, * Nnu;

    virtual void clearBasisCache();

    void BuildShape_atQua( const int &quaindex,
        const double &hx, const double &hy, const double &hz,
        const double * const &ctrl_x, const double * const &ctrl_y,
        const double * const &ctrl_z, const double * const &ctrl_w );


    void resize_container();

    void normalize_3d_vector(double &x, double &y, double &z, double &len) const
    {
      len = sqrt(x*x + y*y + z*z);
      x = x/len;
      y = y/len;
      z = z/len;
    }
};


#endif
